import logging
import pickle
from bisect import bisect_right
from itertools import accumulate
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
import torch
from filelock import FileLock
from qqtools import qData, qDictDataset
from qqtools.qm import ensureAtomNumber

# eV per atom, fitted without an intercept on all v4 Main records.
REF_ENERGIES = {
    1: -16.317913385676807,
    6: -1036.239848400066,
    7: -1488.6023276581857,
    8: -2045.266547983185,
    9: -2714.676260744514,
    15: -9288.809490307502,
    16: -10831.817895018525,
    17: -12518.89928730063,
    35: -70038.33428031162,
    53: -8103.095336224553,
}

ENERGY_FIELDS = ("offsets", "z", "pos", "energy")
CACHE_FIELDS = (*ENERGY_FIELDS, "ref_energy", "force")


def pack_chunk(records, first_index):
    """Store contiguous atom arrays with molecule offsets for memory-mapped reads."""
    atoms, positions, forces = [], [], []
    energies, references = [], []
    offsets = [0]
    for index, record in enumerate(records, first_index):
        elements = record["elements"]
        if not isinstance(elements, (list, tuple, np.ndarray)):
            elements = [elements]
        # qqtools 1.3.28 accepts Python ints but rejects NumPy integer scalars.
        z = ensureAtomNumber([int(el) if isinstance(el, np.integer) else el for el in elements])
        pos = np.asarray(record.get("coordinates", record.get("position")), dtype=np.float64)
        raw_force = record.get("forces", record.get("force"))
        force = np.asarray(raw_force, dtype=np.float64) if raw_force is not None else np.zeros_like(pos)
        if not z or pos.shape != (len(z), 3) or force.shape != pos.shape:
            raise ValueError(f"Invalid molecular geometry at record {index}")
        reference = sum(REF_ENERGIES[number] for number in z)
        atoms.extend(z)
        positions.append(pos)
        forces.append(force)
        energies.append(float(record["energy"]) - reference)
        references.append(reference)
        offsets.append(len(atoms))

    return dict(
        z=torch.tensor(atoms, dtype=torch.long),
        pos=torch.from_numpy(np.concatenate(positions)),
        force=torch.from_numpy(np.concatenate(forces)),
        energy=torch.tensor(energies, dtype=torch.float64),
        ref_energy=torch.tensor(references, dtype=torch.float64),
        offsets=torch.tensor(offsets),
    )


def legacy_cache_path(cache):
    return cache.with_suffix("").with_suffix(".tensors.pt")


def ensure_chunk_cache(raw, cache, first_index, force_reload):
    if cache.exists() and not force_reload:
        return
    # Coordinate all users of this shared cache, including independent jobs.
    with FileLock(str(cache) + ".lock", timeout=-1):
        if cache.exists() and not force_reload:
            return
        legacy = legacy_cache_path(cache)
        if legacy.exists() and not force_reload:
            old = torch.load(legacy, mmap=True, weights_only=True, map_location="cpu")
            chunk = {key: old[key] for key in CACHE_FIELDS}
            del old
        else:
            with Path(raw).open("rb") as handle:
                chunk = pack_chunk(pickle.load(handle), first_index)
        with NamedTemporaryFile(dir=cache.parent, delete=False) as tmp:
            temp_path = Path(tmp.name)
        try:
            torch.save(chunk, temp_path)
            temp_path.replace(cache)
        finally:
            temp_path.unlink(missing_ok=True)


class QO2MolDataset(qDictDataset):
    main_file_names = [f"QO2Mol26_Main_chunk{i}.pkl" for i in range(10)]
    ood_file_names = ["QO2Mol26_OOD_medium.pkl", "QO2Mol26_OOD_large.pkl"]
    atom_numbers = list(REF_ENERGIES)

    @property
    def raw_file_names(self):
        scopes = {"main": self.main_file_names, "ood_medium": self.ood_file_names[:1],
                  "ood_large": self.ood_file_names[1:]}
        if self.scope not in scopes:
            raise ValueError(f"Unknown scope {self.scope!r}; expected one of {tuple(scopes)}")
        return scopes[self.scope]

    @property
    def processed_file_names(self):
        return [Path(name).stem + ".v4-ref-ev.core.pt" for name in self.raw_file_names]

    def __init__(self, root, scope="main", force_reload=False, double=False, *, energy_only=False):
        super().__init__()
        self.root = Path(root)
        self.scope = scope
        self.energy_only = energy_only
        self.dtype = torch.float64 if double else torch.float32
        self.processed_dir.mkdir(parents=True, exist_ok=True)
        self.chunks = []
        self._identifiers = {}
        first_index = 0
        paths = [(Path(raw), Path(cache)) for raw, cache in zip(self.raw_paths, self.processed_paths)]
        for raw, cache in paths:
            if not raw.is_file() and (force_reload or not (cache.is_file() or legacy_cache_path(cache).is_file())):
                raise FileNotFoundError(f"Missing {self.scope} input: {raw}")
        for number, (raw, cache) in enumerate(paths, 1):
            action = "Reusing" if cache.is_file() and not force_reload else "Preparing"
            logging.getLogger(__name__).info("[%d/%d] %s %s", number, len(paths), action, cache.name)
            ensure_chunk_cache(raw, cache, first_index, force_reload)
            chunk = torch.load(cache, mmap=True, weights_only=True, map_location="cpu")
            if energy_only:
                chunk = {key: chunk[key] for key in ENERGY_FIELDS}
            self.chunks.append(chunk)
            first_index += len(chunk["energy"])
        self.ends = list(accumulate(len(chunk["energy"]) for chunk in self.chunks))

    def len(self):
        return self.ends[-1]

    def __iter__(self):
        for true_index in self.indices():
            yield self.get(true_index)

    def _chunk_identifiers(self, chunk_index):
        if chunk_index not in self._identifiers:
            raw = Path(self.raw_paths[chunk_index])
            if raw.exists():
                start = self.ends[chunk_index - 1] if chunk_index else 0
                with raw.open("rb") as handle:
                    records = pickle.load(handle)
                identifiers = [(str(r.get("confid", i)), str(r.get("inchikey", "")))
                               for i, r in enumerate(records, start)]
            else:
                legacy = legacy_cache_path(Path(self.processed_paths[chunk_index]))
                if not legacy.exists():
                    raise FileNotFoundError(f"Full sample identifiers require {raw}")
                chunk = torch.load(legacy, mmap=True, weights_only=True, map_location="cpu")
                identifiers = list(zip(chunk["name"], chunk["inchi_key"]))
            self._identifiers[chunk_index] = identifiers
        return self._identifiers[chunk_index]

    def get(self, index):
        index = range(self.len())[index]
        chunk_index = bisect_right(self.ends, index)
        chunk = self.chunks[chunk_index]
        local = index - (self.ends[chunk_index - 1] if chunk_index else 0)
        start, end = chunk["offsets"][local:local + 2].tolist()
        sample = qData(num_nodes=end - start, z=chunk["z"][start:end],
                       pos=chunk["pos"][start:end].to(self.dtype), energy=chunk["energy"][local].to(self.dtype))
        if not self.energy_only:
            name, inchi_key = self._chunk_identifiers(chunk_index)[local]
            sample.update(atom_count=end - start, force=chunk["force"][start:end].to(self.dtype),
                          ref_energy=chunk["ref_energy"][local].to(self.dtype),
                          name=name, inchi_key=inchi_key)
        return sample

    def get_norm_factor(self):
        energy = torch.cat([chunk["energy"] for chunk in self.chunks]).to(self.dtype).numpy()
        indices = self.indices()
        if isinstance(indices, range) and indices.step == 1:
            energy = energy[indices.start:indices.stop]
        else:
            energy = energy[np.asarray(indices, dtype=np.intp)]
        return float(energy.mean()), float(energy.std())

    def get_elements(self):
        return sorted(set().union(*(torch.unique(chunk["z"]).tolist() for chunk in self.chunks)))
