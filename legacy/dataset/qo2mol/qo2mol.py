import sys
from itertools import chain
from pathlib import Path

import torch
from qqtools import Timer, load_pickle
from torch_geometric.data import InMemoryDataset

proj_root = Path(__file__).resolve().parent.parent.parent.parent
if __name__ == "__main__":
    sys.path.insert(0, str(proj_root))

from legacy.dataset.qo2mol.datasetbase import BaseQMDataset
from legacy.dataset.qo2mol.utils import convert_std_input

mp_atomNumber2Idx = {1: 0, 6: 1, 7: 2, 8: 3, 9: 4, 15: 5, 16: 6, 17: 7, 35: 8, 53: 9}
REF_ENERGIES = [
    -377.1312596544013,
    -23896.71978240409,
    -34329.56617414761,
    -47167.507219805506,
    -62605.90249199728,
    -214201.9540501211,
    -249791.03533796355,
    -288696.57944521913,
    -1615128.7577386901,
    -186866.01419930425,
]


class QO2MolDataset(BaseQMDataset):
    """
    The QO2Mol Dataset
    """

    atom_numbers = [1, 6, 7, 8, 9, 15, 16, 17, 35, 53]
    main_file_names = [f"QO2Mol_main_chunk_{i}.pkl" for i in range(10)]
    other_file_names = ["QO2Mol_set_b.pkl", "QO2Mol_set_c.pkl"]

    def __init__(self, root, scope="main", force_reload=False, double=False):
        with Timer("Initializing dataset...", prefix="[Dataset] "):
            self.scope = scope
            self.force_reload = force_reload
            super().__init__(root)

            # The locally processed legacy cache contains PyG Data objects.
            self.data_dict = torch.load(self.processed_dir / self.data_filepath, weights_only=False)
            self._data, self.slices = self.data_dict["data"], self.data_dict["slices"]
            self.check_double(double=double)

    @property
    def raw_file_names(self):
        return self.main_file_names + self.other_file_names[:1]

    @property
    def processed_file_names(self):
        return ["QO2Mol_main_processed.pt", "QO2Mol_b_processed.pt"]

    @property
    def data_filepath(self):
        if self.scope == "main":
            return self.processed_file_names[0]
        elif self.scope == "b":
            return self.processed_file_names[1]

    def _process(self):

        # cache check
        if all(Path(fp).exists() for fp in self.processed_paths):
            return

        main_file_paths = [self.raw_dir / fn for fn in self.main_file_names]
        set_b_file = self.raw_dir / self.other_file_names[0]
        for path in [*main_file_paths, set_b_file]:
            if not path.is_file():
                raise FileNotFoundError(path)
        self.processed_dir.mkdir(parents=True, exist_ok=True)

        self._process_files(main_file_paths, self.processed_paths[0], REF_ENERGIES)
        self._process_files([set_b_file], self.processed_paths[1], REF_ENERGIES)

    def _process_files(self, file_paths, processed_fp, reference_energies):
        with Timer("Reading and converting molecular records..."):
            records = chain.from_iterable(map(load_pickle, file_paths))
            data_list = convert_std_input(records, mp_atomNumber2Idx, reference_energies)

        with Timer("Collating graphs..."):
            data, slices = InMemoryDataset.collate(data_list)
            del data_list

        data_dict = {"data": data, "slices": slices}

        with Timer(f"Saving {processed_fp}..."):
            torch.save(data_dict, processed_fp)


if __name__ == "__main__":

    root_dir = proj_root / "download"
    print(f"data root: {root_dir}")
    dataset = QO2MolDataset(root=root_dir)
    print(dataset)
