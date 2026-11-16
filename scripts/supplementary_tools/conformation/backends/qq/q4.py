"""Seeded combination sampling and coordinate-based, rotation-invariant deduplication."""

import copy
from dataclasses import asdict
from rdkit.Chem import rdMolTransforms
import math
import random
import numpy as np
from ...options import Q4Options
from .geometry import (
    positions,
    moving_component,
    ContactCheck,
    aligned_rmsd,
    rdkit_reference,
    set_torsions,
)


def normalize_minima(molecule, minima):
    if not isinstance(minima, dict) or not minima:
        raise ValueError("Q4 requires a nonempty torsion-to-minima mapping")
    normalized, centers = {}, set()
    for term, values in minima.items():
        try:
            atoms = tuple(int(s) for s in term.split("-"))
        except (AttributeError, ValueError) as exc:
            raise ValueError("Expected a four-index torsion key") from exc
        if (
            len(atoms) != 4
            or len(set(atoms)) != 4
            or any(i < 0 or i >= len(molecule.Atoms) for i in atoms)
        ):
            raise ValueError("Q4 requires four distinct valid atom indices")
        if any(j not in molecule.Atoms[i].connect for i, j in zip(atoms, atoms[1:])):
            raise ValueError("Q4 torsion atoms must form a bonded path")
        b, c = atoms[1:3]
        if molecule.Atoms[b].bond_type[molecule.Atoms[b].connect.index(c)] != "1":
            raise ValueError("qq Q4 only rotates single central bonds")
        if moving_component(molecule, b, c) is None:
            raise ValueError("Q4 central bond is in a ring")
        center = tuple(sorted((b, c)))
        if center in centers:
            raise ValueError("Multiple torsion definitions for the same central bond")
        centers.add(center)
        if not isinstance(values, list):
            raise ValueError("Each minimum-angle collection must be a list")
        raw_angles = []
        for angle in values:
            if (
                isinstance(angle, bool)
                or not isinstance(angle, (int, float))
                or not math.isfinite(angle)
            ):
                raise ValueError("Minimum angles must be finite numbers in degrees")
            value = (float(angle) + 180.0) % 360.0 - 180.0
            raw_angles.append(value)
        angles = []
        for value in sorted(raw_angles):
            if not any(
                abs((value - old + 180.0) % 360.0 - 180.0) <= 1e-8 for old in angles
            ):
                angles.append(value)
        canonical = min(atoms, atoms[::-1])
        normalized[canonical] = sorted(angles)
    return dict(sorted(normalized.items()))


def combination_indices(total, limit, seed):
    """Uniform without replacement, O(limit) memory even for enormous products."""
    if total <= limit:
        return list(range(total))
    rng = random.Random(seed)
    selected = set()
    # Floyd's algorithm; unlike random.sample(range(total)), no ssize_t limit.
    for j in range(total - limit, total):
        candidate = rng.randrange(j + 1)
        selected.add(j if candidate in selected else candidate)
    return sorted(selected)


def decode_combination(index, values):
    result = []
    for choices in reversed(values):
        index, remainder = divmod(index, len(choices))
        result.append(choices[remainder])
    return tuple(reversed(result))


def graph_signature(molecule):
    return [
        (a.elem, a.formal_charge, a.isotope, list(a.connect), list(a.bond_type))
        for a in molecule.Atoms
    ]


def generate_q4(molecule, minima, *, limit=64, optimizer=None, options=None):
    if type(limit) is not int or limit < 1:
        raise ValueError("Q4 limit must be a positive integer")
    settings = options or Q4Options()
    normalized = normalize_minima(molecule, minima)
    terms, values = list(normalized), list(normalized.values())
    fixed_terms = []
    for term, angles in normalized.items():
        if not angles:
            if settings.empty_minima == "error":
                raise ValueError("no_minima for torsion " + "-".join(map(str, term)))
            rd = rdkit_reference(molecule, positions(molecule))
            angle = float(rdMolTransforms.GetDihedralDeg(rd.GetConformer(), *term))
            if not math.isfinite(angle):
                raise ValueError("Undefined reference angle for empty minima")
            angles.append((angle + 180) % 360 - 180)
            fixed_terms.append(dict(torsion=list(term), reference_angle=angles[0]))
    total = math.prod(len(v) for v in values)
    sampled = combination_indices(total, limit, settings.seed)
    reference = positions(molecule)
    reference_rd = rdkit_reference(molecule, reference)
    contact = ContactCheck(molecule, reference, settings.clash_scale)
    rmsd_indices = [
        i
        for i, a in enumerate(molecule.Atoms)
        if settings.rmsd_atoms == "all" or a.elem != "H"
    ]
    if not rmsd_indices:
        raise ValueError("No atoms selected for RMSD")
    retained = [("reference", reference)]
    candidates, plan = {}, []
    for ordinal, index in enumerate(sampled, 1):
        angles = decode_combination(index, values)
        name = f"combination-{ordinal}"
        entry = dict(
            name=name,
            product_index=index,
            angles={"-".join(map(str, t)): v for t, v in zip(terms, angles)},
        )
        plan.append(entry)
        try:
            xyz = set_torsions(reference_rd, terms, angles)
            violation = contact.violation(xyz)
            if violation is not None:
                entry.update(
                    status="rejected_contact_before_optimization", contact=violation
                )
                continue
        except (ValueError, RuntimeError, FloatingPointError) as exc:
            entry.update(status="rejected_geometry", error=str(exc))
            continue
        if optimizer is not None:
            candidate = copy.deepcopy(molecule)
            signature = graph_signature(candidate)
            for atom, point in zip(candidate.Atoms, xyz):
                atom.coor = point.tolist()
            optimizer(candidate)
            if graph_signature(candidate) != signature:
                raise ValueError("Q4 optimizer changed atom order or molecular graph")
            xyz = positions(candidate)
            violation = contact.violation(xyz)
            if violation is not None:
                entry.update(
                    status="rejected_contact_after_optimization", contact=violation
                )
                continue
        distances = [aligned_rmsd(xyz, other, rmsd_indices) for _, other in retained]
        nearest = int(np.argmin(distances))
        entry.update(nearest_rmsd=distances[nearest], nearest=retained[nearest][0])
        if distances[nearest] <= settings.rmsd_threshold:
            entry["status"] = "rejected_rmsd_duplicate"
            continue
        candidates[name] = xyz
        retained.append((name, xyz))
        entry["status"] = "accepted"
    return candidates, dict(
        fixed_reference_torsions=fixed_terms,
        candidate_plan=plan,
        total_combinations=total,
        generated_count=len(sampled),
        accepted_count=len(candidates),
        limit=limit,
        settings=asdict(settings),
        normalized_minima={"-".join(map(str, t)): v for t, v in normalized.items()},
        enumeration="all" if total <= limit else "seeded_uniform_without_replacement",
        optimization="none" if optimizer is None else "external_callback",
        deduplication="proper-rotation aligned RMSD, fixed atom correspondence; no symmetry remapping",
    )
