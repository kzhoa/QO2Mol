"""Periodic sampled-curve minima; no inference of unconstrained stationary points."""

import math
from dataclasses import asdict
from ...options import ScanOptions


def _analyze_records(records, *, options=None, return_report=False):
    settings = options or ScanOptions()
    groups = {}
    for r in records:
        torsion = r["torsion"]
        if (
            len(torsion) != 4
            or len(set(torsion)) != 4
            or any(type(i) is not int or i < 0 for i in torsion)
        ):
            raise ValueError("Four distinct zero-based torsion indices required")
        key = (str(r["id"]), "-".join(map(str, torsion)))
        groups.setdefault(key, []).append(r)
    result, reports = {}, []
    for (name, term), rows in groups.items():
        # Missing metadata is allowed only consistently; the caller is then
        # responsible for common calculation settings and energy units.
        for field in ("level", "charge", "multiplicity", "energy_unit"):
            values = [r.get(field) for r in rows]
            if any(value != values[0] for value in values):
                raise ValueError(f"Inconsistent {field} in scan {name}/{term}")
        points, failures = [], 0
        for r in rows:
            if "converged" in r and type(r["converged"]) is not bool:
                raise ValueError("converged must be a boolean")
            if r.get("converged") is False:
                failures += 1
                continue
            angle, energy = float(r["angle"]), float(r["energy"])
            if not math.isfinite(angle) or not math.isfinite(energy):
                raise ValueError("Finite scan energies/angles required")
            points.append(((angle + 180) % 360 - 180, energy))
        points.sort()
        unique = []
        for angle, energy in points:
            match = next(
                (
                    p
                    for p in unique
                    if abs((angle - p[0] + 180) % 360 - 180) <= settings.angle_tolerance
                ),
                None,
            )
            if match is not None:
                if abs(energy - match[1]) > settings.energy_tolerance:
                    raise ValueError(
                        "Conflicting energies at a periodic duplicate angle"
                    )
                match[1] = min(match[1], energy)
            else:
                unique.append([angle, energy])
        n = len(unique)
        if n < 3 or any(
            abs((unique[(i + 1) % n][0] - unique[i][0]) % 360 - settings.expected_step)
            > settings.angle_tolerance
            for i in range(n)
        ):
            raise ValueError(
                f"Incomplete or irregular periodic scan {name}/{term}; expected step {settings.expected_step} degrees"
            )
        # Connected near-equal runs, including runs crossing the periodic seam.
        unseen, plateaus = set(range(n)), []
        while unseen:
            seed = min(unseen)
            unseen.remove(seed)
            component, pending = [seed], [seed]
            while pending:
                i = pending.pop()
                for j in ((i - 1) % n, (i + 1) % n):
                    if (
                        j in unseen
                        and abs(unique[i][1] - unique[j][1])
                        <= settings.energy_tolerance
                    ):
                        unseen.remove(j)
                        component.append(j)
                        pending.append(j)
            plateaus.append(component)
        selected = []
        for component in plateaus:
            energies = [unique[i][1] for i in component]
            outside = {
                j
                for i in component
                for j in ((i - 1) % n, (i + 1) % n)
                if j not in component
            }
            if (
                outside
                and max(energies) - min(energies) <= settings.energy_tolerance
                and all(
                    unique[j][1] > max(energies) + settings.energy_tolerance
                    for j in outside
                )
            ):
                representative = min(
                    component, key=lambda i: (unique[i][1], unique[i][0])
                )
                selected.append(unique[representative])
        selected.sort(key=lambda p: (p[1], p[0]))
        result.setdefault(name, {})[term] = [p[0] for p in selected]
        reports.append(
            dict(
                id=name,
                torsion=term,
                unique_points=n,
                failed_points=failures,
                duplicate_points=len(points) - n,
                minima=selected,
                convergence_assumed=any("converged" not in r for r in rows),
            )
        )
    report = dict(
        settings=asdict(settings),
        curves=reports,
        interpretation="minima on a complete periodic constrained scan; saved labels are ignored",
    )
    return dict(minima=result, report=report) if return_report else result


def minima_from_scan_records(records, *, options=None, return_report=False):
    """Collect mode excludes the entire molecule if any of its curves fails."""
    settings = options or ScanOptions()
    if settings.error_policy == "collect" and not return_report:
        raise ValueError("collect mode requires return_report=True to retain failures")
    molecules = {}
    for row in records:
        # An unidentifiable record cannot safely be isolated to one molecule.
        if "id" not in row:
            raise ValueError("Scan record requires id for failure isolation")
        molecules.setdefault(str(row["id"]), []).append(row)
    minima, curves, failed, provenance, calculation_settings = {}, [], {}, {}, {}
    for name, rows in molecules.items():
        try:
            for field in (
                "level",
                "charge",
                "multiplicity",
                "energy_unit",
                "provenance",
            ):
                if any(r.get(field) != rows[0].get(field) for r in rows):
                    raise ValueError(
                        f"Inconsistent {field} across molecule scan curves"
                    )
            saved = rows[0].get("provenance")
            if saved is not None:
                if not isinstance(saved, dict) or str(saved.get("id")) != name:
                    raise ValueError("Invalid scan provenance identity")
                for field in ("ordered_graph_sha256", "reference_sha256"):
                    value = saved.get(field)
                    if (
                        not isinstance(value, str)
                        or len(value) != 64
                        or any(c not in "0123456789abcdef" for c in value)
                    ):
                        raise ValueError(f"Invalid provenance {field}")
                for field in ("level", "charge", "multiplicity", "energy_unit"):
                    if (
                        field in saved
                        and field in rows[0]
                        and saved[field] != rows[0][field]
                    ):
                        raise ValueError(f"Conflicting provenance {field}")
                if "expected_torsions" in saved:
                    expected = {tuple(t) for t in saved["expected_torsions"]}
                    actual = {tuple(r["torsion"]) for r in rows}
                    if expected != actual:
                        raise ValueError("Missing or unexpected torsion curves")
            outcome = _analyze_records(rows, options=settings, return_report=True)
            minima.update(outcome["minima"])
            curves.extend(outcome["report"]["curves"])
            calculation_settings[name] = {
                field: rows[0][field]
                for field in ("level", "charge", "multiplicity", "energy_unit")
                if field in rows[0]
            }
            if saved is not None:
                provenance[name] = saved
        except (ValueError, TypeError, KeyError, IndexError) as exc:
            if settings.error_policy == "strict":
                raise ValueError(f"Scan molecule {name}: {exc}") from exc
            failed[name] = str(exc)
    report = dict(
        settings=asdict(settings),
        curves=curves,
        failed_molecules=failed,
        provenance=provenance,
        calculation_settings=calculation_settings,
        status="partial_failure" if failed and minima else "errors" if failed else "ok",
        interpretation="complete periodic constrained scans; failure isolation is per molecule",
    )
    return dict(minima=minima, report=report) if return_report else minima
