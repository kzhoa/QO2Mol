# Figure 2: one command for products and a coloured PDF

From the `supplementary_tools` directory, install the package dependencies once:

```bash
python -m pip install -r requirements.txt
```

This includes CairoSVG for figure rendering. CairoSVG also needs the Cairo system library (e.g. `libcairo2` on Debian/Ubuntu); see the CairoSVG installation instructions for other platforms.

Then run:

```bash
python examples/figure2/run.py
```

The default input is `examples/figure2/input/figure2_3702543.smi`, the default backend is **f147**, and the result directory is `examples/figure2/output/`. These defaults are resolved relative to the script, so it can be launched from any working directory.

The output contains the standard fragmentation products (`parents.jsonl`, `occurrences.jsonl`, `candidates.jsonl`, `groups.jsonl`, `excluded.jsonl`, `errors.jsonl`, `candidates.sqlite`, `run.json`) plus:

- `figure2.pdf`: coloured vector PDF, suitable for enlargement and paper layout.
- `figure2.svg`: editable vector figure.
- `figure2.png`: preview image.
- `figure2.json`: labels, displayed SMILES, identities and elementary-block membership for auditing.

The parent is coloured by elementary-fragment membership. Red marks identify cut bonds. The following panels show independently capped elementary fragments and larger combinations, retaining the same colours. For figure2_3702543 with f147, these are 4 elementary fragments and 5 larger non-parent combinations; complete machine-readable results contain 21 occurrences and 10 unique structures, including a parent-equivalent TF product.



To use another single-parent input (SMILES, SDF/MOL, JSON or JSONL, optionally gzipped):

```bash
python examples/figure2/run.py --input examples/figure2/input/figure2_3702543.sdf --output examples/figure2/output_sdf
```



## Candidate inputs

All three candidates are in `input/`, each supplied as `.smi`, `.sdf`, `.mol`, `.json` and `.jsonl`. Names and titles use only these example identifiers; source-database annotations are not included.

| Example | Heavy atoms | Structural features |
| --- | --- | --- |
| `figure2_3702543` (default) | 21 | Ester, aminothiophene, phenyl and amide regions |
| `figure2_2356088` | 22 | Amide, pyrazole, fluorophenyl and alcohol regions |
| `figure2_1321483` | 19 | Sulfonamide, amide, cyclopropyl and fluorophenyl regions |

From `supplementary_tools`, render the alternatives with:

```bash
python examples/figure2/run.py --input examples/figure2/input/figure2_2356088.smi --output examples/figure2/output_2356088
python examples/figure2/run.py --input examples/figure2/input/figure2_1321483.smi --output examples/figure2/output_1321483
```
