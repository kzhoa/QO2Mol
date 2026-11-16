<div style="
  position: relative;
  width: 100%;
  padding-top: 56.25%;
  margin-bottom: 20px;
  background: #f0f0f0 url('assets/images/banner_960.jpg') center/contain no-repeat;
  background-size: cover;
">
  <img src="assets/images/banner_1.5k.png"
       alt=""
       style="
         position: absolute;
         top: 0;
         left: 0;
         width: 100%;
         height: 100%;
         opacity: 0;
       "
       loading="lazy">
</div>

<h1 align="center">🌟QO2Mol Dataset</h1>

<p align="center">
<img alt="" src="https://img.shields.io/badge/python->=3.11-blue" style="max-width: 100%;">
<a href="https://developer.nvidia.com/cuda-downloads"><img src="https://img.shields.io/badge/GPU-CUDA%20Recommended-76B900?style=flat-square&logo=nvidia" alt="CUDA"/></a>
<!-- <img alt="" src="https://img.shields.io/badge/license-CC_BY--NC--SA_4.0-blue" style="max-width: 100%;"> -->
</p>




 > The latest version **`v2.0.0`**, featuring more reasonable conformations and more precise structures.
 >
>For datasets up to v1.2.0, follow the instructions in [legacy/readme.md](legacy/readme.md).

QO2Mol is a dataset designed for advanced molecular research with:
- A diverse collection of molecular structures, over 20 million configurations.
- Comprehensive elemental coverage, including 10 essential elements (C, H, O, N, S, P, F, Cl, Br, I) to facilitate broad studies.
- Key quantum mechanical properties, optimized and computed at the B3LYP/def2-SVP level.
- An ideal platform for AI-driven molecular research and the exploration of structure-property relationships.
- Open-source license.



# Quick Navigation

- [Prepare Data](#download)
- [Run Demo](#usage-demo)
- [FAQ](#faq)
- [License](#-license)


# Prepare Data

<a id="download"></a>

## Download Files

1. Download the `.zip` files from the [☁️ Zenodo Record](https://zenodo.org/records/23995922).
2. For baseline training, extract the ten Main chunks (0–9). Extract the OOD archives if you plan to use those sets.
3. Verify the MD5 checksums of the extracted `.pkl` files against [CHANGELOG.md](CHANGELOG.md).
4. Put the extracted `.pkl` files under `./download/raw/`.

For version details and file checksums, please refer to [CHANGELOG.md](CHANGELOG.md).


## Setup Environment

```bash
# Clone this repository
git clone https://github.com/kzhoa/QO2Mol.git
cd QO2Mol

# Install dependencies
pip install -r requirements.txt
```

The training entry requires PyTorch and qqtools 1.3.29 or later.


## 🐛 Process the Data

From the repository root, you can prepare the ten Main chunks before starting training. Training also creates missing Main caches on first use.

```bash
python -m src.dataset.process
```

> ⚠️ The first run might take 30–60 minutes — perfect time for a coffee break ☕️!

The script reads one pickle at a time from `download/raw/` and writes a tensor cache per chunk under `download/processed/`. Unpickling a chunk requires substantially more RAM than its file size. If interrupted, rerun the same command to reuse completed chunks.

To prepare a specific OOD set, pass `--scope ood_medium` or `--scope ood_large`. To prepare Main and both OOD sets, run `python -m src.dataset.process --scope all`. Use `--root /path/to/data` for a different directory containing `raw/`.



# Data File Format

Within each `.pkl` file, keys are listed as follows:

| Field               | Description                                                                                                                                          |
| :------------------ | :--------------------------------------------------------------------------------------------------------------------------------------------------- |
| inchikey            | String, the identity of the conformer.                                                                                                               |
| confid              | String, the identity of the conformer.                                                                                                               |
| atom_count          | Integer, the number of atoms in the molecule.                                                                                                        |
| bond_count          | Integer, the number of bonds in the molecule.                                                                                                        |
| elements            | List, length equal to the number of atoms. Each value indicates the atomic number in the periodic table.                                             |
| coordinates         | List, length equal to the number of atoms. Each element is a 3-tuple representing the 3D coordinates (x, y, z) of the corresponding atom.   Unit: Å. |
| edge_list           | List, length equal to the number of bonds multiplied by 2. Each element (i, j) represents an edge from atom i to atom j.                             |
| edge_attr           | List, length equal to the number of bonds multiplied by 2. Each value represents a bond type. '1': single bond, '2': double bond, '3': triple bond.  |
| energy              | Float, the calculated potential energy of the molecule. Unit: eV.                                                                                    |
| forces              | Cartesian forces in stored atom order, shape (N, 3). Each row gives the x, y, and z force components for one atom in eV/Å.                           |
| net_charge          | Scalar, the overall charge of a molecule.                                                                                                            |
| formal_charge       | List, length equal to the number of atoms. Each element represents the formal charge of the corresponding atom.                                      |
| zpe_correction      | Float, correction to zero-point vibrational energy.               Unit: Ha.                                                                          |
| thermal_corrections | Dictionary containing the thermal corrections to energy, enthalpy, and gibbs. These are correction terms rather than corrected totals. Unit: Ha.     |
| iso_polar           | Float, the static isotropic polarizability. Unit: Bohr³.                                                                                             |
| ese                 | Float, the electronic spatial extent, ⟨R²⟩. Unit: Bohr² (atomic units).                                                                              |
| dipole              | Tuple of three floats (μx, μy, μz), representing the dipole-moment components in the Gaussian input orientation. Unit: Debye.                        |


## Note on Units
We adopt the listed constants to convert units:
```python
1 Hartree = 27.21134 eV
1 Hartree/Bohr = 51.422 eV/Å
```
With these constants, you can easily switch back to atomic units and then convert to any other units you desire.

<a id="usage-demo"></a>

# Example Usage


## Train baseline

From the repository root, run the command for your model.

| Model | Training command |
| --- | --- |
| DimeNet++ | `sh scripts/dimenetpp/qm_benchmark/tr_baseline_e.sh` |
| SphereNet | `sh scripts/spherenet/tr_baseline_e.sh` |

The optional training code lives in `baselines/`. Both scripts run `python -m baselines.train` with a model configuration from `baselines/configs/` and forward command-line overrides.
The default data root is `download/`. Each baseline builds tensor caches from the ten Main chunks on first use.
Each chunk cache stores `offsets`, `z`, `pos`, `force`, `energy` (residual eV), and `ref_energy` (summed atomic reference eV). Their two energy fields sum to the raw SCF energy. Energy-only training uses the tensor cache directly. Full samples also expose conformer identifiers, read from the raw pickle on first access to that chunk.
Main chunk 9 is the test set. From chunks 0–8, the baseline randomly selects 10% of conformers (rounded down) for validation using the configured seed and uses the rest for training. Energy standardization uses the training split. The fixed atomic reference energies were fitted on all ten Main chunks.

To report test energy MAE in eV, use the same config and dataset as the training run. For example, to evaluate a DimeNet++ checkpoint, run

```bash
python -m baselines.train --config baselines/configs/dimenetpp.yml --test --ckp /path/to/best.pt
```

<a id="faq"></a>

# FAQ

**1️⃣ Check Data Integrity (Super Important!)**

To make sure everything's as it should be, please **check the MD5 hash** against the one in [CHANGELOG.md](CHANGELOG.md). Most of the email questions we get are actually solved by this simple check!

**2️⃣ Verify File Paths**

Ensure all `.pkl` files are placed in the correct directory.

```bash
# The correct directory structure should be:
QO2Mol/
├── download/
│   └── raw/           # <-- .pkl files go here!
│       ├── file1.pkl
│       ├── file2.pkl
│       └── ...
├── src/
└── ...
```

**3️⃣ Check Your Python Environment**
```bash
# Check Python version
python --version
# We strongly recommend python>=3.11.9

# Check NumPy version
python -c "import numpy; print(numpy.__version__)"
# Should be >=2.1.1
```

If all checks pass but it still fails, please contact us with full error messages.



<!-- <a id="cite"></a>

# 🔗 How to Cite

```bash

``` -->


# License

### Code

The code in this repository is licensed under the MIT License.
See [LICENSE](LICENSE) for the full license text.

### Data

The datasets are hosted on Zenodo DOI:  [10.5281/zenodo.23995922](https://doi.org/10.5281/zenodo.23995922).
Please refer to the corresponding Zenodo record for the data license and terms of use.
