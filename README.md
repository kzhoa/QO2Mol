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
<img alt="" src="https://img.shields.io/badge/license-CC_BY--NC--SA_4.0-blue" style="max-width: 100%;">
</p>


Welcome to the Quantum Open Organic Molecular (QO2Mol)​ database!

**Project History:**
Established in 2024, this database is continuously refined with annual updates to its data and the removal of unreasonable structures.

🥳 In Feb. 2026, we released the latest version **`v1.3.0`**, featuring more reasonable conformations and more precise structures.

>We recommend that users of older versions upgrade to the latest data release. However, if you have a strong dependency on an older version, all historical code has been preserved in the `old_versions` directory, allowing you to continue running them.

✨ QO2Mol is a curated dataset designed for advanced molecular research with:
- A diverse collection of molecular structures, over 20 million configurations.
- Comprehensive elemental coverage, including 10 essential elements (C, H, O, N, S, P, F, Cl, Br, I) to facilitate broad studies.
- Key quantum mechanical properties, optimized and computed at the B3LYP/def2-SVP level.
- An ideal platform for AI-driven molecular research and the exploration of structure-property relationships.
- Open-source license.



# 📑 Quick Navigation

- [⏬ Prepare Data](#download)
- [💫 Run Demo](#usage-demo)
- [❓ FAQ](#faq)
- [🤝 Contact & Contribute](#contact-us)
- [📄 License](#-license)
<!-- - [🔗 How to Cite](#cite) -->


# ⏬ Prepare Data

<a id="download"></a>

## ☁️ Download Files

 Download the dataset files from our ☁️Google Drive: [Google Drive](https://drive.google.com/drive/folders/1-4FrnNrVBlL2RaBuXpalgNCk1q79VHtc?usp=drive_link)

>💡 Tip:
Check out [CHANGELOG.md](CHANGELOG.md) for the latest updates and version details!
>
>The latest version now is `v1.3.0`.

After download , put all `*.pkl` files under `./download/raw/` directory.
Make md5 hash check with values in `CHANGELOG.md`.


## ⚙️ Setup Your Environment

```bash
# Clone this repository
git clone https://github.com/kzhoa/QO2Mol.git
cd QO2Mol

# Install dependencies
pip install -r requirements.txt
```

>Note that `torch_geometric` may need to be installed separately follow the instruction on [PyG Documentation](https://pytorch-geometric.readthedocs.io/en/latest/).


## 🐛 Process the Data

Make sure you're in `/your/path/QO2Mol`

Then run the data processing script.
```bash
sh process_data.sh
```


>⚠️This process might take 30-60 minutes — perfect time for a coffee break☕️!


>💡 If you encounter an error, please check the [FAQ](#faq) first.



# Data File Format

Within each `.pkl` file, keys are listed as follows:

| Field         | Description                                                                                                                                          |
| :------------ | :--------------------------------------------------------------------------------------------------------------------------------------------------- |
| inchikey      | String, the identity of the conformer.                                                                                                               |
| confid        | String, the identity of the conformer.                                                                                                               |
| atom_count    | Integer, the number of atoms in the molecule.                                                                                                        |
| bond_count    | Integer, the number of bonds in the molecule.                                                                                                        |
| elements      | List, length equal to the number of atoms. Each value indicates the atomic number in the periodic table.                                             |
| coordinates   | List, length equal to the number of atoms. Each element is a 3-tuple representing the 3D coordinates (x, y, z) of the corresponding atom.   Unit: Å. |
| edge_list     | List, length equal to the number of bonds multiplied by 2. Each element (i, j) represents an edge from atom i to atom j.                             |
| edge_attr     | List, length equal to the number of bonds multiplied by 2. Each value represents a bond type. '1': single bond, '2': double bond, '3': triple bond.  |
| energy        | Float, the calculated potential energy of the molecule. Unit: eV.                                                                                    |
| force         | List[List[float]], shape (N, 3). Each element represents the force component (x, y, z) of an atom. Unit: eV/Å​ .                                     |
| net_charge    | Float, the overall charge of a molecule.                                                                                                             |
| formal_charge | List, length equal to the number of atoms. Each element represents the formal charge of the corresponding atom.                                      |
| ZPVE          | Float, zero-point vibrational energy.               Unit: Ha.                                                                                        |


## Note on Units
We adopt the listed constants to convert units:
```python
1 Hartree = 27.21134 eV
1 Hartree/Bohr = 51.422 eV/Å
```
With these constants, you can easily switch back to atomic units and then convert to any other units you desire.

<a id="usage-demo"></a>

# 💫 Example Usage


## Train a DimeNet Model
We provide a demo for training `dimenet` model on `QO2Mol`.


First Navigate to the repository root.
Execute the training script using the following command:
```bash
sh scripts/dimenetpp/qm_benchmark/tr_baseline_e.sh
```

Upon successful completion, a training progress bar will appear in your terminal.


## Explore SphereNet Alternative

For users interested in exploring an alternative model, a comparable training script for SphereNet is available. This script allows you to train a SphereNet model on the `QO2MOL` dataset, similar to the DimeNet example.

Navigate to the repository root directory.
Execute the SphereNet training script using the following command:
```bash
sh scripts/spherenet/tr_baseline_e.sh
```
We encourage you to experiment with this script to compare performance or explore different model architectures.


<a id="faq"></a>

# ❓ FAQ

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

If all checks pass but it still fails, please [contact us](#contact-us) with full error messages.



<!-- <a id="cite"></a>

# 🔗 How to Cite

```bash

``` -->


<a id="contact-us"></a>

# 🤝 Get Help & Contribute


We are excited about your interest in our project! Whether you have a question, need assistance, or want to contribute, we're here to help and collaborate.

**For any questions, suggestions, or to discuss potential collaborations, please feel free to reach out to the first author.** Your input is highly valued, and we look forward to hearing from you!


**First Author Contact:**
*   **Email:** [liuwq20@fudan.edu.cn](mailto:liuwq20@fudan.edu.cn)
*   **Preferred Title:** `Mr. qq`


# 📄 License

This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. The images or other third party material in this article are included in the article’s Creative Commons license, unless indicated otherwise in the credit line; if the material is not included under the Creative Commons license, users will need to obtain permission from the license holder to reproduce the material. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
