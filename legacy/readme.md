# Legacy QO2Mol datasets (through v1.2.0)

This directory preserves the PyTorch Geometric data loader, models, and training code for QO2Mol dataset releases up to and including **v1.2.0**. If you have files named `QO2Mol_main_chunk_0.pkl` through `QO2Mol_main_chunk_9.pkl`, this is the loader for that naming scheme. For the current release, follow the [main README](../README.md).

Run all commands below from the root of the current repository.

## Install PyTorch Geometric (PyG)

Install the repository dependencies with `python -m pip install -r requirements.txt`. Legacy preprocessing and training also require PyG. Follow the [official PyG installation guide](https://pytorch-geometric.readthedocs.io/en/latest/install/installation.html) to install `torch-geometric` (2.7.0 works with NumPy 2) and, for training, the `torch-cluster`, `torch-scatter`, and `torch-sparse` extensions. Choose extension wheels that match your PyTorch version and CPU or CUDA build.

## Process legacy data

Place the ten Main pickle files and `QO2Mol_set_b.pkl` under `download/raw/`, then run:

```bash
python legacy/dataset/qo2mol/qo2mol.py
```

The loader writes `QO2Mol_main_processed.pt` and `QO2Mol_b_processed.pt` under `download/processed/`. See the [changelog](../CHANGELOG.md) for the files updated in each historical release.

## Train a legacy energy baseline

Choose the model you want to train. Both commands use the legacy Main data and its original random train/validation/test split.

```bash
python -m legacy.entry.main_dimenet --config legacy/configs/dimenetpp.yml
# Or train SphereNet:
python -m legacy.entry.main_spherenet --config legacy/configs/spherenet.yml
```

## Optional historical checkout

The preprocessing and training commands above run in the current checkout. To inspect the source as it stood at the v1.2.0 dataset update, use the commit hash `a31b08d`. This command creates a separate directory named `QO2Mol-history` next to the repository.

```bash
git worktree add --detach ../QO2Mol-history a31b08d
cd ../QO2Mol-history
```
