import numpy as np

from .qo2mol import QO2MolDataset


def prepare_train_dataset(root_dir, split_valid=True, seed=None, *, energy_only=False):
    dataset = QO2MolDataset(root_dir, energy_only=energy_only)
    # ends contains exclusive global offsets in Main chunk order.
    test_start = dataset.ends[-2]
    test = dataset[test_start:]
    if not split_valid:
        return dataset[:test_start], test
    rng = np.random if seed is None else np.random.RandomState(seed)
    indices = rng.permutation(test_start)
    valid_size = test_start // 10
    return dataset[indices[valid_size:]], dataset[indices[:valid_size]], test
