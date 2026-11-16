import argparse
import logging

from .qo2mol import QO2MolDataset


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare QO2Mol tensor caches.")
    parser.add_argument("--root", default="download")
    scopes = ("main", "ood_medium", "ood_large")
    parser.add_argument("--scope", choices=(*scopes, "all"), default="main",
                        help="Dataset to prepare (default: main, used by both baselines)")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    for scope in scopes if args.scope == "all" else (args.scope,):
        dataset = QO2MolDataset(args.root, scope=scope, energy_only=True)
        print(f"{scope}: {len(dataset):,} samples ready in {dataset.processed_dir}", flush=True)
