import torch
from qqtools import qDictDataset
from qqtools.plugins.qpipeline import prepare_dataloder, qTaskBase

from src.dataset.load_dataset import prepare_train_dataset


def collate_energy(samples):
    fields = ("num_nodes", "z", "pos", "energy")
    return qDictDataset.collate_graph_samples(
        [{key: sample[key] for key in fields} for sample in samples],
        key_types={"node": ["z", "pos"], "edge": [], "graph": ["energy"]},
    )


class QO2MolEnergyTask(qTaskBase):
    def __init__(self, args):
        super().__init__()
        datasets = prepare_train_dataset(args.task.data_root, seed=args.seed, energy_only=True)
        self.mean, self.std = datasets[0].get_norm_factor() if args.task.standardize else (0., 1.)
        if self.std == 0:
            raise ValueError("Cannot standardize a constant energy target")
        self.train_loader, self.val_loader, self.test_loader = prepare_dataloder(
            *datasets, **args.task.dataloader, distributed=args.distributed, collate_fn=collate_energy,
        )

        if not args.test and len(self.train_loader) == 0:
            raise ValueError("Training split has fewer samples than one training batch")

    def batch_forward(self, model, batch_data):
        energy = model(batch_data["z"], batch_data["pos"], batch_data["batch"])
        return {"energy": energy.flatten()}

    def batch_loss(self, out, batch_data, loss_fn):
        target = (batch_data["energy"] - self.mean) / self.std
        return {"loss": (loss_fn(out["energy"], target), target.numel())}

    def batch_metric(self, out, batch_data):
        target = batch_data["energy"]
        energy = out["energy"].detach() * self.std + self.mean
        return {"e_mae": (torch.nn.functional.l1_loss(energy, target), target.numel())}

    def post_metrics_to_value(self, result, *, stage):
        return float(result["e_mae"])

    def state_dict(self):
        return {"mean": self.mean, "std": self.std}

    def load_state_dict(self, state):
        self.mean, self.std = state["mean"], state["std"]
