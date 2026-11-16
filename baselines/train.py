import math

import torch
from qqtools.plugins.qpipeline.entry import create_pipeline

from baselines.models.dimenetpp import DimeNetPlusPlus
from baselines.models.spherenet import SphereNet
from baselines.task import QO2MolEnergyTask


def prepare_model(args):
    config = dict(args.model)
    model = {"dimenetpp": DimeNetPlusPlus, "spherenet": SphereNet}[config.pop("name")]
    return model(**config)


def configure_epoch_schedule(args):
    params = args.optim.scheduler_params
    if args.optim.scheduler != "lambda" or "warmup_fraction" not in params:
        return
    epochs = args.runner.max_epochs
    warmup = int(epochs * params.warmup_fraction)
    initial, minimum = params.initial_factor, params.min_factor

    def lr_factor(epoch):
        if epoch < warmup:
            return initial + (1 - initial) * epoch / warmup
        phase = math.pi * min(epoch, epochs) / epochs
        return minimum + (1 - minimum) * (1 + math.cos(phase)) / 2

    params.lr_lambda = lr_factor


def main():
    pipe = create_pipeline(prepare_model, QO2MolEnergyTask)
    if pipe.mode == "test":
        if pipe.args.ckp_file:
            checkpoint = torch.load(pipe.args.ckp_file, map_location="cpu", weights_only=True, mmap=True)
            if "task_state_dict" in checkpoint:
                pipe.task.load_state_dict(checkpoint["task_state_dict"])
            del checkpoint
        result = pipe.evaluate_once(pipe.task.test_loader)
        if result is not None:
            score = result.target_value("test_metric")
            pipe.logger.info(f"Test energy MAE: {score:.6f} eV")
    else:
        configure_epoch_schedule(pipe.args)
        pipe.fit()


if __name__ == "__main__":
    main()
