#!/bin/sh
exec python -m baselines.train --config baselines/configs/spherenet.yml "$@"
