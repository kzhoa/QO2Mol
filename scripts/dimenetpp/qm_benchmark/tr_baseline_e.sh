#!/bin/sh
exec python -m baselines.train --config baselines/configs/dimenetpp.yml "$@"
