#!/bin/sh
set -eu
: "${1:?Pass the checkpoint .pt file as the first argument}"
checkpoint=$1
shift
exec python -m baselines.train --test --config baselines/configs/dimenetpp.yml --ckp "$checkpoint" "$@"
