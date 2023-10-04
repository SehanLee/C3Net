# C3Net

[[arXiv]](https://arxiv.org/abs/2309.15334)

Continuum-medium continuous convolutional neural network for interaction potential in heterogeneous systems including pure solvents, 1-octation/water partitioning, PAMPA.

## Prediction

1. `cd` to [`prediction`](prediction).
2. Edit `solvent` and `solute` variables in [`main.py`](prediction/main.py).
3. Run [`main.py`](prediction/main.py) and the result will be written in `prediction.txt` under `output_dir`.

In the output `prediction.txt`, values under `Solvent_ID` column uniquely mean different `solvent`. For example:

`Solvent_ID` | `solvent`
--- | ---
... | ...
102 | water
103 | logP
104 | PAMPA

which follows the indices of `solv_prop_normalized` in [`element.py`](module/input/element.py).
