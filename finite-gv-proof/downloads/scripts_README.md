# Proof Scripts

This directory contains the scripts used for the finite-degree GV
certification in Appendices D and E.

## Core certification script

- `certify_finite_gv_triples.py`

This is the main rigorous computer-assisted proof script. It uses
validated numerics based on interval arithmetic and adaptive subdivision to
certify the finite-domain negativity that remains after the analytic
first-moment reductions.

It produces:

- `results/finite_gv_certificates.json`

The JSON file records, for each certified triple:

- the certified upper bound `delta_bar` used on the HA side,
- the HA small-input constant `lambda_z`,
- the certified HA compact-strip margin `epsilon_z`,
- the MN small-support base `B_x`,
- the certified MN large-support margin `epsilon_x`.

## Table export script

- `export_finite_gv_proof_tables.py`

This script converts `results/finite_gv_certificates.json` into two CSV tables
matching the constants used in the paper:

- `results/finite_gv_ha_constants.csv`
- `results/finite_gv_mn_constants.csv`

## Optional plotting scripts

These are not needed for the proof itself, but were used to visualize the scan
window and the linear-distance proxies:

- `plot_balanced_delta_gv_curve.py`
- `plot_balanced_side_gv_jz_k.py`
- `plot_balanced_delta_lin_curves.py`

## Requirements

The proof scripts require:

- `mpmath`

The minimal pip requirements are listed in:

- `requirements-proof.txt`

The plotting scripts additionally require:

- `matplotlib`

## Minimal usage

From the repository root:

```bash
python3 scripts/certify_finite_gv_triples.py
python3 scripts/export_finite_gv_proof_tables.py
```

After these commands, the certificate JSON and the two CSV tables will be
available under `results/`.
