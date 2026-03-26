#!/usr/bin/env python3
"""Export the finite-degree GV proof constants to CSV tables.

This script reads the rigorous certification output in
`results/finite_gv_certificates.json` and writes two CSV files matching the
constants used in Appendix D (HA side) and Appendix E (MN side).
"""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
INPUT_JSON = RESULTS / "finite_gv_certificates.json"
HA_CSV = RESULTS / "finite_gv_ha_constants.csv"
MN_CSV = RESULTS / "finite_gv_mn_constants.csv"


def load_targets() -> list[dict[str, object]]:
    data = json.loads(INPUT_JSON.read_text())
    if not isinstance(data, dict) or "targets" not in data or not isinstance(data["targets"], list):
        raise ValueError(f"Unexpected certificate format in {INPUT_JSON}")
    return data["targets"]


def format_tuple(triple: list[int]) -> str:
    return f"({triple[0]},{triple[1]},{triple[2]})"


def build_ha_rows(targets: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for target in targets:
        triple = target["triple"]
        ha = target["ha"]
        delta_upper = target["delta_upper_certificate"]
        rows.append(
            {
                "tuple": format_tuple(triple),
                "j_z": triple[0],
                "j_x": triple[1],
                "k": triple[2],
                "beta_z": target["beta_z"],
                "delta_bar": delta_upper["delta_gv_upper"],
                "lambda_z": ha["small_input_lambda"],
                "epsilon_z": -ha["compact_strip_certificate"]["worst_upper_bound"],
            }
        )
    rows.sort(key=lambda row: (row["k"], row["j_z"], row["j_x"]))
    return rows


def build_mn_rows(targets: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for target in targets:
        triple = target["triple"]
        mn = target["mn"]
        rows.append(
            {
                "tuple": format_tuple(triple),
                "j_z": triple[0],
                "j_x": triple[1],
                "k": triple[2],
                "beta_x": target["beta_x"],
                "B_x": mn["small_support_base"],
                "epsilon_x": -mn["large_support_certificate"]["worst_upper_bound"],
            }
        )
    rows.sort(key=lambda row: (row["k"], row["j_z"], row["j_x"]))
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    targets = load_targets()
    ha_rows = build_ha_rows(targets)
    mn_rows = build_mn_rows(targets)
    write_csv(HA_CSV, ha_rows)
    write_csv(MN_CSV, mn_rows)
    print(f"wrote {HA_CSV} ({len(ha_rows)} rows)")
    print(f"wrote {MN_CSV} ({len(mn_rows)} rows)")


if __name__ == "__main__":
    main()
