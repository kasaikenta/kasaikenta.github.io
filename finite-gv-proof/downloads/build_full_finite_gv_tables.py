#!/usr/bin/env python3
"""Build full finite-degree GV constant tables for all certified triples in T_scan."""

from __future__ import annotations

import csv
from pathlib import Path

from certify_finite_gv_triples import TripleCertifier, TripleSpec


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
STATUS_CSV = RESULTS / "balanced_side_gv_jz_k_status.csv"
HA_CSV = RESULTS / "finite_gv_ha_constants.csv"
MN_CSV = RESULTS / "finite_gv_mn_constants.csv"
HA_ROWS_TEX = RESULTS / "finite_gv_ha_rows.tex"
MN_ROWS_TEX = RESULTS / "finite_gv_mn_rows.tex"


def load_certified_specs() -> list[TripleSpec]:
    rows = list(csv.DictReader(STATUS_CSV.open()))
    specs: list[TripleSpec] = []
    for row in rows:
        if row["ha_certified"] != "True" or row["mn_certified"] != "True":
            continue
        specs.append(
            TripleSpec(
                int(row["j_z"]),
                int(row["j_x"]),
                int(row["k"]),
                beta_z=float(row["ha_beta_z"]),
                beta_x=float(row["mn_beta_x"]),
            )
        )
    specs.sort(key=lambda spec: (spec.k, spec.j_z, spec.j_x))
    return specs


def latex_sci(x: float) -> str:
    mantissa, exp = f"{x:.4e}".split("e")
    exponent = int(exp)
    return rf"${mantissa}\times10^{{{exponent}}}$"


def tuple_label(spec: TripleSpec) -> str:
    return rf"$({spec.j_z},{spec.j_x},{spec.k})$"


def build_payloads(specs: list[TripleSpec]) -> list[dict[str, object]]:
    payloads: list[dict[str, object]] = []
    for spec in specs:
        cert = TripleCertifier(spec)
        payloads.append(cert.certify())
    return payloads


def write_csvs(payloads: list[dict[str, object]]) -> None:
    ha_rows: list[dict[str, object]] = []
    mn_rows: list[dict[str, object]] = []
    for payload in payloads:
        triple = payload["triple"]
        ha = payload["ha"]
        mn = payload["mn"]
        delta_upper = payload["delta_upper_certificate"]
        ha_rows.append(
            {
                "tuple": f"({triple[0]},{triple[1]},{triple[2]})",
                "j_z": triple[0],
                "j_x": triple[1],
                "k": triple[2],
                "beta_z": payload["beta_z"],
                "delta_bar": delta_upper["delta_gv_upper"],
                "lambda_z": ha["small_input_lambda"],
                "epsilon_z": -ha["compact_strip_certificate"]["worst_upper_bound"],
            }
        )
        mn_rows.append(
            {
                "tuple": f"({triple[0]},{triple[1]},{triple[2]})",
                "j_z": triple[0],
                "j_x": triple[1],
                "k": triple[2],
                "beta_x": payload["beta_x"],
                "B_x": mn["small_support_base"],
                "epsilon_x": -mn["large_support_certificate"]["worst_upper_bound"],
            }
        )

    with HA_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ha_rows[0].keys()))
        writer.writeheader()
        writer.writerows(ha_rows)

    with MN_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(mn_rows[0].keys()))
        writer.writeheader()
        writer.writerows(mn_rows)


def write_tex_rows(payloads: list[dict[str, object]]) -> None:
    ha_lines: list[str] = []
    mn_lines: list[str] = []
    for payload in payloads:
        spec = TripleSpec(
            payload["triple"][0],
            payload["triple"][1],
            payload["triple"][2],
            beta_z=payload["beta_z"],
            beta_x=payload["beta_x"],
        )
        ha = payload["ha"]
        mn = payload["mn"]
        delta_upper = payload["delta_upper_certificate"]["delta_gv_upper"]
        ha_lines.append(
            " & ".join(
                [
                    tuple_label(spec),
                    f"{payload['beta_z']:.2f}",
                    f"{delta_upper:.8f}",
                    f"{ha['small_input_lambda']:.9f}",
                    latex_sci(-ha["compact_strip_certificate"]["worst_upper_bound"]),
                ]
            )
            + r" \\"
        )
        mn_lines.append(
            " & ".join(
                [
                    tuple_label(spec),
                    f"{payload['beta_x']:.2f}",
                    latex_sci(mn["small_support_base"]),
                    latex_sci(-mn["large_support_certificate"]["worst_upper_bound"]),
                ]
            )
            + r" \\"
        )

    HA_ROWS_TEX.write_text("\n".join(ha_lines) + "\n")
    MN_ROWS_TEX.write_text("\n".join(mn_lines) + "\n")


def main() -> None:
    specs = load_certified_specs()
    payloads = build_payloads(specs)
    write_csvs(payloads)
    write_tex_rows(payloads)
    print(f"wrote {HA_CSV} and {MN_CSV}")
    print(f"wrote {HA_ROWS_TEX} and {MN_ROWS_TEX}")


if __name__ == "__main__":
    main()
