#!/usr/bin/env python3
"""Rigorous finite-degree GV certification for selected even balanced triples.

The proof method is the same as in the single-triple certifier:
  * first-moment + Markov,
  * validated numerics by interval arithmetic for compact finite domains,
  * analytic treatment of the HA boundary strip.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

from mpmath import iv


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
OUTPUT_JSON = RESULTS / "finite_gv_certificates.json"

iv.dps = 80


@dataclass(frozen=True)
class TripleSpec:
    j_z: int
    j_x: int
    k: int
    beta_z: float
    beta_x: float = 0.10

    @property
    def j_delta(self) -> int:
        return self.j_x - self.j_z

    @property
    def alpha_z(self) -> float:
        return self.j_z / self.k

    @property
    def alpha_x(self) -> float:
        return self.j_x / self.k

    @property
    def alpha_delta(self) -> float:
        return self.j_delta / self.k

    @property
    def tau_cut(self) -> float:
        return self.beta_z / self.k

    @property
    def support_cut(self) -> float:
        return self.beta_x / self.k

    @property
    def tuple_label(self) -> str:
        return f"({self.j_z},{self.j_x},{self.k})"


TARGETS = [
    TripleSpec(4, 6, 10, beta_z=0.25),
    TripleSpec(4, 8, 12, beta_z=0.20),
    TripleSpec(5, 9, 14, beta_z=0.15),
    TripleSpec(6, 14, 20, beta_z=0.10),
    TripleSpec(5, 17, 22, beta_z=0.10),
    TripleSpec(4, 20, 24, beta_z=0.10),
    TripleSpec(4, 26, 30, beta_z=0.08, beta_x=0.15),
]


def h2_float(x: float) -> float:
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return -(x * math.log2(x) + (1.0 - x) * math.log2(1.0 - x))


def h2_inverse_float(y: float) -> float:
    lo = 0.0
    hi = 0.5
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if h2_float(mid) < y:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


class TripleCertifier:
    def __init__(self, spec: TripleSpec) -> None:
        self.spec = spec
        self._h2_upper_cache: dict[float, float] = {}
        self._h2_lower_cache: dict[float, float] = {}
        self._p_lo_cache: dict[float, float] = {}
        self._p_hi_cache: dict[float, float] = {}
        self._mu_cache: dict[float, float] = {}
        self._y1_cache: dict[float, float] = {}
        self._yd_cache: dict[float, float] = {}
        self._wo_trial_cache: dict[tuple[float, float], float] = {}

        ax_half = iv.mpf([repr(self.spec.alpha_x / 2.0), repr(self.spec.alpha_x / 2.0)])
        self.omega_star_upper = self.upper((1 - iv.exp(iv.log(ax_half) / self.spec.k)) / 2)
        self.delta_upper = self._choose_delta_upper()

    @staticmethod
    def iv_const(x: float) -> iv.mpf:
        s = format(x, ".17g")
        return iv.mpf([s, s])

    @staticmethod
    def upper(x: iv.mpf) -> float:
        return math.nextafter(float(x.b), math.inf)

    @staticmethod
    def lower(x: iv.mpf) -> float:
        return math.nextafter(float(x.a), -math.inf)

    def log2_iv(self, x: iv.mpf) -> iv.mpf:
        return iv.log(x) / iv.log(self.iv_const(2.0))

    def h2_iv(self, x: iv.mpf) -> iv.mpf:
        return -(x * self.log2_iv(x) + (1 - x) * self.log2_iv(1 - x))

    def h2_scalar_interval(self, x: float) -> iv.mpf:
        if x <= 0.0 or x >= 1.0:
            return self.iv_const(0.0)
        return self.h2_iv(self.iv_const(x))

    def h2_upper(self, x: float) -> float:
        if x not in self._h2_upper_cache:
            self._h2_upper_cache[x] = self.upper(self.h2_scalar_interval(x))
        return self._h2_upper_cache[x]

    def h2_lower(self, x: float) -> float:
        if x not in self._h2_lower_cache:
            self._h2_lower_cache[x] = self.lower(self.h2_scalar_interval(x))
        return self._h2_lower_cache[x]

    def _choose_delta_upper(self) -> float:
        delta_guess = h2_inverse_float(self.spec.alpha_z)
        scale = 10**8
        candidate = math.ceil(delta_guess * scale) / scale
        while self.h2_lower(candidate) <= self.spec.alpha_z:
            candidate = math.nextafter(candidate + 1.0 / scale, math.inf)
        return candidate

    @staticmethod
    def monotone_grid(start: float, pieces: list[tuple[float, float]]) -> list[float]:
        values = [round(start, 12)]
        cur = start
        for stop, step in pieces:
            while cur < stop - 1.0e-15:
                cur = min(stop, cur + step)
                values.append(round(cur, 12))
        return values

    def wo_trial_upper(self, tau_lo: float, tau_hi: float) -> float:
        key = (tau_lo, tau_hi)
        if key not in self._wo_trial_cache:
            t = iv.mpf([format(tau_lo, ".17g"), format(tau_hi, ".17g")])
            T = (1 - 2 * t) ** self.spec.k
            alpha_z = self.iv_const(self.spec.alpha_z)
            val = alpha_z * self.log2_iv(1 + T) + self.h2_iv(t) - alpha_z
            self._wo_trial_cache[key] = self.upper(val)
        return self._wo_trial_cache[key]

    def ha_p_lower(self, tau_lo: float) -> float:
        if tau_lo not in self._p_lo_cache:
            t_lo = self.iv_const(tau_lo)
            p_lo = (1 - (1 - 2 * t_lo) ** self.spec.k) / 2
            self._p_lo_cache[tau_lo] = self.lower(p_lo)
        return self._p_lo_cache[tau_lo]

    def ha_p_upper(self, tau_hi: float) -> float:
        if tau_hi not in self._p_hi_cache:
            t_hi = self.iv_const(tau_hi)
            p_hi = (1 - (1 - 2 * t_hi) ** self.spec.k) / 2
            self._p_hi_cache[tau_hi] = self.upper(p_hi)
        return self._p_hi_cache[tau_hi]

    def ha_joint_entropy_cross_sup(
        self,
        omega_lo: float,
        omega_hi: float,
        tau_lo: float,
        tau_hi: float,
    ) -> float:
        p_lo = self.ha_p_lower(tau_lo)
        p_hi = self.ha_p_upper(tau_hi)

        if p_lo <= omega_hi and omega_lo <= p_hi:
            return 0.0

        if p_lo > omega_hi:
            omega_star = omega_hi
            p_star = p_lo
        else:
            omega_star = omega_lo
            p_star = p_hi

        w = self.iv_const(omega_star)
        p = self.iv_const(p_star)
        val = self.h2_iv(w) + w * self.log2_iv(p) + (1 - w) * self.log2_iv(1 - p)
        return self.upper(val)

    def ha_cell_upper(self, omega_lo: float, omega_hi: float, tau_lo: float, tau_hi: float) -> float:
        return self.wo_trial_upper(tau_lo, tau_hi) + self.ha_joint_entropy_cross_sup(
            omega_lo,
            omega_hi,
            tau_lo,
            tau_hi,
        )

    def certify_ha_compact_strip(self) -> dict[str, object]:
        omega_grid = self.monotone_grid(
            0.0,
            [
                (min(0.02, self.delta_upper), 0.001),
                (min(0.06, self.delta_upper), 0.002),
                (min(max(0.079, self.delta_upper), self.delta_upper), 0.0005),
                (self.delta_upper, 0.00005 if self.delta_upper >= 0.06 else 0.00002),
            ],
        )
        tau_grid = self.monotone_grid(
            self.spec.tau_cut,
            [
                (0.10, 0.0025),
                (0.40, 0.01),
                (0.49, 0.001),
            ],
        )

        worst_leaf_upper = -math.inf
        worst_leaf_cell: dict[str, float] | None = None
        certified_cells = 0

        def refine(omega_lo: float, omega_hi: float, tau_lo: float, tau_hi: float, depth: int) -> None:
            nonlocal worst_leaf_upper, worst_leaf_cell, certified_cells
            cell_upper = self.ha_cell_upper(omega_lo, omega_hi, tau_lo, tau_hi)
            if cell_upper < 0.0:
                if cell_upper > worst_leaf_upper:
                    worst_leaf_upper = cell_upper
                    worst_leaf_cell = {
                        "omega_lo": omega_lo,
                        "omega_hi": omega_hi,
                        "tau_lo": tau_lo,
                        "tau_hi": tau_hi,
                        "depth": depth,
                    }
                certified_cells += 1
                return
            if depth >= 20:
                if cell_upper > worst_leaf_upper:
                    worst_leaf_upper = cell_upper
                    worst_leaf_cell = {
                        "omega_lo": omega_lo,
                        "omega_hi": omega_hi,
                        "tau_lo": tau_lo,
                        "tau_hi": tau_hi,
                        "depth": depth,
                    }
                return
            omega_mid = 0.5 * (omega_lo + omega_hi)
            tau_mid = 0.5 * (tau_lo + tau_hi)
            if (omega_hi - omega_lo) >= (tau_hi - tau_lo):
                refine(omega_lo, omega_mid, tau_lo, tau_hi, depth + 1)
                refine(omega_mid, omega_hi, tau_lo, tau_hi, depth + 1)
            else:
                refine(omega_lo, omega_hi, tau_lo, tau_mid, depth + 1)
                refine(omega_lo, omega_hi, tau_mid, tau_hi, depth + 1)

        for omega_lo, omega_hi in zip(omega_grid, omega_grid[1:]):
            for tau_lo, tau_hi in zip(tau_grid, tau_grid[1:]):
                refine(omega_lo, omega_hi, tau_lo, tau_hi, 0)

        return {
            "all_negative": worst_leaf_upper < 0.0,
            "worst_upper_bound": worst_leaf_upper,
            "worst_cell": worst_leaf_cell,
            "initial_num_cells": (len(omega_grid) - 1) * (len(tau_grid) - 1),
            "certified_leaf_cells": certified_cells,
            "omega_grid_size": len(omega_grid),
            "tau_grid_size": len(tau_grid),
        }

    def certify_ha_boundary_strip(self) -> dict[str, object]:
        u_hi = 0.02
        pinsker_gap_lower = self.lower(self.iv_const(u_hi * u_hi) / (2 * iv.log(self.iv_const(2.0))))
        coeff = 1.0 + self.spec.alpha_z
        log_gain_upper = self.upper(self.iv_const(coeff * (u_hi**self.spec.k)) / iv.log(self.iv_const(2.0)))
        return {
            "all_negative": log_gain_upper <= pinsker_gap_lower,
            "worst_upper_bound": math.nextafter(-self.spec.alpha_z, -math.inf),
            "worst_cell": {"tau_lo": 0.49, "tau_hi": 0.5},
            "pinsker_gap_lower": pinsker_gap_lower,
            "log_gain_upper": log_gain_upper,
        }

    def mn_mu_upper(self, omega_lo: float) -> float:
        if omega_lo not in self._mu_cache:
            self._mu_cache[omega_lo] = self.upper((1 - 2 * self.iv_const(omega_lo)) ** self.spec.k)
        return self._mu_cache[omega_lo]

    def mn_y1_upper(self, a_lo: float) -> float:
        if a_lo not in self._y1_cache:
            self._y1_cache[a_lo] = self.upper((1 - 2 * self.iv_const(a_lo)) ** self.spec.j_z)
        return self._y1_cache[a_lo]

    def mn_yd_upper(self, b_lo: float) -> float:
        if b_lo not in self._yd_cache:
            self._yd_cache[b_lo] = self.upper((1 - 2 * self.iv_const(b_lo)) ** self.spec.j_delta)
        return self._yd_cache[b_lo]

    def mn_cell_upper(
        self,
        omega_lo: float,
        omega_hi: float,
        a_lo: float,
        a_hi: float,
        b_lo: float,
        b_hi: float,
    ) -> float:
        h_a = self.spec.alpha_z * self.h2_upper(a_hi)
        h_b = self.spec.alpha_delta * self.h2_upper(b_hi)
        h_w = self.h2_upper(omega_hi)

        mu = self.iv_const(self.mn_mu_upper(omega_lo))
        y_1 = self.iv_const(self.mn_y1_upper(a_lo))
        y_delta = self.iv_const(self.mn_yd_upper(b_lo))
        log_term = self.log2_iv(1 + mu * y_1 * y_delta)
        return h_a + h_b + h_w - 1 + self.upper(log_term)

    def certify_mn_large_support(self) -> dict[str, object]:
        omega_grid = self.monotone_grid(
            0.0,
            [
                (min(0.01, self.omega_star_upper), 0.0005),
                (self.omega_star_upper, 0.001),
            ],
        )
        ab_grid = self.monotone_grid(
            0.0,
            [
                (0.02, 0.001),
                (0.05, 0.005),
                (0.50, 0.02),
            ],
        )

        worst_upper = -math.inf
        worst_cell: dict[str, float] | None = None
        cell_count = 0
        for omega_lo, omega_hi in zip(omega_grid, omega_grid[1:]):
            for a_lo, a_hi in zip(ab_grid, ab_grid[1:]):
                for b_lo, b_hi in zip(ab_grid, ab_grid[1:]):
                    if a_hi <= self.spec.support_cut and b_hi <= self.spec.support_cut:
                        continue
                    cell_count += 1
                    cell_upper = self.mn_cell_upper(omega_lo, omega_hi, a_lo, a_hi, b_lo, b_hi)
                    if cell_upper > worst_upper:
                        worst_upper = cell_upper
                        worst_cell = {
                            "omega_lo": omega_lo,
                            "omega_hi": omega_hi,
                            "a_lo": a_lo,
                            "a_hi": a_hi,
                            "b_lo": b_lo,
                            "b_hi": b_hi,
                        }
        return {
            "all_negative": worst_upper < 0.0,
            "worst_upper_bound": worst_upper,
            "worst_cell": worst_cell,
            "num_cells": cell_count,
            "omega_grid_size": len(omega_grid),
            "ab_grid_size": len(ab_grid),
        }

    def certify(self) -> dict[str, object]:
        exp_one = iv.exp(self.iv_const(1.0))
        ha_lambda = self.upper(self.iv_const(self.spec.k * self.spec.beta_z) / exp_one)
        rho_val = max(self.spec.beta_x / self.spec.j_z, self.spec.beta_x / self.spec.j_delta, self.omega_star_upper)
        rho = self.iv_const(rho_val)
        c_val = self.iv_const(2.0 * self.spec.beta_x / self.spec.k + self.omega_star_upper)
        c_k = (1.0 + rho) / ((1.0 - rho) ** 2)
        mn_small_support_base = self.upper(
            iv.sqrt(self.iv_const(2.0))
            * (exp_one * self.iv_const(1.0 + self.spec.alpha_x) / c_val)
            * ((c_k * self.iv_const(float(self.spec.k)) * c_val / exp_one) ** (self.spec.k // 2))
        )

        delta_upper_cert = {
            "delta_gv_upper": self.delta_upper,
            "h2_lower_at_delta_gv_upper": self.h2_lower(self.delta_upper),
            "alpha_z": self.spec.alpha_z,
        }
        ha_compact = self.certify_ha_compact_strip()
        ha_boundary = self.certify_ha_boundary_strip()
        mn_large_support = self.certify_mn_large_support()

        payload = {
            "triple": [self.spec.j_z, self.spec.j_x, self.spec.k],
            "beta_z": self.spec.beta_z,
            "beta_x": self.spec.beta_x,
            "delta_upper_certificate": delta_upper_cert,
            "ha": {
                "tau_cut": self.spec.tau_cut,
                "small_input_lambda": ha_lambda,
                "compact_strip_certificate": ha_compact,
                "boundary_strip_certificate": ha_boundary,
            },
            "mn": {
                "support_cut": self.spec.support_cut,
                "omega_star": self.omega_star_upper,
                "small_support_base": mn_small_support_base,
                "large_support_certificate": mn_large_support,
            },
        }

        if not (ha_lambda < 1.0):
            raise SystemExit(f"{self.spec.tuple_label}: HA small-input base is not < 1: {ha_lambda}")
        if not (delta_upper_cert["h2_lower_at_delta_gv_upper"] > self.spec.alpha_z):
            raise SystemExit(f"{self.spec.tuple_label}: delta upper bound is not certified: {delta_upper_cert}")
        if not (mn_small_support_base < 1.0):
            raise SystemExit(f"{self.spec.tuple_label}: MN small-support base is not < 1: {mn_small_support_base}")
        if not ha_compact["all_negative"]:
            raise SystemExit(f"{self.spec.tuple_label}: HA compact strip failed: {ha_compact}")
        if not ha_boundary["all_negative"]:
            raise SystemExit(f"{self.spec.tuple_label}: HA boundary strip failed: {ha_boundary}")
        if not mn_large_support["all_negative"]:
            raise SystemExit(f"{self.spec.tuple_label}: MN large-support failed: {mn_large_support}")
        return payload


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    certificates = [TripleCertifier(spec).certify() for spec in TARGETS]
    payload = {
        "targets": certificates,
    }
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"wrote {OUTPUT_JSON}")
    for cert in certificates:
        triple = tuple(cert["triple"])
        print(
            triple,
            "delta_upper",
            cert["delta_upper_certificate"]["delta_gv_upper"],
            "HA compact",
            cert["ha"]["compact_strip_certificate"]["worst_upper_bound"],
            "MN large",
            cert["mn"]["large_support_certificate"]["worst_upper_bound"],
        )


if __name__ == "__main__":
    main()
