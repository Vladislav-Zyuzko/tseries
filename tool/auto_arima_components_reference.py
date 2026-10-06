#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, Vladislav Zyuzko
"""Independent reference numbers for the pure-Dart components of autoArima
(gate A of the autoArima specification, §9).

statsmodels and numpy are used as EXTERNAL TOOLS: only their numeric output
is written to the fixtures; no code of theirs is part of the package.

Writes, under test/fixtures/auto_arima_components/:

  kpss_reference.txt     name|transform|T|lag|eta
      statsmodels.tsa.stattools.kpss(x, regression='c', nlags=lag)[0] with
      lag = floor(4*(T/100)**0.25) (the KPSS 1992 "l4" rule), on the six
      series of test/fixtures/auto_arima/, their first and second differences
      and, for the monthly series, the seasonal difference and its first
      difference.

  decomposition_reference.txt   name|period|F_S|i_0,...,i_{s-1}
      statsmodels.tsa.seasonal.seasonal_decompose(x, model='additive',
      period=s) (centred moving average, 2 x s for even s, seasonal indices
      centred to zero); F_S = max(0, 1 - Var(R)/Var(S+R)) over the positions
      where the trend is defined. i_j = the seasonal component at position j.

  roots_reference.txt    a_1,...,a_m|min_modulus_mpmath|min_modulus_numpy
      10^4 random polynomials a(z) = 1 - sum a_j z^j (5000 of degree <= 5,
      5000 of degree <= 2; a third with random coefficients, a third with
      root moduli spread over [0.9, 1.2], a third with root moduli within
      1e-5 of 1.01, the default margin), coefficients rounded to 9 decimals, and the smallest
      root modulus of the rounded polynomial from numpy.roots and from
      mpmath.polyroots at 60 significant digits. numpy is the reference the
      spec names; on clustered roots near the boundary it is itself only
      ~1e-8 accurate, so the 60-digit mpmath value is the arbiter.

Usage (from packages/tseries):
    python tool/auto_arima_components_reference.py
Tested with statsmodels 0.14.6 / numpy 2.5.3 / mpmath 1.4.1 / Python 3.12.
"""

import math
import pathlib
import warnings

import mpmath
import numpy as np
import statsmodels
from statsmodels.tsa.seasonal import seasonal_decompose
from statsmodels.tsa.stattools import kpss

ROOT = pathlib.Path(__file__).resolve().parent.parent
SERIES = ROOT / "test" / "fixtures" / "auto_arima"
OUT = ROOT / "test" / "fixtures" / "auto_arima_components"

NAMES = [
    ("air_passengers_log", 12),
    ("wwwusage", 1),
    ("lynx", 1),
    ("lake_huron", 1),
    ("nile", 1),
    ("us_acc_deaths", 12),
]


def load(name):
    return np.array(
        [float(l) for l in (SERIES / f"series_{name}.csv").read_text().split()]
    )


def l4(t):
    return int(math.floor(4 * (t / 100) ** 0.25))


def kpss_lines():
    lines = []
    for name, s in NAMES:
        y = load(name)
        transforms = [("level", y), ("diff1", np.diff(y)), ("diff2", np.diff(y, 2))]
        if s > 1:
            sd = y[s:] - y[:-s]
            transforms += [("sdiff", sd), ("sdiff_diff1", np.diff(sd))]
        for label, x in transforms:
            lag = l4(len(x))
            with warnings.catch_warnings():
                # statsmodels warns when eta is outside its p-value table; the
                # statistic itself is what we take.
                warnings.simplefilter("ignore")
                eta = kpss(x, regression="c", nlags=lag)[0]
            lines.append(f"{name}|{label}|{len(x)}|{lag}|{float(eta)!r}")
    return lines


def fs(seasonal, resid):
    mask = ~np.isnan(resid)
    r = resid[mask]
    sr = seasonal[mask] + r
    return max(0.0, 1 - np.var(r, ddof=1) / np.var(sr, ddof=1))


def decomposition_lines():
    lines = []
    cases = [("air_passengers_log", 12), ("us_acc_deaths", 12), ("nile", 5), ("lynx", 10)]
    for name, s in cases:
        y = load(name)
        d = seasonal_decompose(y, model="additive", period=s)
        idx = ",".join(repr(float(v)) for v in d.seasonal[:s])
        lines.append(f"{name}|{s}|{float(fs(d.seasonal, d.resid))!r}|{idx}")
    return lines


def poly_from_roots(roots):
    # prod (1 - z/r) has constant term 1; return a_j with 1 - sum a_j z^j.
    c = np.array([1.0 + 0j])
    for r in roots:
        c = np.convolve(c, np.array([1.0, -1.0 / r]))
    c = np.real(c)
    return [-v for v in c[1:]]


def random_roots(rng, m, lo, hi):
    roots = []
    while len(roots) < m:
        mod = rng.uniform(lo, hi)
        if m - len(roots) >= 2 and rng.random() < 0.6:
            ang = rng.uniform(0.05, math.pi - 0.05)
            r = mod * complex(math.cos(ang), math.sin(ang))
            roots += [r, r.conjugate()]
        else:
            roots.append(mod * (1 if rng.random() < 0.5 else -1))
    return roots


def roots_lines():
    rng = np.random.default_rng(20261005)
    lines = []
    for max_deg, count in [(5, 5000), (2, 5000)]:
        for i in range(count):
            m = int(rng.integers(1, max_deg + 1))
            kind = i % 3
            if kind == 0:
                a = list(rng.uniform(-1.5, 1.5, size=m))
            elif kind == 1:
                a = poly_from_roots(random_roots(rng, m, 0.9, 1.2))
            else:
                a = poly_from_roots(random_roots(rng, m, 1.01 - 1e-5, 1.01 + 1e-5))
            a = [round(float(v), 9) for v in a]
            while a and a[-1] == 0.0:
                a.pop()
            if not a:
                continue
            coeffs = [-v for v in reversed(a)] + [1.0]
            mod = float(np.min(np.abs(np.roots(coeffs))))
            with mpmath.workdps(60):
                exact = [mpmath.mpf(repr(v)) for v in coeffs]
                rts = mpmath.polyroots(exact, maxsteps=2000, extraprec=600)
                hp = float(min(abs(x) for x in rts))
            lines.append(",".join(repr(v) for v in a) + f"|{hp!r}|{mod!r}")
    return lines


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    head = (
        f"# statsmodels {statsmodels.__version__} | numpy {np.__version__} | "
        f"mpmath {mpmath.__version__} | "
        "tool/auto_arima_components_reference.py\n"
    )
    (OUT / "kpss_reference.txt").write_text(
        head + "# name|transform|T|lag|eta\n" + "\n".join(kpss_lines()) + "\n"
    )
    (OUT / "decomposition_reference.txt").write_text(
        head + "# name|period|F_S|seasonal indices\n"
        + "\n".join(decomposition_lines()) + "\n"
    )
    (OUT / "roots_reference.txt").write_text(
        head + "# a_1,...,a_m|min root modulus of 1 - sum a_j z^j: "
        "mpmath (60 digits)|numpy.roots\n"
        + "\n".join(roots_lines()) + "\n"
    )


if __name__ == "__main__":
    main()
