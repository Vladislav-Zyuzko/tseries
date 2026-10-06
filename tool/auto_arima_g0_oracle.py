# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, Vladislav Zyuzko
"""Gate G0 of the autoArima specification (§9): every model form fitted by
autoArima has its log-likelihood re-evaluated by statsmodels at the
parameters tseries returned.

Input: the G0 export of tool/auto_arima_journal.dart (`--g0 <file>`), one line
per fitted candidate:
    series|mode|index|p|d|q|P|D|Q|s|c|verdict|status|loglik|sigma2|ar|ma|sar|sma|intercept|drift
MA/SMA in ctsa's 1 - theta*B convention; c = mean when d + D = 0, drift when
d + D = 1; sigma2 = profile MLE S/T.

How the reference is built (the same model tseries fits):
* fitSarimax is a regression with (S)ARIMA errors on the DIFFERENCED series:
  ctsa differences y (and any regressor) and evaluates the exact Gaussian
  likelihood of the differenced data, with sigma2 concentrated out (S/T,
  T = n - d - s*D). statsmodels: SARIMAX(..., simple_differencing=True,
  concentrate_scale=True), whose `loglike` is the same concentrated exact
  likelihood of the differenced series.
* Mean (d + D = 0): a regression on a column of ones (y - mu is ARMA). NOT
  statsmodels' trend='c', which is an intercept of the ARMA recursion
  (mean = c / (phi(1) * Phi(1))); the two are equivalent only after that
  conversion, which the script also checks (column `ll_sm_trend`).
* Drift (d + D = 1): fitSarimax adds the regressor t = 1..n and differences it
  with the series, i.e. the differenced model has the constant beta (and
  s*beta for a seasonal difference). statsmodels with simple_differencing
  differences exog too; the exog column t is passed and differenced by it.
* Signs: statsmodels MA = -ctsa MA (both regular and seasonal).
* Tolerance: |ll_ctsa - ll_sm| <= 1e-6 * max(1, |ll_sm|).
Fits with status 10/12 (CSS-only results) are not part of G0, nor are fits
without a reported log-likelihood (status != 1, e.g. 4); other rejected fits
are evaluated for information only.

Usage:
    python tool/auto_arima_g0_oracle.py <g0.txt> <out.tsv>
(run from packages/tseries; needs numpy and statsmodels).
"""
import math
import sys
import warnings

import numpy as np
from statsmodels.tsa.statespace.sarimax import SARIMAX

warnings.simplefilter("ignore")
TOL = 1e-6


def series(name):
    with open(f"test/fixtures/auto_arima/series_{name}.csv") as f:
        return np.array([float(l) for l in f if l.strip()])


def floats(field):
    return [float(v) for v in field.split(",")] if field and field != "-" else []


def sm_loglike(y, p, d, q, P, D, Q, s, ar, ma, sar, sma, mean, drift, trend_form):
    n = len(y)
    exog = None
    beta = []
    trend = "n"
    if mean is not None:
        if trend_form:
            trend = "c"
        else:
            exog = np.ones((n, 1))
            beta = [mean]
    if drift is not None:
        exog = np.arange(1, n + 1, dtype=float).reshape(-1, 1)
        beta = [drift]
    seasonal = (P, D, Q, s) if (P + D + Q) > 0 else (0, 0, 0, 0)
    mod = SARIMAX(
        y, exog=exog, order=(p, d, q), seasonal_order=seasonal, trend=trend,
        simple_differencing=True, concentrate_scale=True,
        # enforce_stationarity=True only selects the exact stationary
        # initialisation (with False statsmodels switches to an approximate
        # diffuse prior and the likelihood is no longer the exact one ctsa
        # computes); loglike() takes the params as given (transformed=True).
        enforce_stationarity=True, enforce_invertibility=False,
    )
    params = []
    if trend == "c":
        phi1 = 1.0 - sum(ar)
        Phi1 = 1.0 - sum(sar)
        params.append(mean * phi1 * Phi1)
    params += beta + ar + [-v for v in ma] + sar + [-v for v in sma]
    assert len(params) == len(mod.param_names), (mod.param_names, params)
    return float(mod.loglike(np.array(params))), mod


def main():
    rows = []
    cache = {}
    for line in open(sys.argv[1], encoding="utf-8"):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("|")
        name, mode, idx = f[0], f[1], int(f[2])
        p, d, q, P, D, Q, s, c = map(int, f[3:11])
        verdict, status = f[11], int(f[12])
        if f[13] == "-":
            # No log-likelihood reported (ctsa status != 1: the shim withholds it).
            rows.append((name, mode, idx, f"({p},{d},{q})({P},{D},{Q})[{s}]", c, verdict, status, None, None, None, None, f"excluded (no logLik, status {status})"))
            continue
        ll = float(f[13])
        ar, ma, sar, sma = floats(f[15]), floats(f[16]), floats(f[17]), floats(f[18])
        mean = None if f[19] == "-" else float(f[19])
        drift = None if f[20] == "-" else float(f[20])
        if status in (10, 12):
            rows.append((name, mode, idx, f"({p},{d},{q})({P},{D},{Q})[{s}]", c, verdict, status, ll, None, None, None, "excluded (CSS-only)"))
            continue
        if name not in cache:
            cache[name] = series(name)
        y = cache[name]
        llsm, _ = sm_loglike(y, p, d, q, P, D, Q, s, ar, ma, sar, sma, mean, drift, False)
        lltr = None
        if mean is not None:
            lltr, _ = sm_loglike(y, p, d, q, P, D, Q, s, ar, ma, sar, sma, mean, drift, True)
        rel = abs(ll - llsm) / max(1.0, abs(llsm))
        scope = "G0" if verdict == "accepted" else "info"
        ok = "PASS" if rel <= TOL else "FAIL"
        rows.append((name, mode, idx, f"({p},{d},{q})({P},{D},{Q})[{s}]", c, verdict, status, ll, llsm, lltr, rel, f"{scope} {ok}"))
    with open(sys.argv[2], "w", encoding="utf-8") as out:
        out.write("series\tmode\tindex\torder\tc\tverdict\tstatus\tll_ctsa\tll_sm\tll_sm_trend\trel\tresult\n")
        for r in rows:
            out.write("\t".join("" if v is None else (f"{v:.17g}" if isinstance(v, float) else str(v)) for v in r) + "\n")
    g0 = [r for r in rows if r[11].startswith("G0")]
    info = [r for r in rows if r[11].startswith("info")]
    worst = max(g0, key=lambda r: r[10])
    print(f"G0 (accepted): {len(g0)} checked, {sum(r[11].endswith('PASS') for r in g0)} pass, "
          f"worst rel {worst[10]:.3e} ({worst[0]} {worst[1]} {worst[3]} c={worst[4]})")
    if info:
        wi = max(info, key=lambda r: r[10])
        print(f"info (rejected, not CSS-only): {len(info)} checked, {sum(r[11].endswith('PASS') for r in info)} within tol, worst rel {wi[10]:.3e} ({wi[0]} {wi[3]} {wi[5]})")
    tr = [abs(r[9] - r[8]) / max(1, abs(r[8])) for r in rows if r[9] is not None]
    if tr:
        print(f"mean as trend='c' (converted) vs regression on ones: max rel diff {max(tr):.3e} over {len(tr)} fits")
    print("excluded (no logLik or status 10/12):", sum(1 for r in rows if r[11].startswith("excluded")))
    for r in rows:
        if r[11].endswith("FAIL"):
            print("FAIL", r)


if __name__ == "__main__":
    main()
