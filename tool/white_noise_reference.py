# SPDX-License-Identifier: BSD-3-Clause
"""Generates the statsmodels reference values for tseries's white-noise tests.

Writes test/fixtures/white_noise_reference.dart. Run from the package root:

    python tool/white_noise_reference.py

Covers the "white-noise" orders (p = q = P = Q = 0): white noise with a mean,
the random walk (0,1,0) with and without drift, (0,2,0), the seasonal random
walk (0,0,0)(0,1,0)[12] and (0,1,0)(0,1,0)[12]. Before tseries 0.8.x these
returned uninitialised forecast standard errors (and, for d + D > 0, point
forecasts that were whatever the caller's buffer held); see PROVENANCE.md,
Local modifications 12.

Reference implementation: statsmodels SARIMAX, the same model tseries fits --
regression with (seasonal) ARIMA errors, exact likelihood of the undifferenced
series:

  * trend='n'; the mean (d = D = 0) and the drift are explicit regression
    columns (a column of ones; t = 1..n), exactly as tseries passes them;
  * use_exact_diffuse=True, so the d + s*D starting values are conditioned on
    exactly (the default approximate diffuse prior leaves ~1e-6 noise in the
    forecast variances).

For these orders the maximum-likelihood estimates have a closed form (the mean
or drift is the sample mean of the differenced series, sigma^2 the mean squared
deviation from it), so the optimiser's answer is replaced by the closed form
and statsmodels is used to evaluate the log-likelihood, forecasts and forecast
standard errors AT those parameters (`smooth`). The fixture also records the
closed-form parameters themselves.

Versions used to produce the committed fixture are printed into its header.
"""

import sys

import numpy as np
import scipy
import statsmodels
from statsmodels.tsa.statespace.sarimax import SARIMAX

H = 24


def _r6(a):
    return np.round(np.asarray(a, dtype=float), 6)


def _difference(y, d, D, s):
    w = np.asarray(y, dtype=float)
    for _ in range(D):
        w = w[s:] - w[:-s]
    for _ in range(d):
        w = w[1:] - w[:-1]
    return w


def build_cases():
    cases = {}

    rng = np.random.default_rng(2026_10_05)
    cases["whiteNoiseMean"] = dict(
        y=_r6(5.0 + rng.normal(0, 0.8, 50)), d=0, D=0, s=0, mean=True,
        drift=False)

    rng = np.random.default_rng(11)
    cases["randomWalk"] = dict(
        y=_r6(100.0 + rng.normal(0, 1.5, 60).cumsum()), d=1, D=0, s=0,
        mean=False, drift=False)

    rng = np.random.default_rng(12)
    cases["randomWalkDrift"] = dict(
        y=_r6(20.0 + (0.4 + rng.normal(0, 1.0, 80)).cumsum()), d=1, D=0, s=0,
        mean=False, drift=True)

    rng = np.random.default_rng(13)
    cases["integratedTwice"] = dict(
        y=_r6(rng.normal(0, 0.3, 45).cumsum().cumsum() + 10.0), d=2, D=0,
        s=0, mean=False, drift=False)

    rng = np.random.default_rng(14)
    t = np.arange(72)
    seas = 3.0 * np.sin(2 * np.pi * t / 12.0)
    base = np.zeros(72)
    e = rng.normal(0, 0.5, 72)
    for i in range(72):
        base[i] = (base[i - 12] if i >= 12 else 50.0 + seas[i]) + e[i]
    cases["seasonalRandomWalk"] = dict(
        y=_r6(base), d=0, D=1, s=12, mean=False, drift=False)

    rng = np.random.default_rng(15)
    w = rng.normal(0, 0.7, 96)
    y = np.zeros(96)
    for i in range(96):
        if i < 13:
            y[i] = 200.0 + 5.0 * np.cos(2 * np.pi * i / 12.0) + w[i]
        else:  # (1 - B)(1 - B^12) y = w
            y[i] = y[i - 1] + y[i - 12] - y[i - 13] + w[i]
    cases["airlineIntegration"] = dict(
        y=_r6(y), d=1, D=1, s=12, mean=False, drift=False)
    return cases


def _reference(c):
    y = np.asarray(c["y"], dtype=float)
    n = len(y)
    d, D, s = c["d"], c["D"], c["s"]
    X = None
    names = []
    if c["mean"]:
        X = np.ones((n, 1))
        names = ["intercept"]
    if c["drift"]:
        X = np.arange(1, n + 1, dtype=float)[:, None]
        names = ["drift"]
    w = _difference(y, d, D, s)
    beta = []
    if X is not None:
        wx = _difference(X[:, 0], d, D, s)
        # Closed-form least squares on the differenced data (one column).
        b = float(np.dot(wx, w) / np.dot(wx, wx))
        beta = [b]
        resid = w - b * wx
    else:
        resid = w
    sigma2 = float(np.mean(resid ** 2))

    model = SARIMAX(
        y, exog=X, order=(0, d, 0),
        seasonal_order=(0, D, 0, s) if D > 0 else (0, 0, 0, 0),
        trend="n", use_exact_diffuse=True)
    params = np.array(beta + [sigma2])
    res = model.smooth(params)
    if X is not None:
        Xf = (np.ones((H, 1)) if c["mean"]
              else np.arange(n + 1, n + H + 1, dtype=float)[:, None])
    else:
        Xf = None
    fc = res.get_forecast(H, exog=Xf)
    return dict(
        names=names, beta=beta, sigma2=sigma2, llf=float(res.llf),
        forecast=np.asarray(fc.predicted_mean, dtype=float),
        stderr=np.asarray(fc.se_mean, dtype=float))


def _dart_list(values, kind="double"):
    return f"<{kind}>[" + ", ".join(repr(float(v)) if kind == "double"
                                     else str(v) for v in values) + "]"


def main(out_path):
    cases = build_cases()
    lines = [
        "// GENERATED by tool/white_noise_reference.py -- do not edit by hand.",
        "//",
        "// Reference: statsmodels SARIMAX (trend='n', use_exact_diffuse=True)",
        "// evaluated at the closed-form MLE of each white-noise order. Produced",
        "// with:",
        f"//   python      {sys.version.split()[0]}",
        f"//   numpy       {np.__version__}",
        f"//   scipy       {scipy.__version__}",
        f"//   statsmodels {statsmodels.__version__}",
        "library;",
        "",
        "class WhiteNoiseReferenceCase {",
        "  const WhiteNoiseReferenceCase({",
        "    required this.y,",
        "    required this.d,",
        "    required this.seasonalD,",
        "    required this.period,",
        "    required this.includeMean,",
        "    required this.includeDrift,",
        "    required this.beta,",
        "    required this.sigma2,",
        "    required this.llf,",
        "    required this.forecast,",
        "    required this.stderr,",
        "  });",
        "  final List<double> y;",
        "  final int d;",
        "  final int seasonalD;",
        "  final int period;",
        "  final bool includeMean;",
        "  final bool includeDrift;",
        "",
        "  /// Closed-form MLE of the mean or drift coefficient (empty if none).",
        "  final List<double> beta;",
        "",
        "  /// Closed-form MLE of sigma^2 (mean squared one-step residual).",
        "  final double sigma2;",
        "",
        "  /// statsmodels' exact log-likelihood at (beta, sigma2).",
        "  final double llf;",
        "",
        "  /// statsmodels' forecasts and their standard errors, h = 1..24.",
        "  final List<double> forecast;",
        "  final List<double> stderr;",
        "}",
        "",
        "const whiteNoiseReference = <String, WhiteNoiseReferenceCase>{",
    ]
    for name, c in cases.items():
        r = _reference(c)
        lines += [
            f"  '{name}': WhiteNoiseReferenceCase(",
            f"    y: {_dart_list(c['y'])},",
            f"    d: {c['d']},",
            f"    seasonalD: {c['D']},",
            f"    period: {c['s']},",
            f"    includeMean: {'true' if c['mean'] else 'false'},",
            f"    includeDrift: {'true' if c['drift'] else 'false'},",
            f"    beta: {_dart_list(r['beta'])},",
            f"    sigma2: {r['sigma2']!r},",
            f"    llf: {r['llf']!r},",
            f"    forecast: {_dart_list(r['forecast'])},",
            f"    stderr: {_dart_list(r['stderr'])},",
            "  ),",
        ]
    lines += ["};", ""]
    with open(out_path, "w", newline="\n") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1
         else "test/fixtures/white_noise_reference.dart")
