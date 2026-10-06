# SPDX-License-Identifier: BSD-3-Clause
"""Generates the statsmodels reference values for tseries's SARIMAX golden tests.

Writes test/fixtures/sarimax_reference.dart. Run from the package root:

    python tool/sarimax_reference.py

Reference implementation: statsmodels SARIMAX, configured to be the SAME model
ctsa's sarimax fits -- regression with (seasonal) ARIMA errors, exact Gaussian
likelihood of the undifferenced series:

  * simple_differencing=False  (exact diffuse likelihood, not a fit to the
                                pre-differenced series);
  * trend='n'                  (no state intercept; an intercept, when wanted,
                                is passed as an explicit column of ones, which
                                is exactly ctsa's `imean` -- a regression mean,
                                not statsmodels' state-equation intercept);
  * standard errors from cov_type='approx': the inverse of a NUMERICAL Hessian
    of the exact log-likelihood at the optimum -- the same quantity ctsa's vcov
    computes (finite-difference Hessian of the sigma2-profiled likelihood; by
    the profile-likelihood identity the coefficient block is identical).
    NOT cov_type='oim': in statsmodels that is Harvey's (1989) analytic
    information formula, which on these fits differs from the numerical
    Hessian by 4-15% for the ARMA terms (checked 2026-10-04: ctsa agrees with
    'approx' and with an independent profile Hessian to ~1e-4, and disagrees
    with 'oim' exactly where 'approx' does).

Conventions to keep in mind when reading the fixture:
  * statsmodels' MA polynomial is 1 + theta*B; ctsa's (and tseries's) is
    1 - theta*B. The fixture stores statsmodels' numbers verbatim; the Dart test
    flips the sign.
  * The synthetic data below are generated once with a fixed seed and rounded to
    6 decimals; the rounded values are what both statsmodels and the Dart test
    fit, so the fixture is self-contained and exact.

Versions used to produce the committed fixture are printed into its header.
"""

import sys

import numpy as np
import scipy
import statsmodels
from statsmodels.tsa.statespace.sarimax import SARIMAX


def _arma_noise(rng, n, ar=(), ma=(), sar=(), sma=(), s=0, scale=1.0):
    """Simulates (seasonal) ARMA innovations-driven noise, burn-in discarded."""
    burn = 200 + 3 * max(s, 1)
    e = rng.normal(0.0, scale, n + burn)
    u = np.zeros(n + burn)
    for t in range(n + burn):
        v = e[t]
        for i, a in enumerate(ar, 1):
            if t - i >= 0:
                v += a * u[t - i]
        for i, m in enumerate(ma, 1):
            if t - i >= 0:
                v += m * e[t - i]
        for i, a in enumerate(sar, 1):
            if t - i * s >= 0:
                v += a * u[t - i * s]
        for i, m in enumerate(sma, 1):
            if t - i * s >= 0:
                v += m * e[t - i * s]
        u[t] = v
    return u[burn:]


def _integrate(u, d=0, D=0, s=0):
    x = u.copy()
    for _ in range(D):
        y = np.zeros_like(x)
        for t in range(len(x)):
            y[t] = x[t] + (y[t - s] if t - s >= 0 else 0.0)
        x = y
    for _ in range(d):
        x = np.cumsum(x)
    return x


def _r6(a):
    return np.round(np.asarray(a, dtype=float), 6)


def _fit(y, X, order, seasonal, h, Xf):
    model = SARIMAX(
        y,
        exog=X,
        order=order,
        seasonal_order=seasonal,
        trend="n",
        simple_differencing=False,
    )
    res = model.fit(
        disp=0, method="lbfgs", maxiter=5000, pgtol=1e-12, factr=10.0,
        cov_type="approx",
    )
    fc = res.get_forecast(h, exog=Xf)
    names = model.param_names
    return {
        "names": list(names),
        "params": [float(v) for v in res.params],
        "bse": [float(v) for v in res.bse],
        "llf": float(res.llf),
        "aic": float(res.aic),
        "forecast": [float(v) for v in fc.predicted_mean],
        "forecastSe": [float(v) for v in fc.se_mean],
    }


def build_cases():
    cases = {}

    # A -- ARIMA(1,1,1) errors + 2 regressors, n = 144, horizon 6.
    rng = np.random.default_rng(20261004)
    n, h = 144, 6
    t = np.arange(n + h)
    pulse = np.zeros(n + h)
    for start in range(5, n + h, 30):
        for k in range(24):
            if start + k < n + h:
                pulse[start + k] += 40.0 * (k / 6.0) * np.exp(1.0 - k / 6.0) / 6.0
    x1 = _r6(pulse)
    x2 = _r6(np.sin(2 * np.pi * t / 37.0) * 3.0 + rng.normal(0, 0.5, n + h))
    u = _integrate(_arma_noise(rng, n + h, ar=(0.5,), ma=(0.3,), scale=0.2), d=1)
    y = _r6(5.0 + u + 0.8 * x1 - 0.5 * x2)
    X = np.column_stack([x1, x2])
    cases["nonSeasonal"] = dict(
        y=y[:n], X=X[:n], Xf=X[n:], order=(1, 1, 1), seasonal=(0, 0, 0, 0),
        h=h, names=["x1", "x2"], includeMean=False,
    )

    # B -- airline-type SARIMA(0,1,1)(0,1,1)[12] errors + 1 regressor.
    rng = np.random.default_rng(12)
    n, h = 144, 12
    x = _r6(rng.normal(0, 1, n + h).cumsum() * 0.3 + rng.normal(0, 1, n + h))
    u = _integrate(
        _arma_noise(rng, n + h, ma=(-0.4,), sma=(-0.55,), s=12, scale=0.5),
        d=1, D=1, s=12,
    )
    y = _r6(100.0 + u + 1.5 * x)
    cases["seasonal"] = dict(
        y=y[:n], X=x[:n, None], Xf=x[n:, None], order=(0, 1, 1),
        seasonal=(0, 1, 1, 12), h=h, names=["x"], includeMean=False,
    )

    # C -- stationary ARMA(1,0,1) errors + intercept + 1 regressor (d = 0).
    rng = np.random.default_rng(7)
    n, h = 120, 6
    x = _r6(np.cos(np.arange(n + h) / 5.0) * 2.0 + rng.normal(0, 0.7, n + h))
    u = _arma_noise(rng, n + h, ar=(0.6,), ma=(0.25,), scale=1.0)
    y = _r6(10.0 + 2.0 * x + u)
    ones = np.ones(n + h)
    cases["stationaryWithMean"] = dict(
        y=y[:n], X=x[:n, None], Xf=x[n:, None], order=(1, 0, 1),
        seasonal=(0, 0, 0, 0), h=h, names=["x"], includeMean=True,
        smX=np.column_stack([ones[:n], x[:n]]),
        smXf=np.column_stack([ones[n:], x[n:]]),
    )

    # D -- ARIMA(0,1,1) with drift, the drift being the column t = 1..n.
    rng = np.random.default_rng(99)
    n, h = 100, 5
    u = _integrate(_arma_noise(rng, n + h, ma=(0.4,), scale=1.0), d=1)
    y = _r6(50.0 + 0.3 * np.arange(1, n + h + 1) + u)
    tt = np.arange(1, n + h + 1, dtype=float)
    cases["drift"] = dict(
        y=y[:n], X=np.zeros((n, 0)), Xf=np.zeros((h, 0)), order=(0, 1, 1),
        seasonal=(0, 0, 0, 0), h=h, names=[], includeMean=False,
        includeDrift=True, smX=tt[:n, None], smXf=tt[n:, None],
    )
    return cases


def _dart_list(values):
    return "<double>[" + ", ".join(repr(float(v)) for v in values) + "]"


def main(out_path):
    cases = build_cases()
    lines = [
        "// GENERATED by tool/sarimax_reference.py -- do not edit by hand.",
        "//",
        "// Reference: statsmodels SARIMAX (simple_differencing=False, trend='n',",
        "// cov_type='approx'). Produced with:",
        f"//   python      {sys.version.split()[0]}",
        f"//   numpy       {np.__version__}",
        f"//   scipy       {scipy.__version__}",
        f"//   statsmodels {statsmodels.__version__}",
        "//",
        "// MA coefficients are statsmodels' (1 + theta*B); tseries reports the",
        "// opposite sign. `params`/`bse` follow `names` (statsmodels' order).",
        "library;",
        "",
        "class SarimaxReferenceCase {",
        "  const SarimaxReferenceCase({",
        "    required this.y,",
        "    required this.exog,",
        "    required this.futureExog,",
        "    required this.order,",
        "    required this.seasonal,",
        "    required this.horizon,",
        "    required this.includeMean,",
        "    required this.includeDrift,",
        "    required this.names,",
        "    required this.params,",
        "    required this.bse,",
        "    required this.llf,",
        "    required this.aic,",
        "    required this.forecast,",
        "    required this.forecastSe,",
        "  });",
        "  final List<double> y;",
        "  final Map<String, List<double>> exog;",
        "  final Map<String, List<double>> futureExog;",
        "  /// (p, d, q)",
        "  final List<int> order;",
        "  /// (P, D, Q, s)",
        "  final List<int> seasonal;",
        "  final int horizon;",
        "  final bool includeMean;",
        "  final bool includeDrift;",
        "  final List<String> names;",
        "  final List<double> params;",
        "  final List<double> bse;",
        "  final double llf;",
        "  final double aic;",
        "  final List<double> forecast;",
        "  final List<double> forecastSe;",
        "",
        "  double param(String name) => params[names.indexOf(name)];",
        "  double se(String name) => bse[names.indexOf(name)];",
        "}",
        "",
    ]
    for key, c in cases.items():
        smX = c.get("smX", c["X"])
        smXf = c.get("smXf", c["Xf"])
        if smX.shape[1] == 0:
            smX, smXf = None, None
        ref = _fit(c["y"], smX, c["order"], c["seasonal"], c["h"], smXf)
        exog = ", ".join(
            f"'{nm}': {_dart_list(c['X'][:, j])}" for j, nm in enumerate(c["names"])
        )
        fexog = ", ".join(
            f"'{nm}': {_dart_list(c['Xf'][:, j])}" for j, nm in enumerate(c["names"])
        )
        print(key, dict(zip(ref["names"], np.round(ref["params"], 6))),
              "llf=%.6f" % ref["llf"], file=sys.stderr)
        lines += [
            f"const sarimaxRef{key[0].upper()}{key[1:]} = SarimaxReferenceCase(",
            f"  y: {_dart_list(c['y'])},",
            f"  exog: {{{exog}}},",
            f"  futureExog: {{{fexog}}},",
            f"  order: <int>[{', '.join(str(v) for v in c['order'])}],",
            f"  seasonal: <int>[{', '.join(str(v) for v in c['seasonal'])}],",
            f"  horizon: {c['h']},",
            f"  includeMean: {'true' if c['includeMean'] else 'false'},",
            f"  includeDrift: {'true' if c.get('includeDrift') else 'false'},",
            "  names: <String>[" + ", ".join(f"'{nm}'" for nm in ref["names"]) + "],",
            f"  params: {_dart_list(ref['params'])},",
            f"  bse: {_dart_list(ref['bse'])},",
            f"  llf: {ref['llf']!r},",
            f"  aic: {ref['aic']!r},",
            f"  forecast: {_dart_list(ref['forecast'])},",
            f"  forecastSe: {_dart_list(ref['forecastSe'])},",
            ");",
            "",
        ]
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "test/fixtures/sarimax_reference.dart")
