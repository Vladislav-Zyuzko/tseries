# SPDX-License-Identifier: BSD-3-Clause
"""Generates the statsmodels reference values for the AR(1) likelihood tests.

Writes test/fixtures/ar1_reference.dart. Run from the package root:

    python tool/ar1_likelihood_reference.py

Why this exists: until the fix recorded as "Local modifications" 10 in
third_party/ctsa/PROVENANCE.md, ctsa's AS 154 likelihood gave the first
observation of a pure AR(1) error process (ip == 1, iq == 0 after the seasonal
polynomials are multiplied out) variance sigma^2 instead of
sigma^2 / (1 - phi^2). The fixture pins, for every model form that reaches
that branch and for its neighbours that must not, the exact-MLE optimum.

Reference: statsmodels SARIMAX, configured as the SAME likelihood ctsa
maximises:

  * the series is differenced HERE (d, then D at lag s), and so are the
    regressors, and statsmodels fits an ARMA to the result -- ctsa's
    likelihood is that of the differenced series with a stationary start,
    which is statsmodels' simple_differencing=True, written out explicitly;
  * trend='n'; an intercept (d = D = 0, includeMean) is a column of ones and
    drift is the column t = 1..n before differencing -- a regression mean, as
    ctsa's `imean`/xreg, not statsmodels' state intercept;
  * concentrate_scale=True: sigma^2 profiled out, so `scale` is S/T and `llf`
    the profiled exact log-likelihood (exact pi, as ctsa's since PROVENANCE.md
    "Local modifications" 14; before that ctsa's was larger by
    T/2 * ln(pi / 3.14159) and the test subtracted it);
  * stationary initialisation (statsmodels' default for a stationary ARMA).

The optimum is the best of several starts and optimisers (see `_fit`), so it
is a reference for ctsa's optimum, not just one more optimiser run.

Synthetic series are generated once with a fixed seed and rounded to 6
decimals; the rounded values are what both statsmodels and the test fit. The
real series are written into the fixture verbatim, so the fixture is
self-contained.

MA coefficients are statsmodels' (1 + theta*B); tseries reports 1 - theta*B.
"""

import pathlib
import sys
import warnings

import numpy as np
import scipy
import statsmodels
from statsmodels.tsa.statespace.sarimax import SARIMAX

ROOT = pathlib.Path(__file__).resolve().parent.parent
SERIES = ROOT / "test" / "fixtures" / "auto_arima"


def _csv(name):
    lines = (SERIES / f"series_{name}.csv").read_text().split()
    return np.array([float(v) for v in lines if v.strip()])


def _r6(a):
    return np.round(np.asarray(a, dtype=float), 6)


def _ar1(rng, n, phi, mu):
    x = np.empty(n)
    x[0] = rng.normal() / np.sqrt(1.0 - phi * phi)
    for t in range(1, n):
        x[t] = phi * x[t - 1] + rng.normal()
    return _r6(mu + x)


def _difference(a, d, D, s):
    a = np.asarray(a, dtype=float)
    for _ in range(d):
        a = a[1:] - a[:-1]
    for _ in range(D):
        a = a[s:] - a[:-s]
    return a


def _fit(w, X, order, seasonal):
    """Best of several statsmodels optimisations of the same likelihood."""
    mod = SARIMAX(
        w, exog=X, order=order, seasonal_order=seasonal, trend="n",
        concentrate_scale=True,
    )
    best = None
    starts = [None]
    k = mod.k_params
    # Individual optimisers may stop on their iteration caps; only the best
    # result matters, so their convergence warnings are noise here.
    warnings.simplefilter("ignore")
    for start in starts:
        for method, kw in (
            ("lbfgs", dict(maxiter=10000, pgtol=1e-12, factr=1.0)),
            ("bfgs", dict(maxiter=10000, gtol=1e-10)),
            ("nm", dict(maxiter=40000, xtol=1e-12, ftol=1e-14)),
        ):
            res = mod.fit(start_params=start, method=method, disp=0, **kw)
            if best is None or res.llf > best.llf:
                best = res
    # Polish from the best point found.
    for method, kw in (
        ("bfgs", dict(maxiter=10000, gtol=1e-11)),
        ("nm", dict(maxiter=40000, xtol=1e-13, ftol=1e-15)),
    ):
        res = mod.fit(start_params=best.params, method=method, disp=0, **kw)
        if res.llf > best.llf:
            best = res
    assert len(best.params) == k
    return mod, best


def build_cases():
    nile = _csv("nile")
    lynx = _csv("lynx")
    huron = _csv("lake_huron")
    www = _csv("wwwusage")
    usacc = _csv("us_acc_deaths")
    airlog = _csv("air_passengers_log")
    reg_nile = _r6(np.sin(np.arange(len(nile)) / 5.0) * 3.0)

    def case(name, y, order, seasonal=(0, 0, 0, 0), mean=False, drift=False,
             exog=None, in_branch=True, sm_order=None, sm_seasonal=None):
        return dict(
            name=name, y=np.asarray(y, dtype=float), order=order,
            seasonal=seasonal, includeMean=mean, includeDrift=drift,
            exog=exog or {}, inBranch=in_branch,
            smOrder=sm_order or order, smSeasonal=sm_seasonal or seasonal,
        )

    cases = [
        # --- reach the ip == 1 && iq == 0 branch ---------------------------
        case("nile_100_mean", nile, (1, 0, 0), mean=True),
        case("nile_100_nomean", nile, (1, 0, 0)),
        case("nile_100_mean_reg", nile, (1, 0, 0), mean=True,
             exog={"x": reg_nile}),
        case("lynx_100_mean", lynx, (1, 0, 0), mean=True),
        case("huron_100_mean", huron, (1, 0, 0), mean=True),
        case("www_110", www, (1, 1, 0)),
        case("www_110_drift", www, (1, 1, 0), drift=True),
        case("usacc_100_010_12", usacc, (1, 0, 0), (0, 1, 0, 12)),
        case("usacc_100_010_12_drift", usacc, (1, 0, 0), (0, 1, 0, 12),
             drift=True),
        case("usacc_110_010_12", usacc, (1, 1, 0), (0, 1, 0, 12)),
        case("airlog_110_010_12", airlog, (1, 1, 0), (0, 1, 0, 12)),
        # A seasonal AR(1) at period 1 multiplies out to ip == 1 as well.
        case("nile_000_100_1_mean", nile, (0, 0, 0), (1, 0, 0, 1), mean=True,
             sm_order=(1, 0, 0), sm_seasonal=(0, 0, 0, 0)),
        # --- neighbours that must NOT reach it ------------------------------
        case("usacc_100_100_12_mean", usacc, (1, 0, 0), (1, 0, 0, 12),
             mean=True, in_branch=False),
        case("airlog_110_011_12", airlog, (1, 1, 0), (0, 1, 1, 12),
             in_branch=False),
        case("nile_001_mean", nile, (0, 0, 1), mean=True, in_branch=False),
        case("nile_001_nomean", nile, (0, 0, 1), in_branch=False),
        case("huron_011", huron, (0, 1, 1), in_branch=False),
        case("lynx_200_mean", lynx, (2, 0, 0), mean=True, in_branch=False),
    ]
    rng = np.random.default_rng(20261005)
    for phi in (-0.8, 0.05, 0.98, 0.995):
        for n in (30, 200):
            tag = f"{'m' if phi < 0 else ''}{abs(phi):g}".replace(".", "")
            cases.append(case(f"syn_phi{tag}_n{n}_mean",
                              _ar1(rng, n, phi, 10.0), (1, 0, 0), mean=True))
            cases.append(case(f"syn_phi{tag}_n{n}_nomean",
                              _ar1(rng, n, phi, 0.0), (1, 0, 0)))
    return cases


def reference(c):
    p, d, q = c["order"]
    P, D, Q, s = c["seasonal"]
    n = len(c["y"])
    cols = []
    names = []
    if c["includeMean"] and d + D == 0:
        cols.append(np.ones(n))
        names.append("intercept")
    for nm, v in c["exog"].items():
        cols.append(np.asarray(v, dtype=float))
        names.append(nm)
    if c["includeDrift"]:
        cols.append(np.arange(1, n + 1, dtype=float))
        names.append("drift")
    w = _difference(c["y"], d, D, s)
    X = None
    if cols:
        X = np.column_stack([_difference(col, d, D, s) for col in cols])
    sp, _, sq = c["smOrder"]
    sP, _, sQ, ss = c["smSeasonal"]
    mod, res = _fit(w, X, (sp, 0, sq), (sP, 0, sQ, ss))
    # statsmodels lists the regression coefficients first; name them as
    # tseries does (statsmodels would call a differenced drift column 'const').
    sm_names = list(mod.param_names)
    out_names = names + sm_names[len(names):]
    if c["smOrder"] != c["order"]:
        # The period-1 seasonal AR fitted as an ordinary AR(1).
        out_names = [{"ar.L1": "ar.S.L1"}.get(nm, nm) for nm in out_names]
    return dict(
        T=len(w), llf=float(res.llf), scale=float(res.scale),
        names=out_names, params=[float(v) for v in res.params],
    )


def _dart_list(values):
    return "<double>[" + ", ".join(repr(float(v)) for v in values) + "]"


def main(out_path):
    lines = [
        "// GENERATED by tool/ar1_likelihood_reference.py -- do not edit by hand.",
        "//",
        "// Reference: statsmodels SARIMAX on the explicitly differenced series",
        "// (trend='n', concentrate_scale=True, stationary initialisation); best",
        "// of several optimisers. Produced with:",
        f"//   python      {sys.version.split()[0]}",
        f"//   numpy       {np.__version__}",
        f"//   scipy       {scipy.__version__}",
        f"//   statsmodels {statsmodels.__version__}",
        "//",
        "// MA coefficients are statsmodels' (1 + theta*B); tseries reports the",
        "// opposite sign. `params` follow `names` (statsmodels' order, with",
        "// tseries' names: intercept, the exog keys, drift, ar.L1, ma.L1,",
        "// ar.S.L<s>, ma.S.L<s>).",
        "library;",
        "",
        "class Ar1ReferenceCase {",
        "  const Ar1ReferenceCase({",
        "    required this.name,",
        "    required this.y,",
        "    required this.exog,",
        "    required this.order,",
        "    required this.seasonal,",
        "    required this.includeMean,",
        "    required this.includeDrift,",
        "    required this.inBranch,",
        "    required this.t,",
        "    required this.llf,",
        "    required this.scale,",
        "    required this.names,",
        "    required this.params,",
        "  });",
        "  final String name;",
        "  final List<double> y;",
        "  final Map<String, List<double>> exog;",
        "",
        "  /// (p, d, q) as fitted by tseries.",
        "  final List<int> order;",
        "",
        "  /// (P, D, Q, s) as fitted by tseries.",
        "  final List<int> seasonal;",
        "  final bool includeMean;",
        "  final bool includeDrift;",
        "",
        "  /// Whether the model reaches ctsa's ip == 1 && iq == 0 branch (the",
        "  /// one the AR(1) fix changed).",
        "  final bool inBranch;",
        "",
        "  /// Observations in the likelihood (after differencing).",
        "  final int t;",
        "",
        "  /// Profiled exact log-likelihood at the optimum (exact pi).",
        "  final double llf;",
        "",
        "  /// sigma^2 at the optimum, S/T.",
        "  final double scale;",
        "  final List<String> names;",
        "  final List<double> params;",
        "",
        "  double param(String name) => params[names.indexOf(name)];",
        "}",
        "",
        "const ar1ReferenceCases = <Ar1ReferenceCase>[",
    ]
    for c in build_cases():
        ref = reference(c)
        print(c["name"], dict(zip(ref["names"], np.round(ref["params"], 6))),
              "llf=%.6f" % ref["llf"], file=sys.stderr)
        exog = ", ".join(f"'{k}': {_dart_list(v)}" for k, v in c["exog"].items())
        lines += [
            "  Ar1ReferenceCase(",
            f"    name: '{c['name']}',",
            f"    y: {_dart_list(c['y'])},",
            f"    exog: {{{exog}}},",
            f"    order: <int>[{', '.join(str(v) for v in c['order'])}],",
            f"    seasonal: <int>[{', '.join(str(v) for v in c['seasonal'])}],",
            f"    includeMean: {'true' if c['includeMean'] else 'false'},",
            f"    includeDrift: {'true' if c['includeDrift'] else 'false'},",
            f"    inBranch: {'true' if c['inBranch'] else 'false'},",
            f"    t: {ref['T']},",
            f"    llf: {ref['llf']!r},",
            f"    scale: {ref['scale']!r},",
            "    names: <String>[" + ", ".join(f"'{nm}'" for nm in ref["names"]) + "],",
            f"    params: {_dart_list(ref['params'])},",
            "  ),",
        ]
    lines += ["];", ""]
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "test/fixtures/ar1_reference.dart")
