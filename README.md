# tseries

A framework-agnostic **time-series analysis SDK for Dart**: ARIMA, SARIMA and
SARIMAX powered by the [ctsa](https://github.com/rafat/ctsa) C library, plus
pure-Dart forecasting primitives.

`tseries` provides a clean, idiomatic Dart API for time-series work. It is a
plain Dart package with **no Flutter dependency**, so it is reusable from
servers, CLIs and tests, and independently testable.

> Status: pure-Dart primitives are live, and the **native ARIMA/SARIMA core
> (ctsa, over FFI) is live and validated** against published reference fits
> (see "Numerical validation" below). Verified platforms are listed under
> [Requirements](#requirements).

## Requirements

- **Dart ≥ 3.10 / Flutter ≥ 3.38.** The native core is compiled by a Dart
  [build hook](https://dart.dev/tools/hooks) (`hook/build.dart`), and build
  hooks / code assets are stable from Dart 3.10 (Flutter 3.38). No
  `--enable-experiment` flag is needed.
- **A C toolchain on the machine that builds the consumer app.** The C code
  is compiled from source on first build — nothing prebuilt is downloaded:
  - Android: the Android NDK (found automatically from the Android SDK);
  - iOS / macOS: Xcode (command-line tools);
  - Windows: MSVC (Visual Studio or Build Tools with the "Desktop
    development with C++" workload) — `native_toolchain_c` looks for
    `cl.exe`, not clang;
  - Linux: clang (`native_toolchain_c` uses gcc only for cross-compiling).
- **Web is not supported**: the core is reached through `dart:ffi`, which
  does not exist on the web.

| Platform | Status |
| --- | --- |
| Windows x64 | Verified: full test suite (MSVC, Dart 3.10.9), on every change. |
| Android | Verified on an API 36 emulator (no physical device): arm64-v8a and x86_64 built by the hook (NDK 29, Flutter 3.38.10, `-ffp-contract=off`); x86_64 run natively, arm64-v8a run under the emulator's ARM translation. checkNativeCore and a SARIMAX fit give the same numbers as on the Windows host on both ABIs; autoArima picks the same (d, D) and the same model as on Windows in all 166 test runs on both ABIs (see **Determinism** under [autoArima](#automatic-order-selection-autoarima)). armeabi-v7a builds and links (no FMA instructions) but has not been run. |
| Linux x64 | Verified: full test suite in Docker (`dart:3.10.9`, Debian, clang 19); `libctsa.so` links `libm.so.6`, exports only the `tseries_*` shim; autoArima: same (d, D) and model as Windows in all 166 test runs. Not tested on a desktop Flutter app. |
| iOS, macOS | Build flags are in place (`-fvisibility=hidden`, `-dead_strip`, libSystem math) but **not verified** — no macOS build host yet. |

## What lives where (by design)

Simple, well-understood arithmetic is implemented in **pure Dart**. Only
genuinely hard, validated numerics are delegated to a **native C core** across
an FFI boundary. Adding a native dependency has to earn its keep.

| Capability | Implementation |
| --- | --- |
| Simple moving average (SMA) | Pure Dart — available |
| Linear (OLS) trend forecast + prediction SE | Pure Dart — available |
| Damped Holt exponential smoothing (caller-supplied alpha/beta/phi) | `holtForecast` — pure Dart — available |
| SES / EWMA | Pure Dart — available as `holtForecast` with `phi == 0` (see its documentation) |
| z-score, rolling statistics | Pure Dart — planned |
| ARIMA / SARIMA fit + forecast | Native C core (ctsa) via FFI — available |
| Robust forecast + confidence band + order ladder | `arimaForecast` (pure-Dart wrapper over the native fit) — available |
| Auto-ARIMA order selection (KPSS for d, seasonal strength for D, stepwise or exhaustive search, AIC/AICc/BIC, search journal) | `autoArima` — pure Dart order search over the native fit — available (see [below](#automatic-order-selection-autoarima)) |
| KPSS level-stationarity test, classical decomposition + seasonal strength F_S, root-location check of lag polynomials, information criteria | Pure Dart — available |
| SARIMAX: regression with (seasonal) ARIMA errors, coefficient SEs/covariance | `fitSarimax` — native C core (ctsa) via FFI — available |
| SARIMAX forecast + confidence band + order ladder | `sarimaxForecast` — available |

## Usage

```dart
import 'package:tseries/tseries.dart';

void main() {
  // Tier 1 — pure-Dart linear (OLS) trend forecast. Fits a straight line to the
  // last `window` points of an already-regular series and extrapolates it.
  // `horizon` is in steps, not time. Returns point forecasts plus per-step
  // prediction standard errors that widen with the horizon.
  final lin = linearForecast(series, window: 4, horizon: 3);
  print(lin.slope);     // trend per step
  print(lin.intercept); // line value at the window's local origin
  print(lin.point);     // 3 point forecasts
  print(lin.se);        // 3 prediction standard errors (grow with h)

  // Fixed-order fit (e.g. the classic airline model on log data):
  final fit = fitSarima(
    logSeries,
    order: const ArimaOrder(
      p: 0, d: 1, q: 1, seasonalP: 0, seasonalD: 1, seasonalQ: 1,
      seasonalPeriod: 12,
    ),
    horizon: 12,
  );
  print(fit.ma);   // MA coefficients
  print(fit.aic);  // model AIC

  // Robust forecast with a confidence band and an order-degradation ladder.
  // `series` must already be REGULAR (evenly spaced, no gaps) — resampling and
  // imputation are the caller's responsibility, not the package's. `horizon` is
  // in steps, not time. The ladder is tried in order; the first order that fits
  // is used and reported back.
  final f = arimaForecast(
    series,
    orders: const [
      ArimaOrder(p: 2, d: 1, q: 1), // preferred
      ArimaOrder(p: 1, d: 1, q: 1), // fallback
      ArimaOrder(p: 0, d: 1, q: 1), // last resort
    ],
    horizon: 6,
    confidenceLevel: 0.95, // z ≈ 1.96
  );
  print(f.order);   // which candidate actually fitted
  print(f.point);   // 6 point forecasts
  print(f.ciLower); // lower band: point − z·stdErr
  print(f.ciUpper); // upper band: point + z·stdErr
}
```

### Regression with ARIMA errors (SARIMAX)

```dart
// y_t = Σ β_j·x_{j,t} + u_t,  u_t ~ ARIMA(1,1,1). Regressors are named columns
// of the same length as the series; the map order is the column order.
final fit = fitSarimax(
  series,
  {'temperature': temperature, 'promo': promo},
  order: const ArimaOrder(p: 1, d: 1, q: 1),
  horizon: 6,
  // The forecast is conditional on the regressors' future values.
  futureExog: {'temperature': tempNext6, 'promo': promoNext6},
);
print(fit.coefficient('temperature'));    // β
print(fit.regressors.first.standardError); // its SE (null if not trustworthy)
print(fit.covariance?.parameters);         // [ar.L1, ma.L1, temperature, promo]
print(fit.forecast);                       // 6 point forecasts

// Same with a fallback ladder and a 95% band, like arimaForecast:
final f = sarimaxForecast(
  series,
  {'temperature': temperature},
  orders: const [ArimaOrder(p: 2, d: 1, q: 1), ArimaOrder(p: 0, d: 1, q: 1)],
  horizon: 6,
  futureExog: {'temperature': tempNext6},
  exogPolicy: ExogPolicy.dropDegenerate, // drop unusable columns, don't throw
);
print(f.order);           // the rung that fitted
print(f.fit.regressors);  // which columns were used / dropped and why
```

* **Model.** Regressors are differenced together with the series (d and D) —
  the standard "regression with ARIMA errors" (statsmodels
  `SARIMAX(trend='n')`, R `arima(xreg=)`). An intercept is estimated only when
  `d + D == 0` and `includeMean` is true (the default); `includeDrift` adds a
  `drift` column `1..n` (allowed for `d + D ≤ 1`).
* **Methods.** `SarimaxMethod.cssMle` (default; with no regressors it is
  bit-for-bit `fitSarima` with `ArimaMethod.mle`), `mle` (exact MLE from a zero
  start), `css`.
* **Regressor screening.** Before the native fit every column is checked on
  its *differenced* values: all-zero, constant after differencing (e.g. a
  constant under `d ≥ 1`), a constant duplicating the intercept, or within 1e-8
  of a linear combination of the intercept and earlier columns. With
  `ExogPolicy.strict` (default) the first such column throws an
  `ExogenousColumnError` naming it; with `ExogPolicy.dropDegenerate` it is
  dropped and reported in `fit.regressors` with a null coefficient.
* **Guards.** Non-finite values (also after differencing), fewer than 2
  residual degrees of freedom, and infeasible seasonal orders (see below) are
  refused with an `ArgumentError` before ctsa is called.
* **Regressor units don't matter.** ctsa optimises β in raw units and on its own
  is not scale-invariant (a regressor ×1e6 gave status 1 at a clearly worse
  optimum). The shim hands each column to ctsa multiplied by an exact power of
  two matched to the differenced series' scale and maps β and its covariance
  back; results are invariant to the regressor's unit from 1e-100 to 1e140
  (tested). Values beyond ±1e150 are refused.
* **Forecast SEs** treat the future regressor values as known.
* **Known limitation (shared with `fitSarima`):** the *series'* own scale and
  level are not normalised. Multiplying it by 1e6, or offsetting it by 1e6 over
  unit-size noise, moves the intercept/AR estimates by ~1e-3 relative on both
  paths (ctsa optimises the intercept in raw units). Rescale such series.

### Infeasible seasonal orders

ctsa's exact likelihood needs memory growing with the fourth power of
`max(p + s·P, q + s·Q + 1)`: a seasonal MA term at `s = 288` would need ~7 GB
per likelihood evaluation, and ctsa does not check its allocations. Every
fixed-order entry point (`fitSarima`, `arimaForecast`, `fitSarimax`,
`sarimaxForecast`) refuses an
order whose state exceeds **64 MiB** with a `ModelTooLargeError` (an
`ArgumentError`) instead; the forecast ladders skip such a rung. Seasonal
*differencing* alone (`D = 1`, `P = Q = 0`) stays cheap at any period.

`arimaForecast` returns only numbers. Any domain interpretation (thresholds,
alerts, probability of crossing a level) is the consumer's job, computed from
`point` and `standardErrors` — the package deliberately knows nothing about it.

The public surface takes and returns only plain Dart types. `List<num>` (and
`Map<String, List<num>>` for regressors) in, `ArimaForecast` /
`SarimaFitResult` / `SarimaxFitResult` / `SarimaxForecast` out. No `Pointer`,
no native detail.

### Automatic order selection (`autoArima`)

`autoArima` is pure Dart over the native fixed-order fit `fitSarimax`. It
follows Hyndman & Khandakar (2008, JSS 27(3), §3.1–3.2), with the seasonal
difference chosen as in *Forecasting: Principles and Practice* §9.1:

1. **D** — by the strength of seasonality F_S = max(0, 1 − Var(R)/Var(S + R))
   of a classical additive decomposition (`seasonalStrength`,
   `classicalDecomposition`), stepped by the number of full seasonal cycles
   k = ⌊N/s⌋ (`DifferencingReport.seasonalCycles`): the seasonal difference
   is chosen automatically only with at least 5 full cycles — with 5–9
   cycles only for strong seasonality (F_S ≥ 0.70), from 10 cycles for
   F_S ≥ 0.35. On 2–4 cycles the measure cannot tell seasonality from
   noise: for white noise over 3 cycles it exceeds 0.35 in 73–95 % of
   series. If the seasonality is known from the domain (monthly sales, a
   daily profile), pass `seasonalD: 1`. On short series the D decision can
   differ from `forecast::auto.arima` (R measures on STL with 0.64) — a
   known difference;
2. **d** — successive KPSS tests of level stationarity at 5 %
   (`kpssLevelTest`, Kwiatkowski et al. 1992) on the seasonally differenced
   series. The lag of the long-run variance is chosen by `kpssLag`, re-applied
   to the length T of each step: by default `KpssLag.short` = ⌊3·√T/13⌋, the
   rule that reproduces the differencing decisions of R's
   `forecast::ndiffs`; `KpssLag.l4` / `KpssLag.l12` (KPSS 1992) and
   `KpssLag.fixed(l)` are options. With the default, log AirPassengers gets
   d = 1, D = 1 and the airline model, as in R;
3. a working series that is **constant** after the chosen differences (a
   constant, a straight line) is refused before the search —
   `AutoArimaNoModelException` with `reason: constantSeries` — because every
   model fitted to it would be degenerate or spurious;
4. **(p, q, P, Q, constant)** at fixed (d, D) — the stepwise search of HK08
   §3.2 with the neighbourhood of FPP3 §9.8 (`ArimaSearch.stepwise`,
   default: from the current model, p and/or q ±1 — including opposite
   signs —, P and/or Q ±1, the constant toggled; up to 17 neighbours, the
   seasonal ones first, the first improvement becomes current;
   `maxOrderSum` does not bound it) or an exhaustive grid with
   p + q + P + Q ≤ 5 (`ArimaSearch.exhaustive`), minimising AICc (or AIC/BIC)
   computed from the fit's log-likelihood (`InformationCriterion`);
5. candidates are rejected when too short for the order, infeasible
   (`ModelTooLargeError`), not converged (status ≠ success), degenerate, or
   with an AR/MA root of modulus < 1.01 (`rootMargin`; a seasonal factor in
   z^s against 1.01^s; HK08 used 1.001) (`allRootsOutside`, a step-down
   test that needs no root finding).

```dart
final r = autoArima(logAirPassengers, period: 12, horizon: 12);
print(r.order);            // ARIMA(p,d,q)(P,D,Q)[12]
print(r.constant);         // drift (d + D = 1) or mean (d + D = 0)?
print(r.criterionValue);   // AICc of the selected model
print(r.forecast);         // 12 point forecasts; r.lower / r.upper: 95 % band
print(r.differencing);     // KPSS steps (T, lag rule, lag, η, p bracket), F_S
for (final c in r.trace) { // every candidate with its verdict
  print(c);
}
```

When no model is selected, `AutoArimaNoModelException` carries its `reason`
(`constantSeries`, `noStartModelAccepted`, `noCandidateAccepted`), the whole
journal and the differencing report in its fields; there is no silent
fallback model. The search is synchronous and deterministic (the journal is
byte-identical between runs); a seasonal search can take seconds — run it on
a background isolate (`Isolate.run`) where that matters. There is no time
limit, only a budget of fits (`maxModels`, default 100: an exhaustive
seasonal search with a constant has 192 grid points — raise it).

Cost (measured on 2 200 synthetic series of 200–500 points and on the
reference series; Windows desktop, JIT, indicative only): the stepwise search
takes a median of 13 fits (95th percentile 22; seasonal series: 14, 95th
percentile 35, at most 57) and never reached the default budget. Log
AirPassengers (144 months): 14 fits, ~0.4 s stepwise; 96 fits, ~5 s
exhaustive. Non-seasonal series of ~100 points: 14–19 fits, 20–200 ms.

**Limitations of this version** (deliberately out of scope):

* no regressors (xreg); no Box–Cox; no missing values; no period detection;
* D is chosen by F_S, not by the Canova–Hansen test — expect different D
  than R's `auto.arima` on some series, especially short ones (no automatic
  D below 5 full cycles; constants `autoArimaSeasonalMinCycles`,
  `autoArimaSeasonalManyCycles`, `autoArimaSeasonalStrengthThreshold`,
  `autoArimaSeasonalStrengthThresholdFewCycles`, all provisional). Example:
  on 4-cycle windows of USAccDeaths (51–59 months) R's `auto.arima` takes
  D = 1 and we do not; our forecast of the next year is 1.06–3.9 times worse
  (MASE) than R's there, and with `seasonalD: 1` it is on par with R — pass
  `seasonalD: 1` when the seasonality is known;
* a missed difference is expensive far ahead: on 84–90-point windows of
  WWWusage KPSS (default lag) keeps d = 0, and the stepwise choice forecasts
  10 steps ~30 % worse (median MASE ratio 1.32) than the same search with
  `d: 1` — pass `d` when the series is known to be integrated;
* d by KPSS only (no ADF/PP); p-values are reported as Table 1 brackets,
  never interpolated. At α = 1 % the decision still compares η with the
  Table 1 critical value, whereas R's `ndiffs` never differences at 1 % (its
  p-value is bounded below) — a deliberate difference. The short default lag
  differences persistent stationary series (e.g. AR(1) φ = 0.7) more often
  than l4: a known size/power trade-off of KPSS, chosen because a missing
  difference costs more in forecasting than an extra one;
* the stepwise search can end on a different local minimum than R's
  `auto.arima` (different paths through the neighbourhood); on some series
  it finds a lower AICc — e.g. LakeHuron: ARIMA(2,1,1), the exhaustive
  optimum (AICc 213.51), where R's stepwise stops at ARIMA(0,1,0) (220.26);
* **root check vs. boundary optima.** Candidates are fitted by exact maximum
  likelihood with a CSS-based start. When the true MA part is close to
  non-invertible, the likelihood maximum can sit exactly on the unit circle
  (|MA root| = 1); such fits are correctly rejected by the root check
  (`rootMargin` = 1.01), but the search does not then look for an interior
  local optimum. As a result, on 8 of 2260 synthetic benchmark series a model
  chosen by R `forecast::auto.arima` was rejected here and a different order
  was selected;
* **local optima.** The optimiser may converge to a local rather than global
  likelihood maximum. Example: for an ARIMA(2,1,2) with drift (series
  `ima_n200_r55` in our benchmark) the fit is 3.5 log-likelihood units below
  the reference optimum, which changes the selected model. Starting from zero
  ARMA coefficients recovers the reference optimum in this case; multi-start
  fitting is planned for a future release;
* the seasonal root margin 1.01^s has been checked against R as a black box
  only for s = 12;
* the coverage of the 95 % forecast interval has not been measured; it is
  nominal;
* every candidate is fitted by exact maximum likelihood (no CSS shortcut);
* σ² is the profile MLE S/T, so bands are ~1–3 % narrower than with a
  degrees-of-freedom correction; absolute AIC values differ from R's by up to
  ~3e-3 (R's approximate diffuse initialisation) — compare models within this
  package, not across packages;
* no parallel search, no progress callback.

**Determinism.** Within one build, autoArima results are bit-for-bit
reproducible. The native core is built without floating-point contraction
(`-ffp-contract=off`, no FMA). On the verified platforms (Windows x64,
Android x86_64, Android arm64) we guarantee: the same differencing orders
(d, D) and KPSS statistic; the same selected model and the same search path;
estimates agreeing within 1e-3·max(1,|c|) for coefficients and 1e-4 for the
log-likelihood. We do **not** guarantee bit-identical numbers across
platforms: system math libraries (UCRT, bionic, …) round elementary
functions differently. Occasionally this flips the verdict of a candidate
that was not selected anyway, near the root-admissibility boundary; in our
verification set it never changed the selected model. When the likelihood is
nearly flat (roots near the unit circle, over-differenced series) the
estimate itself is sensitive to rounding-level perturbations: determinism
gives the same answer everywhere, it does not make such an estimate
well-conditioned. iOS and macOS are not verified yet.

The measurements behind this (per platform, binary sizes, the static FMA
check and how to re-run it) are in
[`doc/platform_verification.md`](doc/platform_verification.md).

## Architecture

The package is layered so that native detail never leaks upward:

1. **Raw bindings** (`lib/src/ffi/bindings/`) — bindings to the native C shim.
   Generated from `ffigen.yaml`, private, never exported. The only place
   `dart:ffi` `Pointer` types live.
2. **Safe wrapper** (`lib/src/ffi/wrapper/`) — owns native memory (`Arena` /
   `using`, freed on every path including throws), null-checks and bounds-checks
   native returns, scans every result for NaN/Inf, and turns C status codes into
   typed Dart exceptions. No `Pointer` crosses this boundary.
3. **Domain API** (`lib/src/domain/`) — the clean, idiomatic surface. This is
   the only layer re-exported from `lib/tseries.dart`.

**Consumers only ever import `package:tseries/tseries.dart`.** They cannot tell
whether a given operation is pure Dart or native — and that is the point.

### Native shim (`native/`)

The Dart bindings target a thin C shim (`native/tseries_ctsa.c`) rather than
ctsa directly, for two reasons: (1) it explicitly exports its symbols
(`__declspec(dllexport)` / `visibility("default")`) so `@Native` lookups work on
Windows as well as Android/iOS; and (2) it keeps every ctsa struct on the C side
and returns only flat scalars/buffers, so no fragile struct layout is ever bound
from Dart. See the header comment in `native/tseries_ctsa.h`.

## Native core, vendoring & attribution

The ARIMA/SARIMA capability is **built on
[ctsa](https://github.com/rafat/ctsa)** by Rafat Hussain, a C library for
univariate time-series analysis. The sources are vendored under
`third_party/ctsa/`, pinned to a specific upstream commit; see
`third_party/ctsa/PROVENANCE.md` for the exact revision and rationale.

**License note.** Most of ctsa is BSD-3-Clause, with MIT, MINPACK and CC0
fragments. The units upstream marks **LGPL-3.0-or-later** (`brent.c`,
`erfunc.c`, two sections of `dist.c`), the MIT-tagged SVD routines of `lls.c`
that name no copyright holder, the polynomial root finder (`polyroot.c`, ACM
Algorithm 419 under a non-commercial licence) and `ppsum()` (from the GPL R
package tseries) **have been removed** from the vendored tree and replaced by
tseries' own BSD-3-Clause code in `third_party/ctsa_replacements/` (written
clean-room from published algorithms, checked against scipy/numpy — see
`tool/cleanroom_accuracy.c` and `tool/polyroot_accuracy.c`). Unreachable code
of restrictive origin (Numerical Recipes, AS 197, R's ARIMA transforms) was
deleted, and so was the whole auto-ARIMA cluster (a C version of GPL R code).
`erf`/`erfc` now come from the platform's C math library. Some
third-party code with non-standard or GPL terms remains and is listed under
[License](#license); the per-file audit is in
`third_party/ctsa/PROVENANCE.md` ("License audit"), together with every local
change to the vendored code.

The native library exports only the `tseries_*` shim functions (hidden
visibility elsewhere) and is linked with dead-code elimination, so only the
ctsa code reachable from the shim ships.

## Numerical validation

The native core is validated against published reference fits — "compiles and
returns a number" is not "correct". `test/arima_golden_test.dart` fits the
canonical multiplicative airline model, `SARIMA(0,1,1)(0,1,1)[12]` on
`log(AirPassengers)` (Box & Jenkins Series G), and asserts agreement with the
widely published R `stats::arima` result:

| Statistic | Reference (R) | tseries/ctsa | Tolerance |
| --- | --- | --- | --- |
| MA₁ (θ) | −0.4018 | 0.3991 (sign flipped) | 0.01 on \|θ\| |
| seasonal MA₁ (Θ) | −0.5569 | 0.5544 (sign flipped) | 0.01 on \|Θ\| |
| σ² | 0.001348 | 0.0013479 | 5e-5 |
| log-likelihood | 244.7 | 244.73 | 0.5 |
| AIC | −483.4 | −483.46 | 1.0 |

ctsa parameterises the MA polynomial as `1 − θB` where R uses `1 + θB`, so MA
estimates carry the opposite sign; magnitudes match to well under 1%.

### What that validation does *not* cover

The agreement above speaks for the **ordinary MLE path only**. The airline model
has `p == 0` and `P == 0`, and ctsa only inspects AR roots when `p > 0` (status
`nonStationaryAr`) or `P > 0` (`nonStationarySeasonalAr`). The golden suite is
therefore **structurally incapable** of entering those branches — where ctsa
returns *before* the MLE step and leaves CSS-only estimates behind.

Those branches are reachable in practice, and `test/arima_fit_status_test.dart`
reproduces both and pins their reported behaviour. But they are **not** validated
against R. Read `status` on any result you intend to trust numerically: only
`ArimaFitStatus.probableSuccess` means the estimator you asked for actually ran
to completion. This is also why `aic` is null on those branches — see below.

### Exact likelihood

The log-likelihood is the exact Gaussian likelihood of the differenced series
(stationary start, σ² profiled out), evaluated at the returned estimate, and
checked against statsmodels and an exact Cholesky computation on AR, MA, ARMA,
seasonal and regression-error models. Pure AR(1) error models —
ARIMA(1,d,0) without seasonal AR/MA terms — were wrong in ctsa (the first
observation got variance σ² instead of σ²/(1 − φ²)) and are fixed in the
vendored copy; see `third_party/ctsa/PROVENANCE.md` ("Local modifications"
10) and `test/ar1_likelihood_test.dart`. Like
R, it counts only the n − d − s·D differenced observations (statsmodels'
exact-diffuse `llf` additionally charges −½·ln 2π per differenced-away one).

### White-noise orders and inputs that are refused

Orders without ARMA terms — the random walk (0,1,0), (0,2,0), white noise with
a mean, seasonal random walks `(0,d,0)(0,D,0)[s]` — forecast through the same
AS 182 Kalman forecaster as every other order (as AR(1) with φ = 0); their
forecasts and standard errors are checked against the closed form and
statsmodels in `test/white_noise_test.dart`. In upstream ctsa their standard
errors were uninitialised memory (PROVENANCE.md, "Local modifications" 12).

A few inputs on which ctsa reads or writes outside its buffers are refused
with an `ArgumentError` instead of being run (PROVENANCE.md, "Local
modifications" 12): `fitSarimax` with `cssMle`/`css` when the differenced
series is not longer than p + s·P; Box-Jenkins on a differenced white-noise
order (nothing to estimate) or when p + s·P or q + s·Q exceeds the differenced
length. Two more length limits exist because ctsa would otherwise read
uninitialised memory and return run-to-run different results; `fitSarima`
refuses them with a `SeriesTooShortError` (and `arimaForecast` skips the
order): MLE/CSS fits whose differenced series is not longer than p + s·P, and
seasonal Box-Jenkins fits with fewer than (P + Q + 1)·s + 2 differenced
points.

### Diagnostics you can trust

`aic` and `logLikelihood` are `double?` and are null whenever ctsa did not
genuinely compute them, rather than being filled with a plausible-looking
number:

| Method | `logLikelihood` | `aic` |
| --- | --- | --- |
| `mle`, status `probableSuccess` | real | real |
| `css`, status `probableSuccess` | real (a CSS log-likelihood) | **null** — ctsa computes none |
| `boxJenkins` | **null** | **null** — ctsa computes neither |
| any method, status 10/12 | **null** | **null** — would be a CSS-based, incomparable AIC |

A fit is never *rejected* on the strength of a status code — the code is
reported, and what to do about it is the caller's decision. A fit **is** rejected
when its residual variance is ≤ 0 (a constant/flatlined series), because such a
model reports zero forecast uncertainty.

### SARIMAX: what is validated and what is not

`test/sarimax_golden_test.dart` compares `fitSarimax` with statsmodels
`SARIMAX(simple_differencing=False, trend='n')` on four fixed datasets
(generated by `tool/sarimax_reference.py`, which records the library versions
in the fixture it writes):

| Case | Coefficients | log-lik | Forecast | Forecast SE | Coefficient SE |
| --- | --- | --- | --- | --- | --- |
| ARIMA(1,1,1) + 2 regressors | 1.4e-5 | 5.9e-5 | 8.0e-6 | 2.4e-5 rel | 2.8e-5 rel |
| `SARIMA(0,1,1)(0,1,1)[12]` + 1 regressor | 6.8e-6 | 7.0e-5 | 6.9e-6 | 4.6e-6 rel | 1.7e-5 rel |
| ARMA(1,0,1) + intercept + 1 regressor (d = 0) | 1.4e-5 | 3.6e-5 | 2.0e-5 | 5.4e-6 rel | 1.4e-5 rel |
| ARIMA(0,1,1) + drift | 5.7e-6 | 3.0e-5 | 2.1e-5 | 7.4e-6 rel | 9.8e-6 rel |

(largest absolute differences, measured 2026-10-04; the tests assert 2e-4 / 2e-3 / 2e-4 / 2e-3 rel /
2e-3 rel, just above the ~1e-5 by which the reference itself moves when
re-polished by a second optimiser). MA signs are flipped (ctsa: `1 − θB`).

The coefficient covariance is the inverse of a **numerical** Hessian of the
profiled likelihood. It matches statsmodels `cov_type='approx'` (table
above) and an independently computed profile Hessian. It does **not** match
statsmodels' default-for-`fit` `cov_type='oim'` for the ARMA terms (4–15% off):
that one is Harvey's analytic approximation, a different estimator — not an
error on either side. Regression-coefficient and intercept SEs agree with all
three to < 2%.

**Not validated against a reference:**
* `SarimaxMethod.mle` (pure MLE start) — only checked to land on the same
  optimum as `cssMle` on a well-posed problem; `SarimaxMethod.css` — only
  sanity-checked. Neither has a covariance (CSS) or was compared numerically.
* Seasonal AR terms (P > 0) with regressors, multiple seasonal terms, and
  periods other than 12.
* Every status other than `probableSuccess`: as on the SARIMA path, ctsa
  returns CSS-only estimates on 4/10/12/15 under `cssMle`, and `aic`,
  `logLikelihood` and `covariance` are null there — reported, never trusted.
* The covariance's off-diagonal entries beyond symmetry (only the diagonal is
  compared).

## Building & testing

The native build is driven by a Dart **build hook** (`hook/build.dart`) using
`package:native_toolchain_c`, producing a `ctsa` code asset. No OS-specific
build files are required. A C toolchain is needed (MSVC on Windows; clang on
macOS/Linux; the NDK for Android).

> Build hooks / code assets are **stable** as of Dart 3.10 (Flutter 3.38). This
> package requires `sdk: ^3.10.0`, so the native build and `dart test` run on the
> standard stable SDK — no `--enable-experiment` flag and no beta channel needed.

The hook compiles with floating-point contraction off (`-ffp-contract=off`;
`/fp:precise` on MSVC) — see **Determinism** above. Gate E0 checks a built
binary for fused multiply-add instructions (expected: 0; needs `llvm-objdump`,
taken from the Android NDK if present); the tool lives in the repository, not
in the published archive, and
[`doc/platform_verification.md`](doc/platform_verification.md) shows how to
run it per platform:

```
dart run tool/fma_count.dart path/to/libctsa.so
```

Regenerate the raw bindings (requires libclang / an LLVM install):

```
dart run ffigen --config ffigen.yaml
```

## License

The `tseries` package itself (Dart code, the C shim in `native/` and the C
replacements in `third_party/ctsa_replacements/`) is licensed under
**BSD-3-Clause** — see [`LICENSE`](LICENSE).

The native core is **[ctsa](https://github.com/rafat/ctsa)** by Rafat Hussain
(BSD-3-Clause), vendored under `third_party/ctsa/`. Its notices, together with
the MIT, BSD-style and CC0 fragments it contains, are reproduced in
[`THIRD_PARTY_NOTICES`](THIRD_PARTY_NOTICES). No LGPL code remains in the package.

**Licensing status.** One open item,
recorded in `THIRD_PARTY_NOTICES` and `third_party/ctsa/PROVENANCE.md` ("License audit",
"Open licensing questions"); it is an engineering note, not a legal
conclusion: the exact-likelihood and forecasting core used by **every** fit
contains translations of Applied Statistics algorithms AS 154 / AS 182 /
AS 75, which the Royal Statistical Society permits to be distributed
"provided that no fee is charged" (`neldermead.c`, AS 47, is linked but never
used).

ctsa's GPL-derived auto-ARIMA code, its unaudited `sarimax_wrapper()` (never
called by tseries), `boxcox.c` and the unreachable `initest.c` (Burg, no
licence found for the original) are not in the package.
