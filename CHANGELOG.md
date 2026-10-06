Formerly `tsforge`. Entries up to and including 0.7.0 were written under the
old name and are kept as they were.

## 0.8.0 — 2026-10-06

- `LICENSE` is now the plain BSD-3-Clause text; the third-party notices for the
  vendored ctsa tree moved to `THIRD_PARTY_NOTICES` (pub.dev license detection
  needs an unmodified license text).
- **New: `autoArima` — automatic (seasonal) ARIMA order selection in pure
  Dart** over the native fixed-order fit, the replacement for the removed
  `autoArimaForecast` (not a drop-in: different name, different algorithm,
  orders can differ). Hyndman & Khandakar (2008) stepwise search with the
  neighbourhood of FPP3 §9.8 — up to 17 neighbours, p and/or q ±1 with
  either sign, P and/or Q likewise (seasonal ones first), the constant
  toggled — (or an
  exhaustive grid, `ArimaSearch.exhaustive`) at (d, D) fixed beforehand — d
  by successive KPSS tests (lag rule `kpssLag`, default `KpssLag.short` =
  ⌊3·√T/13⌋ re-applied at every step, which reproduces R's `ndiffs`), D by
  the strength of seasonality of a classical decomposition instead of
  Canova–Hansen, stepped by full cycles k = ⌊N/s⌋ (no automatic D below 5
  cycles, F_S ≥ 0.70 for 5–9, F_S ≥ 0.35 from 10; reported as
  `DifferencingReport.seasonalCycles` / `seasonalStrength` /
  `seasonalThreshold`); a series that is constant after
  differencing (constant, straight line) is refused before the search; AICc/AIC/BIC computed
  from the fit's log-likelihood; candidates rejected for length, size,
  non-convergence, degeneracy or AR/MA roots within 1.01 of the unit circle
  (`rootMargin`).
  Returns `AutoArimaResult` with the fit, forecast and band, a
  `DifferencingReport` and the full search journal (`CandidateRecord` per
  candidate with its `CandidateVerdict`); throws `AutoArimaNoModelException`
  (`reason`: `constantSeries` / `noStartModelAccepted` /
  `noCandidateAccepted`; journal in its fields) when no model is selected. Deterministic within a
  build (the journal is byte-identical run to run); across the verified
  platforms the same (d, D), search path and selected model, with estimates
  equal within stated tolerances, not bit for bit (README, "Determinism").
  No time limit, a fit
  budget (`maxModels`). Written clean-room from the published
  methodology (`doc/auto_arima_clean_room.md`). Known limitations
  (D on short series, boundary MA optima rejected by the root check, local
  likelihood optima, the seasonal margin checked only at s = 12, unmeasured
  interval coverage, cross-platform ties) are listed in the README,
  "Automatic order selection (`autoArima`)" → "Limitations of this version".
- **New pure-Dart primitives** used by `autoArima` and exported on their own:
  `kpssLevelTest` / `KpssResult` / `KpssLag` (`short`, `l4`, `l12`,
  `fixed`) / `kpssLevelCriticalValues` / `constantSeriesTolerance` (KPSS 1992
  level-stationarity test; η agrees with statsmodels to 1e-8 and with R's
  `urca::ur.kpss` at lags 0–8; at α = 1 % the table value is used, unlike
  R's `ndiffs`), `classicalDecomposition` / `seasonalStrength`
  (classical additive decomposition and F_S; agrees with statsmodels'
  `seasonal_decompose`), `allRootsOutside` / `minRootModulus` (whether every
  root of a lag polynomial lies outside a circle — step-down / Schur–Cohn in
  double-double arithmetic, no root finding), and `InformationCriterion`
  (AIC, AICc, BIC from a log-likelihood).
- **Build: the native core is compiled without floating-point contraction**
  (`-ffp-contract=off` for clang/GCC on Android, Linux, macOS and iOS;
  explicit `/fp:precise` on MSVC). clang contracted `a * b + c` into FMA on
  arm64 by default (497 fused instructions in `libctsa.so`), which made one
  of 166 autoArima verification runs on Android arm64 select a different
  model than Windows; without contraction all 166 agree on every verified
  platform. New development gate E0, `tool/fma_count.dart`: counts fused
  multiply-add instructions in a built binary, expected 0 on every ABI
  (`doc/platform_verification.md`).

- **Fix: process crash (access violation) in `fitSarimax` under the Dart
  JIT.** The FFI wrapper no longer reads or writes native buffers through
  `Pointer.asTypedList` views; it accesses them element by element through
  the `Pointer`. Cause — a Dart VM defect (3.10.9, JIT only: `dart run`,
  `dart test`, Flutter debug builds; AOT/release never deoptimizes and was not
  affected): the optimizer scalar-replaces a non-escaping `asTypedList` view,
  and when the optimized code is deoptimized the VM rebuilds the view with its
  data address stored as a tagged integer — the view then points to twice the
  real address. In `fitSarimax` the first fit with a regressor (e.g. the
  drift column) after enough fits without one lazily deoptimized the wrapper
  while its diagnostics view was live, and the next read faulted
  (0xC0000005). Found by the autoArima length sweep (case 1948 after 1947
  others); proven with the VM's deoptimization trace (`null Field @
  offset(8) <- <address>`) and the faulting address (exactly 2 × the
  buffer). Results are bit-identical (6 787-fit dump). Tests:
  `test/no_typed_data_views_test.dart` (runs the crashing sequence in a child
  VM with flags that make it deterministic, and keeps `asTypedList` out of
  `lib/`). Child processes of the tests (this one and the heavy gates E/F)
  now run through `test/support/child_vm.dart` — kernel + bare VM with the
  native-assets mapping `dart test` already wrote — instead of `dart run`,
  which failed to re-bundle `ctsa.dll` while the test runner had it loaded.
  Gate F's timeout is 60 min (the sweep alone takes ~19 min).
- **Fix: heap overflow in ctsa's SARIMAX regression set-up refused.** With
  seasonal differencing (D > 0), regular differencing d and r regressors in
  use, ctsa wrote past a heap buffer whenever d·r exceeded the differenced
  length n − d − s·D (e.g. 3 regressors under (0,2,0)(0,1,0)[4] on 11
  points). `fitSarimax` now throws `ArgumentError` for those inputs
  (PROVENANCE.md, "Local modifications" 16; test in
  `test/short_series_guard_test.dart`). Nothing else changes.

- **Fix: short series on which ctsa was non-deterministic are refused, with a
  new typed `SeriesTooShortError`** (an `ArgumentError`, exported). ctsa read
  uninitialised memory — so status and estimates could differ from one call
  to the next, also for fits reported as successful — in two cases, now
  refused by `fitSarima` before the native call (and by the shim):
  - MLE/CSS when the differenced series is not longer than p + s·P (e.g.
    (1,0,0)(1,0,0)[12] on ≤ 13 points);
  - seasonal Box-Jenkins when the differenced series has fewer than
    (P + Q + 1)·s + 2 points (ctsa's start values need autocovariances up to
    lag (P + Q)·s; one season of margin on top of the measured bound).
  `arimaForecast` skips such an order with a note in `attempts` ("below
  structural minimum …") and moves on; it throws `ArgumentError` only when no
  order is long enough. Orders and lengths outside these cases are unchanged
  bit for bit; Sweet Limit's ladder is far from them. (PROVENANCE.md, "Local
  modifications" 15; test `test/short_series_guard_test.dart`.)
- **Fix: log-likelihood and AIC use the exact π.** ctsa evaluated
  ln(2·3.14159), so every `logLikelihood` was too high (and every `aic` too
  low) by ½·T·ln(π/3.14159) ≈ 4.2e-7·T, T = n − d − s·D (6e-5 at T = 144).
  Now exact (PROVENANCE.md, "Local modifications" 14): `logLikelihood`
  decreases by exactly that amount, `aic` increases by twice it; nothing else
  changes (coefficients, σ², forecasts, standard errors bit-identical).
  Comparisons between models on the same series are unaffected. Absolute
  values now agree with statsmodels' `loglike` at the same parameters to
  ~1e-11 relative without any correction.
- **Fix: σ², log-likelihood and AIC are now those of the returned estimate
  (Windows/MSVC builds).** ctsa copied them from its last objective
  evaluation, which happens inside the finite-difference Hessian — on MSVC
  builds a point one step (≈ 6e-6 relative) away from the estimate in the
  last parameter. Now the objective is re-evaluated at the estimate
  (PROVENANCE.md, "Local modifications" 13). Affected: every `fitSarima` and
  `fitSarimax` fit (all methods) and the `sigma2`/`logLikelihood`/`aic` of
  `arimaForecast`, **on Windows desktop builds only** — clang/GCC builds
  (Android, iOS, macOS, Linux) already evaluated the estimate last and are
  unchanged bit for bit. On ordinary fits log-likelihood/σ²/AIC move by
  ~1e-5 (reference fits ≤ 5.6e-5 in ll, ≤ 2.3e-5 relative in σ²; Sweet
  Limit's ladder ≤ 4.2e-5 / 1.2e-4); on degenerate fits (a level series
  without a mean, a near-unit MA root) by much more, because the objective is
  steep there. **Coefficients, their covariance and standard errors, point
  forecasts and forecast standard errors do not change** (bit-identical;
  the forecaster estimates its own σ² from its one-step residuals). After the
  fix σ² and the log-likelihood equal the same objective evaluated
  independently at the estimate bit for bit, agree with an exact
  Cholesky-based Gaussian likelihood to ≤ 2e-13 / 6e-13 relative, and are
  bit-identical across MSVC, clang and GCC for identical estimates.
- **Fix: forecasts of white-noise orders (p = q = P = Q = 0) — the random
  walk (0,1,0), (0,2,0), white noise with a mean, seasonal random walks
  (0,d,0)(0,D,0)[s], with or without drift/regressors.** ctsa's forecaster
  (AS 182 `forkal`) refuses these orders and returned without writing
  anything, so through `fitSarima`, `fitSarimax`, `arimaForecast` and
  `sarimaxForecast`:
  - the **forecast standard errors were uninitialised memory** — a different
    value on every call, sometimes NaN (then `TseriesNumericException`), and
    so were the confidence bands;
  - for a differenced order (d + D > 0) the **point forecast was 0** (the
    contents of the zeroed output buffer) instead of the last observation
    (plus drift). With d = D = 0 the forecast (the mean) was right.
  Fixed in the shim (no vendored change; PROVENANCE.md, "Local
  modifications" 12): the shim now forecasts through `forkal` itself and
  passes white noise as the AR(1) model with φ = 0, which is the same model.
  Forecasts and SEs now equal the closed form (σ²·Σψ²) to ≤ 4e-15 relative
  and statsmodels to the same (5e-10 where a mean/drift is estimated); σ² is the one-step residual variance, as for every other
  order. Not affected: every other order — bit-identical forecasts — and
  Sweet Limit's (2,1,1)/(1,1,1)/(0,1,1) ladder. New tests:
  `test/white_noise_test.dart`, fixture `test/fixtures/white_noise_reference.dart`
  from `tool/white_noise_reference.py`.
- **Fix: process crashes from out-of-bounds memory access in ctsa.** Found
  with AddressSanitizer/MemorySanitizer on a full order × length × method
  matrix; all fixed or refused in the shim (PROVENANCE.md, entry 12):
  - every forecast of a model with d + s·D larger than (number of
    coefficients)² — all differenced white-noise orders, and e.g.
    (0,0,1)(0,1,0)[12] — read up to d + s·D doubles past a heap block. The
    values were never used, but the read faults when the block ends at the
    edge of committed memory: this is the intermittent crash seen in long
    white-noise runs. No longer happens; results unchanged.
  - `fitSarimax` with `cssMle`/`css` on a series whose differenced length is
    not larger than p + s·P (e.g. (0,0,0)(1,1,0)[12] on 15 points) wrote past
    a heap block — heap corruption that could crash the process at any later
    point. Now an `ArgumentError` (series too short); such fits never
    produced a usable model (the CSS sum is empty).
  - `fitSarima` with `boxJenkins` crashed the process (heap corruption, not
    `exit()` as previously documented) for a differenced white-noise order
    (nothing to estimate) — now an `ArgumentError` naming the problem — and
    read out of bounds when p + s·P or q + s·Q exceeds the differenced length
    — now an `ArgumentError` (series too short). White noise with a mean
    still fits by Box-Jenkins.
  - **Who was affected** (Dart stress runs against the unpatched library;
    the VM died with an access violation, "GetAndValidateThreadStackBounds
    failed"): Box-Jenkins on a differenced white-noise order crashed every
    run within its first such fit; seasonal SARIMAX CSS fits on 15–20 points
    crashed every run within ~1 000 fits; forecasting white-noise orders on
    ordinary lengths (n ≥ 30, no Box-Jenkins) crashed 1 run in 6 (~190 000
    fits) — the faulting instruction is the over-read in
    `sarimax_predict()`. All other orders on n ≥ 30 without Box-Jenkins:
    no crash in 40 000 fits. Sweet Limit (non-seasonal (2,1,1)/(1,1,1)/
    (0,1,1), MLE, long series, no white noise) was not affected.
  - SARIMA fits whose differenced series is not longer than p + s·P, and
    seasonal Box-Jenkins fits of short series, read uninitialised values and
    gave run-to-run-different results — now refused (next entry).
- **BREAKING: `autoArimaForecast` and `ForecastResult` removed.** The native
  automatic order search they ran was a C version of GPL-licensed R code (the
  R packages forecast, tseries and urca, R's STL, Friedman's super smoother)
  and could not ship in a BSD-3-Clause package. Removed together with it: the
  shim entry point `tseries_auto_arima_forecast` (the library now exports four
  `tseries_*` functions), the vendored `autoutils.c`, `unitroot.c`,
  `seastest.c`, `stl.c` (+ headers), the auto-ARIMA functions of `ctsa.c`/
  `ctsa.h` and `supsmu` of `talg.c`/`talg.h` (see PROVENANCE.md "Local
  modifications" entry 9). Its replacement is `autoArima` (order search in
  pure Dart over the package's own fit; first entry above). Sweet Limit
  never called `autoArimaForecast`.
  - The fixed-order path is unchanged bit for bit: 635 fits (all methods,
    mean/drift on and off, 6 series, 9 orders) return byte-identical results,
    and every surviving native function has the same size on arm64, armv7 and
    x86_64.
  - Library size (stripped, NDK 29): arm64 245 624 → 185 800 B, armv7
    194 940 → 138 016 B, x86_64 320 000 → 243 920 B; Windows DLL
    552 448 → 470 528 B.
- **Fix: exact likelihood of pure AR(1) error models — numbers change for
  those models.** ctsa's likelihood gave the first (differenced) observation
  of an AR(1) error process variance σ² instead of σ²/(1 − φ²), so the
  maximised objective was not the exact likelihood, and φ, σ², the
  log-likelihood/AIC, and through them the forecasts and their standard
  errors were off. Fixed in the vendored `emle.c` (one line in each of three
  places; PROVENANCE.md, "Local modifications" 10).
  - **Affected:** MLE and CSS-MLE fits (`fitSarima` with `mle`, `fitSarimax`
    with `cssMle`/`mle`, `arimaForecast`, `sarimaxForecast`) of
    **ARIMA(1,d,0) without seasonal AR or MA terms** — any d, any seasonal
    differencing D (e.g. (1,1,0)(0,1,0)[12]), with or without intercept, drift
    or regressors. Not affected: every other order (including Sweet Limit's
    (2,1,1)/(1,1,1)/(0,1,1)), CSS-only and Box-Jenkins fits, and fits whose
    status is 4/10/12 (ctsa stopped before the MLE step).
  - **Size of the error before the fix** (34 reference cases against
    statsmodels): log-likelihood off by up to 2.6 with an intercept or after
    differencing (typically 0.01–0.5), and by much more for an undifferenced
    AR(1) without an intercept with φ near 1 (Nile: 16.4; a synthetic
    φ = 0.995 series: 134); σ² typically 0.5–4% too large in the first group
    (11% for (1,0,0)(0,1,0)[12] on USAccDeaths) and up to ×3.9 in the second;
    φ off by up to 0.07. AIC comparisons between AR(1)
    and other orders were biased by the same amounts.
  - **After the fix:** the likelihood function agrees with statsmodels to
    6e-12 relative; fitted optima within 1.6e-5 in the coefficients and
    1 ± 2.3e-5 in σ²; forecasts and SEs at the fitted parameters equal
    statsmodels' to ~1e-15 relative.
    Every other fit is unchanged bit for bit (6 787-fit dump). New tests:
    `test/ar1_likelihood_test.dart`, fixture `test/fixtures/ar1_reference.dart`
    from `tool/ar1_likelihood_reference.py`.
  - Unchanged: every log-likelihood still carries a +4.2e-7·n offset (ctsa's
    π ≈ 3.14159).
- **Removed from the vendored tree: ctsa's `sarimax_wrapper*` group and
  `boxcox.c`/`boxcox.h`.** Never called by tseries (the shim uses
  `sarimax_init`/`exec`/`predict`/`free`), never audited, and already dropped
  by the linker; `boxcox.c` lost its last caller with it. No behaviour change:
  the fit dump is identical, library sizes and exports unchanged
  (PROVENANCE.md, "Local modifications" 11).
- **Build-hook dependencies raised:** `hooks` ^2.0.0 (was ^1.0.1),
  `code_assets` >=1.2.1 <3.0.0 (was ^1.0.0), `native_toolchain_c` ^0.19.2
  (was ^0.17.4). The lower bounds are the newest that Flutter 3.38 can
  resolve (its SDK pins `meta` 1.17.0); a pure-Dart package resolves the
  latest (`hooks` 2.2.0, `code_assets` 2.1.0, `native_toolchain_c` 0.19.5).
  The test suite passes on both ends. No change to `hook/build.dart` was
  needed.
- **Removed from the vendored tree: `initest.c`/`initest.h`** (Burg,
  Yule-Walker and Hannan-Rissanen estimators; the Burg routine is a
  translation of C++ code whose licence was never found) and their callers in
  `ctsa.c` (`ar`, `ar_exec`, `ar_estimate`, `model_estimate`, `pacf`,
  `pacf_opt`, `yw`, `burg`, `hr`). Never reachable from tseries; the stripped
  Android libraries are byte-identical before and after (PROVENANCE.md,
  "Local modifications" 17).
- **Renamed the package `tsforge` → `tseries`** (decided 2026-10-05, ahead of a
  pub.dev release under that name). No behaviour change; the version stays
  0.7.0 until a release is cut. Breaking for consumers:
  - import `package:tseries/tseries.dart` instead of
    `package:tsforge/tsforge.dart`;
  - exceptions renamed: `TsforgeException` → `TseriesException`,
    `TsforgeNumericException` → `TseriesNumericException`,
    `TsforgeInternalException` → `TseriesInternalException`,
    `TsforgeNativeUnavailableException` → `TseriesNativeUnavailableException`;
  - the native code asset id is now `package:tseries/ctsa` (library file name
    `ctsa` unchanged).
  - Internals renamed to match: shim `native/tseries_ctsa.{c,h}`, C symbols
    `tseries_*`, macros `TSERIES_*`; bindings regenerated with ffigen.
- **Added `LICENSE`.** The package (Dart code and C shim) is BSD-3-Clause; the
  file also carries the third-party notices for the vendored ctsa core
  (BSD-3-Clause, Rafat Hussain) and its MIT/CC0 fragments. (At the time the
  tree still contained LGPL units; see the next entry.)
- **LGPL code removed from the vendored ctsa tree.** `brent.c`, `brent.h`,
  `erfunc.c` and the two LGPL-3.0-or-later sections of `dist.c` (`r8_max`,
  `betainv`) are deleted; tseries' own BSD-3-Clause replacements with the same
  signatures live in `third_party/ctsa_replacements/` (`brent_local_min`,
  `erfinv`/`erfcinv`, `betainv`, `r8_max`), written clean-room from published
  algorithms (Brent 1973; Acklam + Halley refinement; Abramowitz & Stegun
  26.5.22 + safeguarded Halley). `erf`/`erfc` now come from libm instead of
  ctsa's same-named definitions; the unused `brent_zero` was dropped. Accuracy
  against scipy: `erfinv`/`erfcinv` <= 5e-16, `betainv` <= 2e-12, `tinv`
  <= 5e-13 relative (`tool/cleanroom_reference.py`,
  `tool/cleanroom_accuracy.c`, fixture `test/fixtures/cleanroom_reference.txt`).
  All golden tests pass at unchanged tolerances.
- **Unattributed MIT routines removed** from `lls.c` (`svd_gr`, `svd_gr2`,
  `minfit`; no copyright holder exists to credit). `minfit` — linked but never
  executed — is replaced by a BSD-3-Clause wrapper over ctsa's own `svd()`.
- `LICENSE` updated: the LGPL status block and LGPL text are gone; added the
  missing notice for `svd_sort` (Pavel Sakov / CSIRO); listed the remaining
  open licensing questions (AS 154/182 "no fee" condition, ACM Algorithm 419,
  code derived from GPL R sources).
- **ACM Algorithm 419 replaced.** `polyroot.c` (CPOLY, used upstream "under
  ACM Software License Agreement for Non-commercial Use") is deleted;
  `polyroot()` — used by the stationarity/invertibility checks of every fit
  and by auto-ARIMA's model screen — is now tseries' own BSD-3-Clause
  Ehrlich–Aberth root finder (`third_party/ctsa_replacements/polyroot.c`,
  Bini 1996 starting points and stopping rule). Checked on 207 polynomials up
  to degree 289 against numpy roots polished in 50-digit arithmetic
  (`tool/polyroot_reference.py`, `tool/polyroot_accuracy.c`, fixture
  `test/fixtures/polyroot_reference.txt`): backward error <= 4.6 (n+1)·eps,
  forward error <= 5.2·cond·n·eps; `archeck()` verdicts identical. The upstream
  code, run through the same harness, was up to 1e11 times less accurate on
  seasonal polynomials and **killed the process** for degree >= ~120 (seasonal
  periods of that size in auto-ARIMA); it also left its outputs uninitialised
  for a zero leading coefficient or NaN input, which `archeck()` then read —
  both fixed by the replacement. `cpoly`/`cpolyroot` had no external callers
  and were dropped.
- **`ppsum()` replaced** (from the GPL R package tseries, © Trapletti): own
  Newey–West/KPSS implementation, bit-identical to the upstream binary on
  20 000 random cases.
- **`archeck()`/`invertroot()` replaced** (the stationarity check and MA
  root reflection run by every fit; upstream's were C versions of R's GPL-2+
  `arCheck()`/`maInvert()`). New BSD-3-Clause implementation
  `third_party/ctsa_replacements/arma_roots.c`, written clean-room from a
  textbook specification by an author who had not read either version; the
  upstream bodies are cut from `talg.c`. Differential test against the
  upstream binary on 20 000 polynomials (`tool/arma_roots_accuracy.c`,
  fixture `test/fixtures/arma_roots_reference.txt`): identical `archeck()`
  verdicts and `invertroot()` statuses on every finite input outside the
  rounding band of |z| = 1; coefficients within 1e-12 or, on 86
  ill-conditioned cases, within the problem's own sensitivity to a 4-ulp
  input change. Two deliberate differences: a non-finite coefficient (or a
  root-finder failure) now makes `archeck()` report non-stationary (fit
  status 10/12) instead of stationary; and roots are reflected as 1/conj(r),
  which fixes upstream's broken output when rounding put only one root of a
  complex pair inside the unit circle (e.g. a seasonal MA coefficient of
  exactly ±1).
- **Unreachable code of restrictive origin deleted from the vendored source**
  (the linker already dropped it, but the source ships): Numerical Recipes
  `tred2`/`tqli`/`pythag` and their caller `eigensystem` (`matrix.c`), AS 197
  `flikam`/`fas197`/`as197` (`emle.c`) with its helper `twacf` (`talg.c`),
  and the C versions of R's `partrans`/`ARIMA_Gradtrans` family (`talg.c`).
- **Licence audit** of every compiled unit, recorded in `PROVENANCE.md`
  ("License audit"). New open findings: auto-ARIMA's order search,
  `autoutils.c` and `seastest.c` follow the GPL-3 R package forecast;
  `archeck`/`invertroot` mirror R's `arCheck`/`maInvert` (since replaced, see
  above); AS 75
  (`inclu2`/`regres`) sits beside AS 154. `LICENSE` now also carries the
  MINPACK notice (`lnsrchmp.c`, `nls.c`) and the second copyright holder of
  the vendored wavelet files.
- Stripped Android sizes now 246/195/320 KB (arm64-v8a / armeabi-v7a /
  x86_64), from 252/201/328 KB.
- **Smaller native library, minimal ABI.** The build hook now compiles with
  `-fvisibility=hidden -ffunction-sections -fdata-sections` and links with
  `--gc-sections` (`-dead_strip` on Apple). The library exports only the five
  `tseries_*` functions instead of ~450 ctsa symbols (which included generic
  names such as `mean`, `var`, `gamma`, `log1p`), and stripped Android sizes
  drop from 526/401/691 KB to 252/201/328 KB (arm64-v8a / armeabi-v7a /
  x86_64). Windows is unchanged (exports were already explicit).
- The C shim includes `<float.h>` before the ctsa headers, removing two MSVC
  C4005 macro-redefinition warnings (`DBL_MAX_EXP`/`DBL_MIN_EXP`) from its
  translation unit.

## 0.7.0

- **New: SARIMAX — regression with (seasonal) ARIMA errors.** `fitSarimax`
  fits `y_t = [μ] + Σ β_j·x_{j,t} + u_t` with `u_t ~ SARIMA(p,d,q)(P,D,Q)[s]`,
  regressors differenced together with the series (statsmodels
  `SARIMAX(trend='n')` / R `arima(xreg=)`). Named regressors
  (`Map<String, List<num>>`), future regressors for the forecast, full seasonal
  orders, optional intercept (`includeMean`, when `d + D == 0`) and drift
  (`includeDrift`, a `drift` column `1..n`, for `d + D ≤ 1`).
  - Result (`SarimaxFitResult`): AR/MA/seasonal coefficients, intercept,
    per-regressor coefficient + SE + defect (`regressors` — the map of which
    columns were used), σ², log-likelihood, AIC, forecast + SEs, ctsa `retval`
    / `ArimaFitStatus`, and the full `CoefficientCovariance` (named parameters).
  - Methods: `SarimaxMethod.cssMle` (default), `mle`, `css`. Mapped explicitly
    in the shim — ctsa numbers its SARIMAX methods differently from SARIMA.
  - Trust gates as on the SARIMA path: `aic`, `logLikelihood` and
    `covariance` are null unless status is `probableSuccess` and the estimator
    actually computed them (no AIC/covariance for CSS; ctsa's CSS Hessian is
    mis-scaled).
  - `sarimaxForecast`: confidence band + order-degradation ladder with
    `attempts`, like `arimaForecast`. The numeric floor counts regressors.
    Strict regressor errors are rethrown, not laddered around.
- **Regressor screening before ctsa sees the data.** On the *differenced*
  columns: all-zero, constant after differencing, constant duplicating the
  intercept, or within 1e-8 of the span of the intercept and earlier columns
  (two-pass Gram-Schmidt). `ExogPolicy.strict` (default) throws
  `ExogenousColumnError` naming the first column and listing every verdict;
  `ExogPolicy.dropDegenerate` drops them and reports them on the result.
- **Regressors are scaled before ctsa sees them.** ctsa is not
  scale-invariant in β (measured: regressor ×1e6 → status 1 at an optimum 52
  log-lik units worse; ×1e-6 → iteration cap). Each column is multiplied by an
  exact power of two matched to the differenced series and β/covariance are
  mapped back; fits are now invariant to regressor units from 1e-100 to 1e140.
  Values beyond ±1e150 are refused (ctsa squares them unscaled).
- **Guards against ctsa's `exit()` and NaN paths:** non-finite input (also after
  differencing, which can overflow finite input), fewer than 2 residual degrees
  of freedom (ctsa's `tinv` would `exit(1)`), unknown method (`exit(-1)`).
- **New: `ModelTooLargeError` on every fixed-order path, SARIMA included.**
  ctsa's exact-likelihood state grows with the 4th power of
  `max(p + s·P, q + s·Q + 1)` and is never checked (a seasonal MA at s = 288
  wants ~7 GB). Orders above 64 MiB are refused before ctsa runs; ladders skip
  such rungs; auto-ARIMA guards its largest permitted order. No change for any
  order under the ceiling. (`ModelTooLargeError` is an `ArgumentError`.)
- **Patched a ctsa defect (PROVENANCE #3):** the SARIMAX collinearity check
  passed a column-major matrix to a row-major `rank()`. It let zero/duplicate
  columns through (NaN or garbage fits) and rejected legitimate regressors held
  for two steps as collinear. Two arguments swapped in `emle.c`; the first half
  is additionally caught by tsforge's own screening.
- **Validated against statsmodels** (`test/sarimax_golden_test.dart`, generator
  `tool/sarimax_reference.py`): non-seasonal + 2 regressors, airline-type
  seasonal + regressor, d = 0 with intercept, drift. Coefficients ≤ 1.4e-5,
  log-likelihood ≤ 7.0e-5, forecasts ≤ 2.1e-5, forecast/coefficient SEs
  ≤ 3e-5 relative. The covariance matches statsmodels `cov_type='approx'`
  (numerical Hessian), not `'oim'` (Harvey's analytic approximation, 4–15% off
  on ARMA terms) — see README. With no regressors `fitSarimax` is
  **bit-identical** to `fitSarima`.
- **Memory: no leaks on any SARIMAX exit.** Host harness (clang, process
  private bytes) through the shim, 20 000 calls per exit — success + forecast +
  covariance, drop policy, strict defect, ctsa's early return before MLE (4),
  CSS, pure MLE, degrees-of-freedom refusal — then a linearity run: private
  bytes rise once by the allocator's working set (≤ 28 KB) and then stay flat
  over a further 80 000 fits (40 000 for the seasonal path). ctsa's own
  retval-7 early return (unscreened, called directly) also measured flat.
  `test/sarimax_leak_test.dart` repeats this over 9 000 calls via RSS and was
  shown to fail when a single `free` is removed from the shim. The one known
  leak, in `rank_c()` on SVD non-convergence, needs a design whose entries span
  ≳ 200 orders of magnitude (measured; see PROVENANCE #3) and is documented,
  not patched.
- `ArimaFitStatus.collinearExogenous` doc updated: it can now surface only as
  the `retval` of a `TsforgeNumericException` from the SARIMAX path.
- Cost (host, clang -O2, n = 144): `cssMle` with 2 regressors ≈ 5 ms per
  fit + forecast, `mle` ≈ 2.7 ms, `css` ≈ 2.2 ms; airline-type seasonal
  (s = 12) with a regressor ≈ 5.8 ms. All guards and screening together cost
  < 0.01 ms.
- Internal: order validation / length floor / z-quantile moved to
  `lib/src/domain/internal/arima_support.dart`, shared by both ladders
  (no behaviour change).
- License: SARIMAX runs ctsa's `regress()` → `tinv()` → `betainv()` (LGPL part
  of `dist.c`) on every fit — no new LGPL unit, but LGPL code is now on the
  fixed-order path. Recorded in PROVENANCE.

## 0.6.0

- **Exposed ctsa's own fit status (`retval`) on the public API.** `SarimaFitResult`,
  `ArimaForecast` and `ForecastResult` now carry `retval` (the raw ctsa code) and
  `status` (an `ArimaFitStatus` enum). `ArimaForecast.attempts` records **every
  rung of the degradation ladder**, in order — each with its order, its ctsa
  code, and whether it was accepted — so callers can build a per-order histogram
  of outcomes and the mix of orders actually used.
  - **Reported, never enforced.** Nothing is rejected on the strength of a
    `retval`, and the control flow is unchanged: a fit with code 4/10/12 is
    still returned exactly as before. This is a measurement step; deciding what
    to reject comes after seeing real numbers.
  - An unrecognised code maps to `ArimaFitStatus.unknown` instead of throwing.
  - `ArimaOrder` gained value equality so it can key an aggregation directly.
- **`aic` and `logLikelihood` are now `double?` and are null unless ctsa really
  computed them.** Previously they were read unconditionally out of a ctsa struct
  that `sarima_init` never initialises, so for methods that do not assign them
  the values were **malloc garbage** presented as real numbers.
  - `aic` is non-null only for `ArimaMethod.mle` with `status ==
    probableSuccess`. ctsa computes no AIC at all for CSS or Box-Jenkins; and on
    codes 10/12 it returns *before* the MLE step, so the AIC it does compute is
    derived from a **CSS** log-likelihood and is not comparable with an MLE
    fit's — which is the entire purpose of an AIC.
  - `logLikelihood` is non-null only for MLE/CSS with `status ==
    probableSuccess`.
  - The nullability is the point: selecting a model by AIC is the obvious next
    step, and a silently incomparable number would corrupt that comparison
    without ever looking wrong.
- **Fixed: the shim no longer rejects a fit based on a field it does not trust.**
  The finiteness gate covered `loglik`/`aic`, which are not the product of the
  fit and are not assigned on every path. It now covers only `mean`/`sigma2`,
  which every estimation method assigns.
- **Degenerate fits are now refused explicitly, for the right reason.** A model
  with residual variance ≤ 0 (a constant/flatlined series) fits its input
  exactly and reports **zero forecast uncertainty** — a confidence band of zero
  width. It is rejected on `sigma2`, a field we trust and actually return, via
  the new `TSFORGE_ERR_DEGENERATE` status → `TsforgeNumericException`. This
  preserves the previous behaviour, which had been rejecting such fits only as
  an accident of `loglik`/`aic` blowing up to ±Inf.
- **Fixed two memory leaks in the vendored ctsa (`emle.c`).** On ctsa's 10/12
  early return, `as154()` leaked `x` and `as154_seas()` leaked `inp2` — the
  latter being the path `sarima_exec` actually calls, and one the app re-enters
  on every re-fit. Measured on the reproducing series: ~640 bytes/fit before,
  ~1.3 bytes/fit (allocator noise) after, over 40 000 fits. See
  `third_party/ctsa/PROVENANCE.md` → "Local modifications".
- **Scope of the "validated against R" claim, stated honestly.** The golden suite
  fits the airline model (0,1,1)(0,1,1)[12], which has `p == 0` and `P == 0`.
  ctsa only inspects AR roots when `p > 0` / `P > 0`, so the golden tests are
  **structurally incapable** of entering the 10/12 branch. The published <1%
  agreement with R therefore speaks **only for the ordinary MLE path**, not for
  the CSS-only fallback that codes 10/12 produce. Those branches now have their
  own coverage (`test/arima_fit_status_test.dart`) proving they are reachable
  and correctly reported, but they are **not** validated against R.

## 0.5.0

- **Fixed: the native core never loaded on Android.** `libctsa.so` was built,
  exported its symbols and was bundled for every ABI — but it was not linked
  against the C math library, so `dlopen` failed at the first native call with
  `cannot locate symbol "log"`. Every native-backed API (`fitSarima`,
  `arimaForecast`, `autoArimaForecast`) therefore failed on-device while the
  pure-Dart ones kept working. `hook/build.dart` now links `m` explicitly on
  every target OS where libm is a separate library (i.e. all but Windows, whose
  CRT provides math and which has no `m.lib`).
  - The Windows host build could not have caught this: math lives in the CRT
    there, so the host suite passed throughout. Verified on an Android device
    (x86_64 emulator: smoke check + the full ARIMA ladder now fit and forecast)
    and, for arm64/armeabi-v7a, at the ELF level (`libm.so` present in
    `DT_NEEDED`, no undefined math symbols).
- Added `checkNativeCore()`, a public health check for the native core. It
  returns `NativeCoreHealth` when ctsa is reachable and functional on this
  device, and throws the new `TsforgeNativeUnavailableException` when the code
  asset is missing or its symbols do not resolve. This makes the failure above
  diagnosable in one call instead of silently degrading, and separates "the
  native build is broken here" from `TsforgeNumericException`'s "this input has
  no good model". It touches no user data and costs microseconds.

## 0.4.0

- Added `linearForecast(...)`, the Tier-1 forecast primitive: a pure-Dart
  ordinary-least-squares (OLS) trend fit over the last `window` points of an
  already-regular series, extrapolated over a `horizon`. Closed-form and
  deterministic — no iteration, no native/FFI dependency (simple arithmetic
  stays in Dart by team agreement).
- It returns a `LinearForecast` carrying the point forecast (`point`), per-step
  prediction standard errors (`se = s·√(1 + 1/N + (x−x̄)²/Sxx)`, which grow with
  the horizon), and the fitted line (`slope`, `intercept`). All lists have length
  equal to `horizon`.
- Input guards raise `ArgumentError`: `window` < 3 (the residual std error needs
  `N−2` degrees of freedom), `horizon` < 1, `series` shorter than `window`, and
  any non-finite (NaN/±Inf) value in the fitting window. A degenerate design
  (`Sxx == 0`) is guarded defensively though it cannot arise for `N ≥ 2` regular
  indices.
- Added `test/linear_forecast_test.dart`: exact-line recovery (`y = 2x+1` →
  slope 2, intercept 1, zero residuals → zero SE), a hand-computed SE reference
  (`y = [1,2,2,5]` → se[0] = 1.5, se[1] = √3.33), SE growth, output shapes, and
  the input guards.

## 0.3.0

- Added `arimaForecast(...)`, a forecast-oriented entry point on top of the
  existing (validated) `fitSarima`. It returns an `ArimaForecast` carrying the
  point forecast, per-step standard errors, a confidence band
  (`ciLower`/`ciUpper` = `point ± z·stdErr` for a configurable
  `confidenceLevel`, default 0.95), the diagnostics (`sigma2`, `aic`,
  `logLikelihood`, coefficients), and the model order that actually produced the
  forecast.
- Added an **order-degradation ladder**: `arimaForecast` takes an ordered list
  of candidate `orders` and returns the first that meets the numeric length
  floor and fits to a finite model, reporting which one worked. The package
  stays generic — the caller supplies the ladder (e.g. `(2,1,1)→(1,1,1)→
  (0,1,1)`); the package hardcodes no specific order.
- Added a per-order **numeric length floor** (≈10 observations per estimated
  parameter, above the native estimability floor): a series shorter than every
  candidate's floor raises a clear `ArgumentError`; a NaN/Inf input or a total
  convergence failure across the ladder raises `TsforgeNumericException`.
- The confidence critical value `z` is computed from `confidenceLevel` via
  Acklam's inverse-normal approximation (so any level in `(0,1)` is supported;
  0.95 → z ≈ 1.96). No native code or FFI symbols were added — the whole feature
  is a pure-Dart wrapper over `fitSarima`.
- Added `test/arima_forecast_test.dart` covering the wrapper contract: numeric
  passthrough (matches a direct `fitSarima`), output shapes, the confidence
  band, ladder degradation on a short series, and input guards.

## 0.2.0

- Added the native ARIMA/SARIMA core: the `ctsa` C library is vendored under
  `third_party/ctsa/` (pinned commit; see `PROVENANCE.md`) and compiled via a
  `hook/build.dart` build hook using `package:native_toolchain_c` into a `ctsa`
  code asset.
- Added a thin, explicitly-exported C shim (`native/tsforge_ctsa.c`) that keeps
  all ctsa structs on the C side and exposes a flat, struct-free ABI.
- Added the three FFI layers: raw bindings (`lib/src/ffi/bindings/`), a
  memory-owning safe wrapper (`lib/src/ffi/wrapper/`) with NaN/Inf and status
  guards, and the idiomatic domain API (`autoArimaForecast`, `fitSarima`,
  `ArimaOrder`, `ForecastResult`, `SarimaFitResult`).
- Added golden-value tests validating the native core against the published R
  airline-model fit for `log(AirPassengers)` (coefficients, σ², log-likelihood,
  AIC).
- New dependencies: `ffi`, `native_toolchain_c`, `hooks`, `code_assets`
  (build/runtime) and `ffigen` (dev).
- Documented that the vendored ctsa tree is mixed-license (BSD-3-Clause plus
  some LGPL-3.0-or-later units) — to be resolved before mobile distribution.

## 0.1.0

- Initial skeleton of the `tsforge` time-series SDK.
- Pure-Dart `simpleMovingAverage` (SMA) domain primitive.
- Layered structure prepared for the future native (ctsa/FFI) core; no native
  code yet.
