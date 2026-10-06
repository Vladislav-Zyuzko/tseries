History before the first public release was kept inside the Sweet Limit
monorepo under the name `tsforge` and is not reproduced here.

## 0.8.0 — 2026-10-06 — first public release

A framework-agnostic time-series SDK for Dart (no Flutter dependency).

- **Fixed-order models (native core).** `fitSarima` — ARIMA/SARIMA fit and
  forecast (MLE, CSS, Box-Jenkins) with ctsa's fit status, σ²,
  log-likelihood/AIC (null when not genuinely computed). `fitSarimax` —
  regression with (seasonal) ARIMA errors: named regressors, intercept,
  drift, regressor screening (`ExogPolicy`), scaling of regressors before the
  fit, coefficient SEs and covariance.
- **Forecast ladders.** `arimaForecast` / `sarimaxForecast` — point forecast,
  standard errors and confidence band over an ordered list of candidate
  orders, reporting every attempt.
- **Automatic order selection.** `autoArima` — pure Dart over the native fit,
  after Hyndman & Khandakar (2008): d by successive KPSS tests (`kpssLag`), D
  by the seasonal strength F_S of a classical decomposition, stepwise
  (`ArimaSearch.stepwise`) or exhaustive (`ArimaSearch.exhaustive`) search,
  AICc/AIC/BIC, AR/MA root margin (`rootMargin`), fit budget (`maxModels`).
  Returns the model, forecast and band, a `DifferencingReport` and the full
  search journal; `AutoArimaNoModelException` carries the reason and the
  journal when no model is selected. Written clean-room from the published
  methodology (`doc/auto_arima_clean_room.md`).
- **Pure-Dart primitives.** `kpssLevelTest` (KPSS level stationarity),
  `classicalDecomposition` / `seasonalStrength`, `allRootsOutside` /
  `minRootModulus` (root location of lag polynomials, no root finding),
  `InformationCriterion`, `holtForecast` (damped-trend exponential
  smoothing), `linearForecast` (OLS trend + prediction SE),
  `simpleMovingAverage`.
- **Health check.** `checkNativeCore()` returns `NativeCoreHealth` or throws
  `TseriesNativeUnavailableException` when the native core cannot be loaded.
- **Typed errors.** `TseriesNumericException`, `TseriesInternalException`,
  `ModelTooLargeError` (models whose exact-likelihood state exceeds 64 MiB),
  `SeriesTooShortError`, `ExogenousColumnError`.
- **Native core.** The [ctsa](https://github.com/rafat/ctsa) C library,
  vendored under `third_party/ctsa/` at a pinned commit, reached through a
  thin C shim that exports four functions (`tseries_smoke`,
  `tseries_model_state_bytes`, `tseries_sarima_fit`, `tseries_sarimax_fit`)
  and keeps every ctsa struct on the C side. Bindings are generated with
  ffigen; no `Pointer` reaches the public API. Compiled from source by a Dart
  build hook (`hook/build.dart`, `native_toolchain_c`) without floating-point
  contraction (`-ffp-contract=off`; `/fp:precise` on MSVC), with hidden
  symbol visibility and dead-code elimination.
- **Defects in ctsa fixed or refused.** Memory-safety defects, uninitialised
  reads, likelihood errors and leaks found in the vendored code are patched or
  refused before the native call; each one is documented in
  `third_party/ctsa/PROVENANCE.md` ("Local modifications").
- **Licences.** The package is BSD-3-Clause (`LICENSE`); the notices of the
  vendored code are in `THIRD_PARTY_NOTICES`; the per-file licence audit is
  in PROVENANCE. No LGPL code remains. One question is open: the Applied Statistics algorithms AS 154 /
  AS 182 (with AS 75) on every fit are distributed by the Royal Statistical
  Society "provided that no fee is charged" — see PROVENANCE, "Open licensing
  questions".
- **Verification.** Fits are checked against R and statsmodels reference
  values (airline model, SARIMAX cases, AR(1) likelihood, white-noise
  forecasts); clean-room replacements against scipy/numpy. Per-platform
  results, the FMA gate and how to re-run them are in
  `doc/platform_verification.md`; known limitations of `autoArima` are listed
  in the README, "Limitations of this version".
- **Platforms.** Verified: Windows x64 (MSVC), Android arm64-v8a and x86_64
  (API 36 emulator), Linux x64 (Docker). Built but not run: Android
  armeabi-v7a. Not verified: iOS, macOS. Web is not supported (`dart:ffi`).
