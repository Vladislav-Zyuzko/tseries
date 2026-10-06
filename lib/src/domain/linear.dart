/// Pure-Dart ordinary-least-squares (OLS) linear forecasting.
///
/// This is domain-layer code: closed-form, deterministic arithmetic that has no
/// business crossing the FFI boundary. Per the package's design rule, only
/// genuinely hard, validated econometrics (ARIMA/SARIMA and friends) are
/// delegated to the native C core; a simple linear trend fit stays here in
/// idiomatic Dart.
///
/// This is Tier 1 of the forecasting methodology: fit a straight line by OLS to
/// the last `window` points of an already-regular series and extrapolate it over
/// a horizon. No iteration, no native dependency — just the normal equations.
library;

import 'dart:math' as math;

/// The result of a linear (OLS) forecast produced by [linearForecast].
///
/// Carries the point forecast, its per-step standard errors, and the fitted line
/// (`slope`, `intercept`). All lists have length equal to the requested horizon.
///
/// This is a purely numeric result: it knows nothing about the domain the series
/// comes from. Any domain interpretation (thresholds, alerts, probability of
/// crossing a level) is the caller's job, computed from [point] and [se].
class LinearForecast {
  const LinearForecast({
    required this.point,
    required this.se,
    required this.slope,
    required this.intercept,
  });

  /// Point forecast, one value per step ahead (length == requested horizon).
  ///
  /// `point[h-1] = intercept + slope * (x_last + h)` for `h` in `1..horizon`,
  /// where `x_last = window - 1` is the index of the most recent point in the
  /// fitting window.
  final List<double> point;

  /// Standard error of each forecast step (same length as [point]).
  ///
  /// This is the standard error of prediction for a new observation:
  /// `se[h] = s · √(1 + 1/N + (x_h − x̄)² / Sxx)`, so it grows the further the
  /// forecast index `x_h` sits from the centre of the fitting window. `s` is the
  /// residual standard deviation and `N` is `window`. When the fit is exact
  /// (all residuals zero) `s = 0` and every `se` is 0.
  final List<double> se;

  /// Slope `b` of the fitted line (change in value per one step).
  final double slope;

  /// Intercept `a` of the fitted line at the window's local origin (`x = 0`,
  /// i.e. the oldest point of the fitting window).
  final double intercept;
}

/// Fits an ordinary-least-squares (OLS) trend line to the last [window] points
/// of [series] and forecasts [horizon] steps ahead.
///
/// This is the Tier-1 forecast primitive: a closed-form linear extrapolation,
/// fully deterministic and free of any native/FFI dependency.
///
/// ## Preconditions (the package does NOT relax these)
///
/// * [series] must already be a **regular** series — evenly spaced, with no gaps
///   or missing samples. Resampling and imputation are the caller's job, exactly
///   as for `arimaForecast`; this package is time-agnostic and never invents
///   data points.
/// * [horizon] is measured in **steps**, not in time. The caller owns the
///   sampling interval Δt.
///
/// ## Model
///
/// Let the last `N = window` points be `y_0 .. y_{N-1}` at local indices
/// `x_i = i` (`0 .. N-1`). With `x̄`, `ȳ` the means and
/// `Sxx = Σ(x_i − x̄)²`:
///
/// * `slope  b = Σ(x_i − x̄)(y_i − ȳ) / Sxx`
/// * `intercept a = ȳ − b·x̄`
/// * `point[h] = a + b·(x_last + h)`, `h = 1..H`, `x_last = N − 1`
/// * residual std dev `s = √(Σ resid² / (N − 2))`, `resid_i = y_i − (a + b·x_i)`
/// * `se[h] = s·√(1 + 1/N + (x_last + h − x̄)² / Sxx)`
///
/// ## Parameters
///
/// * [window] (`N`) — how many of the most recent points enter the regression
///   (a hyperparameter; default 4). Must be **≥ 3**: the residual standard
///   error uses `N − 2` degrees of freedom, so it is undefined for fewer points.
///   [series] must have at least [window] points.
/// * [horizon] (`H`) — number of steps ahead to forecast; must be **≥ 1**.
///
/// ## Failure modes (all [ArgumentError])
///
/// * [window] < 3 — the residual std error is not defined.
/// * [horizon] < 1 — nothing to forecast.
/// * `series.length` < [window] — not enough points to fill the window.
/// * [series] contains a non-finite value (NaN / ±Inf).
/// * degenerate design (`Sxx == 0`) — cannot happen for `N ≥ 2` regular
///   indices, but guarded defensively so a pathological input surfaces cleanly
///   instead of producing NaN/Inf.
LinearForecast linearForecast(
  List<num> series, {
  int window = 4,
  required int horizon,
}) {
  if (window < 3) {
    throw ArgumentError.value(
      window,
      'window',
      'must be at least 3 (residual std error uses N-2 degrees of freedom)',
    );
  }
  if (horizon < 1) {
    throw ArgumentError.value(horizon, 'horizon', 'must be at least 1');
  }
  if (series.length < window) {
    throw ArgumentError.value(
      series.length,
      'series.length',
      'must be at least window ($window)',
    );
  }

  final n = window;
  // Take exactly the last `window` points, promoting to double and rejecting
  // any non-finite value up front (never trust the input blindly).
  final y = List<double>.filled(n, 0);
  final start = series.length - n;
  for (var i = 0; i < n; i++) {
    final v = series[start + i].toDouble();
    if (!v.isFinite) {
      throw ArgumentError.value(
        series[start + i],
        'series[${start + i}]',
        'must be finite (no NaN/Inf)',
      );
    }
    y[i] = v;
  }

  // Local indices x_i = 0..N-1. Means and the classic OLS sums.
  final xBar = (n - 1) / 2.0;
  var yBar = 0.0;
  for (var i = 0; i < n; i++) {
    yBar += y[i];
  }
  yBar /= n;

  var sxx = 0.0; // Σ(x_i − x̄)²
  var sxy = 0.0; // Σ(x_i − x̄)(y_i − ȳ)
  for (var i = 0; i < n; i++) {
    final dx = i - xBar;
    sxx += dx * dx;
    sxy += dx * (y[i] - yBar);
  }

  // Cannot occur for N ≥ 2 distinct regular indices, but guard so a degenerate
  // design never leaks NaN/Inf to the caller.
  if (sxx == 0) {
    throw ArgumentError.value(
      window,
      'window',
      'degenerate design: zero variance in the regression indices',
    );
  }

  final slope = sxy / sxx;
  final intercept = yBar - slope * xBar;

  // Residual sum of squares → residual std dev with N-2 d.o.f.
  var rss = 0.0;
  for (var i = 0; i < n; i++) {
    final resid = y[i] - (intercept + slope * i);
    rss += resid * resid;
  }
  final s = math.sqrt(rss / (n - 2));

  final xLast = (n - 1).toDouble();
  final point = List<double>.filled(horizon, 0);
  final se = List<double>.filled(horizon, 0);
  for (var h = 1; h <= horizon; h++) {
    final xh = xLast + h;
    point[h - 1] = intercept + slope * xh;
    final dx = xh - xBar;
    se[h - 1] = s * math.sqrt(1 + 1 / n + (dx * dx) / sxx);
  }

  return LinearForecast(
    point: point,
    se: se,
    slope: slope,
    intercept: intercept,
  );
}
