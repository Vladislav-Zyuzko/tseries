/// Pure-Dart damped-trend exponential smoothing (Gardner–McKenzie "damped Holt").
///
/// This is domain-layer code: closed-form, deterministic arithmetic that has no
/// business crossing the FFI boundary. Per the package's design rule, only
/// genuinely hard, validated econometrics (ARIMA/SARIMA and friends) are
/// delegated to the native C core; a recursive smoother whose whole definition
/// is three lines of algebra stays here in idiomatic Dart.
///
/// This is a Tier-1 forecast primitive, and its value is precisely that it is
/// **closed form and deterministic**: there is no parameter optimisation inside.
/// `alpha`, `beta` and `phi` are supplied by the caller (typically swept over a
/// grid) — fitting them per-origin would trade away the auditability that makes
/// a Tier-1 baseline worth having, and optimising a 1-step SSE systematically
/// inflates `alpha` relative to the multi-step loss we actually care about.
///
/// ## Why damped, and why one primitive
///
/// Damping is not an alternative to Holt's linear trend — it is a **superset**.
/// `phi` is an ordinary parameter, so a single primitive yields three baselines
/// and one free correctness test:
///
/// * `phi == 0` → `ŷ = ℓ_t` → plain **SES/EWMA** (does the local trend carry any
///   information at all?);
/// * `0 < phi < 1` → **damped Holt** (the cure for trend overshoot);
/// * `phi == 1` → classical **Holt** linear trend;
/// * `alpha == 1`, any `beta`, `phi == 0` → `ŷ_{t+h} ≡ y_t` → **persistence**,
///   exactly. This configuration must reproduce a persistence baseline bit for
///   bit; if it does not, the implementation is broken.
library;

/// The result of a damped-Holt forecast produced by [holtForecast].
///
/// This is a purely numeric result: it knows nothing about the domain the series
/// comes from. Any domain interpretation (thresholds, alerts, probability of
/// crossing a level) is the caller's job, computed from [point].
class HoltForecast {
  const HoltForecast({required this.point});

  /// Point forecast, one value per step ahead (length == requested horizon).
  ///
  /// `point[h-1] = ℓ_t + b_t · Σ_{i=1..h} phi^i` for `h` in `1..horizon`, where
  /// `ℓ_t` and `b_t` are the level and trend after the recursion has consumed
  /// the whole series.
  ///
  /// Note there is no standard-error companion (unlike `LinearForecast.se`):
  /// exponential smoothing only yields prediction intervals under an explicit
  /// state-space error model, which this closed-form primitive deliberately does
  /// not assume.
  final List<double> point;
}

/// Forecasts [horizon] steps ahead from [series] with the damped-trend
/// exponential smoothing method (Gardner–McKenzie).
///
/// Fully deterministic and free of any native/FFI dependency: the same input
/// always produces a bit-for-bit identical result.
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
/// With `y_0 .. y_{N-1}` the whole of [series]:
///
/// * `ℓ_t = alpha·y_t + (1 − alpha)·(ℓ_{t−1} + phi·b_{t−1})`
/// * `b_t = beta·(ℓ_t − ℓ_{t−1}) + (1 − beta)·phi·b_{t−1}`
/// * `point[h-1] = ℓ_t + b_t · Σ_{i=1..h} phi^i`, `h = 1..H`, `t = N − 1`
///
/// where the damped horizon sum has a closed form that **must** be branched on
/// `phi`:
///
/// * `Σ_{i=1..h} phi^i = phi·(1 − phi^h)/(1 − phi)` for `phi < 1`
/// * `Σ_{i=1..h} phi^i = h` for `phi == 1` — the general expression divides by
///   `1 − phi` and blows up, so this case is computed explicitly.
///
/// For `phi == 0` the sum is 0 and the forecast is flat at `ℓ_t` (i.e. SES).
///
/// ## There is no `window` parameter — on purpose
///
/// Damped Holt has no window hyperparameter: its memory *is* `alpha`, with
/// weights decaying geometrically as `(1 − alpha)^k` (at `alpha = 0.2` a point
/// 144 steps back contributes ~1e-14). The recursion therefore consumes the
/// whole of [series], and callers who want to bound the history slice it
/// themselves before calling — which also lets them hand this primitive exactly
/// the same window as a competing model, so the comparison stays honest.
///
/// ## Initialisation — deliberately not `b_0 = 0`
///
/// The level and trend are seeded by an ordinary-least-squares fit over the
/// first `min(12, N)` points: `ℓ_0` is the fitted value at the first point and
/// `b_0` is the fitted slope. A zero initial slope is *not* acceptable here: at
/// a small `beta` the trend's memory is ~`1/beta` steps (20+ at `beta = 0.05`),
/// so a zero start does not reliably wash out and leaks into the metric.
///
/// ## Parameters
///
/// * [alpha] — level smoothing, in `(0, 1]`.
/// * [beta] — trend smoothing, in `[0, 1]`. `beta == 0` is legitimate (the slope
///   is then never re-estimated from the data, only damped), which is why it is
///   admitted rather than rejected.
/// * [phi] — damping, in `[0, 1]`. See the vocabulary above for the meaning of
///   the two endpoints.
/// * [horizon] (`H`) — number of steps ahead to forecast; must be **≥ 1**.
///
/// ## Failure modes (all [ArgumentError])
///
/// * [horizon] < 1 — nothing to forecast.
/// * `series.length` < 2 — the OLS seed needs at least two points to define a
///   slope. (The floor is 2 rather than `linearForecast`'s 3: that primitive
///   needs a third point only because its residual standard error consumes
///   `N − 2` degrees of freedom, and this one reports no standard error.)
/// * [series] contains a non-finite value (NaN / ±Inf).
/// * [alpha], [beta] or [phi] is non-finite or outside its range. Non-finiteness
///   is checked explicitly: every comparison against NaN is false, so a bare
///   range check would wave NaN straight through.
HoltForecast holtForecast(
  List<double> series, {
  required double alpha,
  required double beta,
  required double phi,
  required int horizon,
}) {
  if (horizon < 1) {
    throw ArgumentError.value(horizon, 'horizon', 'must be at least 1');
  }
  // NaN fails every comparison, so `alpha <= 0 || alpha > 1` would not catch it.
  // Check finiteness first, then the range.
  if (!alpha.isFinite || alpha <= 0 || alpha > 1) {
    throw ArgumentError.value(alpha, 'alpha', 'must be in (0, 1]');
  }
  if (!beta.isFinite || beta < 0 || beta > 1) {
    throw ArgumentError.value(beta, 'beta', 'must be in [0, 1]');
  }
  if (!phi.isFinite || phi < 0 || phi > 1) {
    throw ArgumentError.value(phi, 'phi', 'must be in [0, 1]');
  }
  if (series.length < 2) {
    throw ArgumentError.value(
      series.length,
      'series.length',
      'must be at least 2 (the OLS seed needs two points to define a slope)',
    );
  }

  final n = series.length;
  // Copy the series, rejecting any non-finite value up front (never trust the
  // input blindly): a single NaN would otherwise poison the whole recursion.
  final y = List<double>.filled(n, 0);
  for (var i = 0; i < n; i++) {
    final v = series[i];
    if (!v.isFinite) {
      throw ArgumentError.value(
        series[i],
        'series[$i]',
        'must be finite (no NaN/Inf)',
      );
    }
    y[i] = v;
  }

  // --- Seed: OLS over the first min(12, N) points -------------------------
  // 12 points = 1 h at a 5-minute cadence: enough for a stable slope estimate
  // without letting the seed dominate a long window.
  final m = n < 12 ? n : 12;
  final xBar = (m - 1) / 2.0;
  var yBar = 0.0;
  for (var i = 0; i < m; i++) {
    yBar += y[i];
  }
  yBar /= m;

  var sxx = 0.0; // Σ(x_i − x̄)²
  var sxy = 0.0; // Σ(x_i − x̄)(y_i − ȳ)
  for (var i = 0; i < m; i++) {
    final dx = i - xBar;
    sxx += dx * dx;
    sxy += dx * (y[i] - yBar);
  }

  // Cannot occur for m ≥ 2 distinct regular indices (sxx ≥ 0.5), but guard so a
  // degenerate design never leaks NaN/Inf to the caller.
  if (sxx == 0) {
    throw ArgumentError.value(
      series.length,
      'series.length',
      'degenerate design: zero variance in the seed regression indices',
    );
  }

  final seedSlope = sxy / sxx;
  // ℓ_0 = fitted value at the first point (x = 0) = the intercept itself.
  var level = yBar - seedSlope * xBar;
  var trend = seedSlope;

  // --- Recursion over the whole series -------------------------------------
  for (var t = 0; t < n; t++) {
    final prevLevel = level;
    level = alpha * y[t] + (1 - alpha) * (prevLevel + phi * trend);
    trend = beta * (level - prevLevel) + (1 - beta) * phi * trend;
  }

  // --- Forecast -------------------------------------------------------------
  // point[h-1] = ℓ_t + b_t · Σ_{i=1..h} phi^i.
  // phi^h is accumulated incrementally rather than re-raised per step, and the
  // phi == 1 branch is taken explicitly: the closed form divides by (1 − phi).
  final point = List<double>.filled(horizon, 0);
  var phiPow = 1.0; // phi^(h-1) at the top of the loop.
  for (var h = 1; h <= horizon; h++) {
    phiPow *= phi; // now phi^h
    final dampedSum = phi == 1.0
        ? h.toDouble()
        : phi * (1 - phiPow) / (1 - phi);
    point[h - 1] = level + trend * dampedSum;
  }

  return HoltForecast(point: point);
}
