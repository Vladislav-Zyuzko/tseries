/// Classical additive seasonal decomposition and the strength of seasonality
/// (pure Dart).
///
/// Sources: the classical decomposition (centred moving-average trend,
/// seasonal indices as per-position means of the detrended series, centred to
/// zero) and the seasonal-strength measure
/// `F_S = max(0, 1 − Var(R) / Var(S + R))` of Hyndman & Athanasopoulos,
/// *Forecasting: Principles and Practice* (3rd ed.), §3.4 and §4.3; the idea
/// of the measure is from Wang, Smith & Hyndman (2006). FPP3 computes F_S on an
/// STL decomposition; here it is computed on the classical one (spec §2.2),
/// whose remainder differs — thresholds calibrated for STL do not transfer
/// automatically.
library;

import 'internal/arima_support.dart';
import 'kpss.dart' show isConstantSeries, maxAbs;

/// The result of [classicalDecomposition]: `y_t = T_t + S_t + R_t`.
class ClassicalDecomposition {
  const ClassicalDecomposition._({
    required this.period,
    required this.trend,
    required this.seasonalIndices,
    required this.seasonal,
    required this.remainder,
    required this.seasonalStrength,
    required this.isDegenerate,
  });

  /// The seasonal period s.
  final int period;

  /// The centred moving-average trend; null at the first and last ⌊s/2⌋
  /// positions, where the moving average is not defined.
  final List<double?> trend;

  /// The s seasonal indices, centred to sum to zero; index `j` belongs to the
  /// positions `t` with `t mod s == j` (t counted from 0 at the first value).
  final List<double> seasonalIndices;

  /// The seasonal component S_t for every position.
  final List<double> seasonal;

  /// The remainder R_t = y_t − T_t − S_t; null where [trend] is null.
  final List<double?> remainder;

  /// F_S = max(0, 1 − Var(R)/Var(S + R)) over the positions where the trend is
  /// defined; in `[0, 1]`. 0 when [isDegenerate].
  final double seasonalStrength;

  /// True when the detrended series S + R is constant by the package's
  /// constancy criterion (`max|x − x̄| ≤ 1e-10·max|y|`, y the input; spec
  /// §2.1/§2.2): the series has no variation around its trend, F_S is
  /// undefined and reported as 0.
  final bool isDegenerate;
}

/// Classical additive decomposition of [series] with seasonal [period]
/// (FPP3 §3.4):
///
/// 1. trend — a centred moving average of order s (a 2×s moving average when
///    s is even);
/// 2. seasonal indices — the mean of the detrended series y − T at each
///    position in the season, over the positions where T is defined, then
///    centred to sum to zero;
/// 3. remainder R = y − T − S where T is defined;
/// 4. [ClassicalDecomposition.seasonalStrength] — F_S (FPP3 §4.3).
///
/// Variances are sample variances; the ratio does not depend on the
/// denominator convention.
///
/// Throws [ArgumentError] for a non-finite value, `period < 2`, or fewer than
/// `2·period` observations.
ClassicalDecomposition classicalDecomposition(
  List<num> series, {
  required int period,
}) {
  final y = toFiniteFloat64(series);
  final s = period;
  if (s < 2) {
    throw ArgumentError.value(period, 'period', 'must be at least 2');
  }
  final n = y.length;
  if (n < 2 * s) {
    throw ArgumentError.value(
      n,
      'series.length',
      'at least 2·period = ${2 * s} observations are required',
    );
  }

  // 1. Centred moving average.
  final half = s ~/ 2;
  final trend = List<double?>.filled(n, null);
  for (var t = half; t < n - half; t++) {
    var acc = 0.0;
    if (s.isOdd) {
      for (var i = t - half; i <= t + half; i++) {
        acc += y[i];
      }
    } else {
      acc = 0.5 * y[t - half] + 0.5 * y[t + half];
      for (var i = t - half + 1; i <= t + half - 1; i++) {
        acc += y[i];
      }
    }
    trend[t] = acc / s;
  }

  // 2. Seasonal indices.
  final sums = List<double>.filled(s, 0);
  final counts = List<int>.filled(s, 0);
  for (var t = 0; t < n; t++) {
    final tr = trend[t];
    if (tr == null) continue;
    sums[t % s] += y[t] - tr;
    counts[t % s]++;
  }
  final raw = [for (var j = 0; j < s; j++) sums[j] / counts[j]];
  var rawMean = 0.0;
  for (final v in raw) {
    rawMean += v;
  }
  rawMean /= s;
  final indices = List<double>.unmodifiable([for (final v in raw) v - rawMean]);
  final seasonal = List<double>.unmodifiable([
    for (var t = 0; t < n; t++) indices[t % s],
  ]);

  // 3. Remainder, and 4. strength over the positions with a trend.
  final remainder = List<double?>.filled(n, null);
  final r = <double>[];
  final sr = <double>[];
  for (var t = 0; t < n; t++) {
    final tr = trend[t];
    if (tr == null) continue;
    final detrended = y[t] - tr;
    remainder[t] = detrended - seasonal[t];
    r.add(detrended - seasonal[t]);
    sr.add(detrended);
  }

  final varSr = _variance(sr);
  final degenerate = isConstantSeries(sr, maxAbs(y)) || !(varSr > 0);
  final strength = degenerate
      ? 0.0
      : (1 - _variance(r) / varSr).clamp(0.0, 1.0).toDouble();

  return ClassicalDecomposition._(
    period: s,
    trend: List<double?>.unmodifiable(trend),
    seasonalIndices: indices,
    seasonal: seasonal,
    remainder: List<double?>.unmodifiable(remainder),
    seasonalStrength: strength,
    isDegenerate: degenerate,
  );
}

/// The strength of seasonality F_S of [series] at [period] on the classical
/// additive decomposition — shorthand for
/// `classicalDecomposition(series, period: period).seasonalStrength`.
double seasonalStrength(List<num> series, {required int period}) =>
    classicalDecomposition(series, period: period).seasonalStrength;

double _variance(List<double> xs) {
  final n = xs.length;
  if (n < 2) return 0;
  var mean = 0.0;
  for (final v in xs) {
    mean += v;
  }
  mean /= n;
  var acc = 0.0;
  for (final v in xs) {
    acc += (v - mean) * (v - mean);
  }
  return acc / (n - 1);
}
