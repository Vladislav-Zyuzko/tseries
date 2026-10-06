/// The KPSS test of level stationarity (pure Dart).
///
/// Source: D. Kwiatkowski, P. C. B. Phillips, P. Schmidt and Y. Shin (1992),
/// "Testing the null hypothesis of stationarity against the alternative of a
/// unit root", *Journal of Econometrics* 54, 159–178 — the level-stationarity
/// statistic η_μ, its Table 1 critical values and the lag rules "l4"/"l12".
/// The default lag rule `short` = ⌊3·√T/13⌋ is the one that reproduces the
/// differencing decisions of R's `forecast::ndiffs` (spec §14.1).
///
/// The null hypothesis is **stationarity**: a large statistic is evidence of a
/// unit root. That is why automatic order selection uses this test to choose
/// the number of differences (Hyndman & Khandakar 2008, §3.1): a test whose
/// null is stationarity does not push towards needless differencing.
library;

import 'internal/arima_support.dart';

/// The significance levels for which Table 1 of KPSS (1992) gives critical
/// values of the level-stationarity statistic η_μ, mapped to those values.
///
/// Only these levels can be used: p-values are never interpolated or
/// extrapolated beyond the table.
final Map<double, double> kpssLevelCriticalValues = Map.unmodifiable({
  0.10: 0.347,
  0.05: 0.463,
  0.025: 0.574,
  0.01: 0.739,
});

/// Where the p-value of a KPSS statistic lies relative to the levels of
/// [kpssLevelCriticalValues]. The p-value itself is not computed: Table 1 only
/// brackets it.
enum KpssPValueBracket {
  /// η ≥ 0.739: p < 0.01.
  below001('<0.01'),

  /// 0.574 ≤ η < 0.739: 0.01 < p ≤ 0.025.
  from001To0025('0.01–0.025'),

  /// 0.463 ≤ η < 0.574: 0.025 < p ≤ 0.05.
  from0025To005('0.025–0.05'),

  /// 0.347 ≤ η < 0.463: 0.05 < p ≤ 0.10.
  from005To010('0.05–0.10'),

  /// η < 0.347: p > 0.10.
  above010('>0.10');

  const KpssPValueBracket(this.label);

  /// Human-readable form, e.g. `0.025–0.05`.
  final String label;

  /// The bracket of a statistic [eta].
  static KpssPValueBracket of(double eta) {
    if (eta >= kpssLevelCriticalValues[0.01]!) return below001;
    if (eta >= kpssLevelCriticalValues[0.025]!) return from001To0025;
    if (eta >= kpssLevelCriticalValues[0.05]!) return from0025To005;
    if (eta >= kpssLevelCriticalValues[0.10]!) return from005To010;
    return above010;
  }
}

/// How the Bartlett-window lag l of the KPSS long-run variance is chosen from
/// the number of observations T of the series being tested.
///
/// In a sequential procedure (`autoArima`) the rule is applied afresh at every
/// step, to the T of the series at that step (spec §2.1, revision 6).
final class KpssLag {
  const KpssLag._(this.name, this._fixed);

  /// `l = ⌊3·√T / 13⌋` (0 for T < 19) — the default. The rule that
  /// reproduces the differencing decisions of R's `forecast::ndiffs` (KPSS,
  /// level) on every series probed (spec §14.1). A short lag: more power
  /// against unit roots, a larger size on persistent stationary series than
  /// [l4].
  static const KpssLag short = KpssLag._('short', null);

  /// `l = ⌊4·(T/100)^{1/4}⌋` — KPSS (1992) "l4".
  static const KpssLag l4 = KpssLag._('l4', null);

  /// `l = ⌊12·(T/100)^{1/4}⌋` — KPSS (1992) "l12".
  static const KpssLag l12 = KpssLag._('l12', null);

  /// A fixed lag [lag] whatever T is.
  const KpssLag.fixed(int lag) : this._('fixed', lag);

  /// The rule's name: `short`, `l4`, `l12` or `fixed`.
  final String name;
  final int? _fixed;

  /// The lag this rule gives for [length] = T observations. Computed in exact
  /// integer arithmetic (`169·l² ≤ 9·T` for [short], `100·l⁴ ≤ c⁴·T` for l4
  /// and l12), so no floating-point rounding moves a boundary.
  int lagFor(int length) {
    if (length < 0) {
      throw ArgumentError.value(length, 'length', 'must be non-negative');
    }
    final f = _fixed;
    if (f != null) return f;
    bool fits(int l) => switch (name) {
      'short' => 169 * l * l <= 9 * length,
      'l4' => 100 * l * l * l * l <= 256 * length,
      _ => 100 * l * l * l * l <= 20736 * length,
    };
    var l = 0;
    while (fits(l + 1)) {
      l++;
    }
    return l;
  }

  /// Text form for reports: `short`, `l4`, `l12`, `fixed(3)`.
  String get description => _fixed == null ? name : 'fixed($_fixed)';

  @override
  bool operator ==(Object other) =>
      other is KpssLag && other.name == name && other._fixed == _fixed;

  @override
  int get hashCode => Object.hash(name, _fixed);

  @override
  String toString() => 'KpssLag.$description';
}

/// The outcome of [kpssLevelTest].
class KpssResult {
  const KpssResult({
    required this.length,
    required this.lagRule,
    required this.lag,
    required this.statistic,
    required this.longRunVariance,
    required this.alpha,
    required this.criticalValue,
  });

  /// Number of observations T the statistic was computed on.
  final int length;

  /// The rule that chose [lag].
  final KpssLag lagRule;

  /// Bartlett-window lag l of the long-run variance.
  final int lag;

  /// The statistic η_μ.
  final double statistic;

  /// The long-run variance estimate s²(l) (its denominator).
  final double longRunVariance;

  /// The significance level of the decision.
  final double alpha;

  /// The KPSS (1992) Table 1 critical value at [alpha].
  final double criticalValue;

  /// Table bracket of the p-value.
  KpssPValueBracket get pValue => KpssPValueBracket.of(statistic);

  /// Whether stationarity is rejected at [alpha], i.e. η > [criticalValue]:
  /// the series should be differenced.
  bool get rejectsStationarity => statistic > criticalValue;

  @override
  String toString() =>
      'KPSS(T=$length, l=$lag ${lagRule.description}): eta=$statistic, '
      'crit($alpha)=$criticalValue, p ${pValue.label} -> '
      '${rejectsStationarity ? 'difference' : 'stationary'}';
}

/// Relative tolerance of the constancy criterion shared by KPSS, the seasonal
/// strength and `autoArima` (spec §2.1, revision 5): a series x is constant
/// when `max|x − x̄| ≤ 1e-10 · scale`, where scale = max|y| of the original
/// input (and always when that scale is 0).
const double constantSeriesTolerance = 1e-10;

/// Whether [x] is constant relative to [scale] (see
/// [constantSeriesTolerance]).
bool isConstantSeries(List<double> x, double scale) {
  if (x.isEmpty || !(scale > 0)) return true;
  var mean = 0.0;
  for (final v in x) {
    mean += v;
  }
  mean /= x.length;
  var dev = 0.0;
  for (final v in x) {
    final a = (v - mean).abs();
    if (a > dev) dev = a;
  }
  return dev <= constantSeriesTolerance * scale;
}

/// max|y| over [y] (0 for an empty list).
double maxAbs(List<double> y) {
  var m = 0.0;
  for (final v in y) {
    if (v.abs() > m) m = v.abs();
  }
  return m;
}

/// The KPSS test of **level** stationarity of [series] (KPSS 1992).
///
/// With e_t = x_t − x̄ and S_t = e_1 + … + e_t:
///
/// ```text
/// s²(l) = T⁻¹ Σ e_t² + 2 T⁻¹ Σ_{j=1..l} (1 − j/(l+1)) Σ_{t=j+1..T} e_t e_{t−j}
/// η_μ   = T⁻² Σ S_t² / s²(l)
/// ```
///
/// [lag] chooses l from T (default [KpssLag.short]); [alpha] must be one of
/// the levels of [kpssLevelCriticalValues] (default 5 %). The decision always
/// compares η with the Table 1 critical value — also at α = 1 %, where R's
/// `ndiffs` (whose p-value is bounded below) never differences: this package
/// deliberately differs there.
///
/// Throws [ArgumentError] for a non-finite value, fewer than 2 observations, a
/// lag outside `[0, T − 1]`, an [alpha] not in the table, or a constant series
/// (`max|x − x̄| ≤ 1e-10·max|x|`): η is undefined there.
KpssResult kpssLevelTest(
  List<num> series, {
  KpssLag lag = KpssLag.short,
  double alpha = 0.05,
}) {
  final x = toFiniteFloat64(series);
  final crit = kpssLevelCriticalValues[alpha];
  if (crit == null) {
    throw ArgumentError.value(
      alpha,
      'alpha',
      'must be one of ${kpssLevelCriticalValues.keys.join(', ')} '
          '(KPSS 1992, Table 1)',
    );
  }
  final t = x.length;
  if (t < 2) {
    throw ArgumentError.value(t, 'series.length', 'must be at least 2');
  }
  final l = lag.lagFor(t);
  if (l < 0 || l >= t) {
    throw ArgumentError.value(l, 'lag', 'must be in [0, ${t - 1}]');
  }
  final stat = kpssStatistic(x, l, scale: maxAbs(x));
  if (stat == null) {
    throw ArgumentError.value(
      series,
      'series',
      'constant series (or non-positive long-run variance): the KPSS '
          'statistic is undefined',
    );
  }
  return KpssResult(
    length: t,
    lagRule: lag,
    lag: l,
    statistic: stat.eta,
    longRunVariance: stat.s2,
    alpha: alpha,
    criticalValue: crit,
  );
}

/// η_μ and s²(l) of [x] with lag [l], or null when [x] is constant relative
/// to [scale] or s²(l) is not a positive finite number (η undefined).
/// Internal: [x] is assumed finite.
({double eta, double s2})? kpssStatistic(
  List<double> x,
  int l, {
  required double scale,
}) {
  if (isConstantSeries(x, scale)) return null;
  final t = x.length;
  var mean = 0.0;
  for (final v in x) {
    mean += v;
  }
  mean /= t;
  final e = [for (final v in x) v - mean];

  var s2 = 0.0;
  for (final v in e) {
    s2 += v * v;
  }
  for (var j = 1; j <= l; j++) {
    var acc = 0.0;
    for (var i = j; i < t; i++) {
      acc += e[i] * e[i - j];
    }
    s2 += 2 * (1 - j / (l + 1)) * acc;
  }
  s2 /= t;
  if (!(s2 > 0) || !s2.isFinite) return null;

  var partial = 0.0;
  var sumSq = 0.0;
  for (final v in e) {
    partial += v;
    sumSq += partial * partial;
  }
  final eta = sumSq / (t.toDouble() * t) / s2;
  if (!eta.isFinite) return null;
  return (eta: eta, s2: s2);
}
