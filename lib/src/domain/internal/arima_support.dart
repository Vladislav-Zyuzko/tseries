/// Helpers shared by the ARIMA and SARIMAX domain APIs.
///
/// Internal to the package: `lib/tseries.dart` does not export this file. It
/// exists so both forecast ladders apply literally the same order validation,
/// numeric length floor and confidence-band arithmetic.
library;

import 'dart:math' as math;
import 'dart:typed_data';

import '../arima.dart' show ArimaMethod, ArimaOrder;

/// Copies [series] into a [Float64List], rejecting any non-finite value with
/// an [ArgumentError] that names [what] and the offending index.
Float64List toFiniteFloat64(List<num> series, {String what = 'series'}) {
  final out = Float64List(series.length);
  for (var i = 0; i < series.length; i++) {
    final v = series[i].toDouble();
    if (!v.isFinite) {
      throw ArgumentError.value(
        series[i],
        '$what[$i]',
        'must be finite (no NaN/Inf)',
      );
    }
    out[i] = v;
  }
  return out;
}

/// Validates the internal consistency of an [ArimaOrder]. A malformed order is
/// a programmer error, so this throws [ArgumentError] (it never contributes to
/// a degradation ladder).
void validateArimaOrder(ArimaOrder o) {
  for (final v in [
    o.p,
    o.d,
    o.q,
    o.seasonalP,
    o.seasonalD,
    o.seasonalQ,
    o.seasonalPeriod,
  ]) {
    if (v < 0) {
      throw ArgumentError.value(v, 'order', 'orders must be non-negative');
    }
  }
  final hasSeasonalTerm = o.seasonalP > 0 || o.seasonalD > 0 || o.seasonalQ > 0;
  if (hasSeasonalTerm && o.seasonalPeriod <= 0) {
    throw ArgumentError.value(
      o.seasonalPeriod,
      'order.seasonalPeriod',
      'must be > 0 when any seasonal term (P/D/Q) is non-zero: $o',
    );
  }
}

/// The numeric minimum series length for an order to be worth fitting.
///
/// This is a NUMERIC floor (below which MLE estimates are not trustworthy), not
/// the methodological minimum a domain might require — the latter is the
/// caller's concern. It requires ≈10 observations per estimated parameter (a
/// customary lower bound for a stable fit), on top of the observations consumed
/// by differencing. The intercept counts as a parameter only for a stationary
/// model (no differencing). [extraParameters] adds further estimated
/// coefficients (regressors), each also at ≈10 observations.
///
/// Examples: ARIMA(2,1,1) → 31, ARIMA(1,1,1) → 21, ARIMA(0,1,1) → 11,
/// airline `(0,1,1)(0,1,1)[12]` → 33. It is always ≥ the native shim's own hard
/// estimability floor, so this friendlier check fires first.
int arimaMinLength(ArimaOrder o, {int extraParameters = 0}) {
  final seasonalSpan = o.seasonalPeriod > 0 ? o.seasonalPeriod : 1;
  final diffConsumed = o.d + o.seasonalD * seasonalSpan;
  final estimatesIntercept = o.d == 0 && o.seasonalD == 0;
  final freeParams =
      o.p +
      o.q +
      o.seasonalP +
      o.seasonalQ +
      (estimatesIntercept ? 1 : 0) +
      extraParameters;
  // ≈10 obs per parameter; keep a small residual-d.o.f. cushion for σ² even for
  // parameter-free orders such as (0,1,0).
  final afterDiff = math.max(10 * freeParams, 8);
  return diffConsumed + afterDiff;
}

/// The smallest series length the native SARIMA fit (`fitSarima`) accepts for
/// [o] with [method], and the rule that sets it — the shim's own size
/// guards (native/tseries_ctsa.c; PROVENANCE.md, Local modifications 12 and
/// 15), mirrored here so the refusal surfaces as a typed
/// `SeriesTooShortError` before the native call. With `Nd = n − d − s·D`:
///
/// * always: `n > d + s·D + p + q + P + Q + 1`;
/// * MLE and CSS: `Nd > p + s·P` — otherwise ctsa's CSS step has no residual
///   left and reads uninitialised memory (results differ run to run);
/// * Box-Jenkins: `Nd ≥ max(p + s·P, q + s·Q)` (out-of-bounds reads), and for a
///   seasonal model (`P + Q > 0`) `Nd ≥ (P + Q + 1)·s + 2` (uninitialised
///   seasonal autocovariances below `(P + Q)·s + 2`, plus one season of
///   margin).
///
/// This is a hard, structural minimum. [arimaMinLength] (≈10 observations per
/// parameter) is the separate numeric floor and is larger for most orders.
({int length, String reason}) sarimaMinimumLength(
  ArimaOrder o,
  ArimaMethod method,
) {
  final s = o.seasonalPeriod > 0 ? o.seasonalPeriod : 0;
  final consumed = o.d + o.seasonalD * s;
  final arLag = o.p + s * o.seasonalP;
  final maLag = o.q + s * o.seasonalQ;
  var length = consumed + o.p + o.q + o.seasonalP + o.seasonalQ + 2;
  var reason = 'n > d + s·D + p + q + P + Q + 1';
  void atLeast(int candidate, String why) {
    if (candidate > length) {
      length = candidate;
      reason = why;
    }
  }

  switch (method) {
    case ArimaMethod.mle:
    case ArimaMethod.css:
      atLeast(
        consumed + arLag + 1,
        'the differenced series must be longer than p + s·P = $arLag',
      );
    case ArimaMethod.boxJenkins:
      atLeast(
        consumed + math.max(arLag, maLag),
        'Box-Jenkins needs at least max(p + s·P, q + s·Q) = '
        '${math.max(arLag, maLag)} differenced observations',
      );
      if (o.seasonalP + o.seasonalQ > 0) {
        final seasonal = (o.seasonalP + o.seasonalQ + 1) * s + 2;
        atLeast(
          consumed + seasonal,
          'seasonal Box-Jenkins needs at least (P + Q + 1)·s + 2 = $seasonal '
          'differenced observations',
        );
      }
  }
  return (length: length, reason: reason);
}

bool allFinite(List<double> xs) {
  for (final x in xs) {
    if (!x.isFinite) return false;
  }
  return true;
}

/// Throws unless [confidenceLevel] is in the open interval (0, 1).
void validateConfidenceLevel(double confidenceLevel) {
  if (!(confidenceLevel > 0 && confidenceLevel < 1)) {
    throw ArgumentError.value(
      confidenceLevel,
      'confidenceLevel',
      'must be in the open interval (0, 1)',
    );
  }
}

/// Two-sided normal critical value `z` such that `P(-z < Z < z) = level`, i.e.
/// `z = Φ⁻¹((1 + level) / 2)`. For level 0.95 this is ≈ 1.959964 (the textbook
/// 1.96). Uses Acklam's rational approximation to the inverse normal CDF
/// (absolute error < 1.2e-9), so any level in (0, 1) is supported.
double normalTwoSidedZ(double level) => _inverseNormalCdf((1 + level) / 2);

double _inverseNormalCdf(double p) {
  const a = <double>[
    -3.969683028665376e+01,
    2.209460984245205e+02,
    -2.759285104469687e+02,
    1.383577518672690e+02,
    -3.066479806614716e+01,
    2.506628277459239e+00,
  ];
  const b = <double>[
    -5.447609879822406e+01,
    1.615858368580409e+02,
    -1.556989798598866e+02,
    6.680131188771972e+01,
    -1.328068155288572e+01,
  ];
  const c = <double>[
    -7.784894002430293e-03,
    -3.223964580411365e-01,
    -2.400758277161838e+00,
    -2.549732539343734e+00,
    4.374664141464968e+00,
    2.938163982698783e+00,
  ];
  const d = <double>[
    7.784695709041462e-03,
    3.224671290700398e-01,
    2.445134137142996e+00,
    3.754408661907416e+00,
  ];
  const pLow = 0.02425;
  const pHigh = 1 - pLow;
  if (p < pLow) {
    final q = math.sqrt(-2 * math.log(p));
    return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q +
            c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  } else if (p <= pHigh) {
    final q = p - 0.5;
    final r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r +
            a[5]) *
        q /
        (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  } else {
    final q = math.sqrt(-2 * math.log(1 - p));
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q +
            c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  }
}
