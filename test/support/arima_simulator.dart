/// A small seeded (S)ARIMA simulator for tests (gate C of the autoArima
/// specification, §9). Deterministic for a given seed within one Dart SDK.
library;

import 'dart:math' as math;

/// Standard normal draws from a seeded [math.Random] (Box–Muller).
class _Normal {
  _Normal(int seed) : _rng = math.Random(seed);
  final math.Random _rng;
  double? _spare;

  double next() {
    final s = _spare;
    if (s != null) {
      _spare = null;
      return s;
    }
    double u;
    do {
      u = _rng.nextDouble();
    } while (u <= 0);
    final v = _rng.nextDouble();
    final r = math.sqrt(-2 * math.log(u));
    _spare = r * math.sin(2 * math.pi * v);
    return r * math.cos(2 * math.pi * v);
  }
}

List<double> _mul(List<double> a, List<double> b) {
  final out = List<double>.filled(a.length + b.length - 1, 0);
  for (var i = 0; i < a.length; i++) {
    for (var j = 0; j < b.length; j++) {
      out[i + j] += a[i] * b[j];
    }
  }
  return out;
}

/// Simulates `n` observations of
/// `φ(B)Φ(B^s)(1−B)^d(1−B^s)^D y_t = c + θ(B)Θ(B^s) e_t`, e_t ~ N(0, σ²),
/// with AR polynomials in the `1 − φB` convention and MA polynomials in the
/// `1 + θB` convention (as R/statsmodels). [constant] is c (with d + D = 1 it
/// is the drift of the differenced series); [mean] is added to the final
/// series (a level for a stationary model). The stationary part starts from
/// zeros and its first [burnIn] values are discarded; integration starts from
/// zero.
List<double> simulateArima({
  required int n,
  required int seed,
  List<double> ar = const [],
  List<double> ma = const [],
  List<double> seasonalAr = const [],
  List<double> seasonalMa = const [],
  int d = 0,
  int seasonalD = 0,
  int period = 1,
  double sigma = 1,
  double constant = 0,
  double mean = 0,
  int burnIn = 200,
}) {
  // Full AR polynomial a(B) = φ(B)Φ(B^s) as 1 − Σ a_j B^j, MA likewise.
  var arPoly = [1.0, for (final v in ar) -v];
  if (seasonalAr.isNotEmpty) {
    final sp = List<double>.filled(seasonalAr.length * period + 1, 0);
    sp[0] = 1;
    for (var i = 0; i < seasonalAr.length; i++) {
      sp[(i + 1) * period] = -seasonalAr[i];
    }
    arPoly = _mul(arPoly, sp);
  }
  var maPoly = [1.0, ...ma];
  if (seasonalMa.isNotEmpty) {
    final sq = List<double>.filled(seasonalMa.length * period + 1, 0);
    sq[0] = 1;
    for (var i = 0; i < seasonalMa.length; i++) {
      sq[(i + 1) * period] = seasonalMa[i];
    }
    maPoly = _mul(maPoly, sq);
  }

  final rng = _Normal(seed);
  final total = n + burnIn;
  final e = [for (var i = 0; i < total; i++) sigma * rng.next()];
  final w = List<double>.filled(total, 0);
  for (var t = 0; t < total; t++) {
    var v = constant + e[t];
    for (var j = 1; j < arPoly.length; j++) {
      if (t - j >= 0) v -= arPoly[j] * w[t - j];
    }
    for (var j = 1; j < maPoly.length; j++) {
      if (t - j >= 0) v += maPoly[j] * e[t - j];
    }
    w[t] = v;
  }
  var x = w.sublist(burnIn);
  for (var k = 0; k < seasonalD; k++) {
    final out = List<double>.filled(x.length, 0);
    for (var t = 0; t < x.length; t++) {
      out[t] = x[t] + (t - period >= 0 ? out[t - period] : 0);
    }
    x = out;
  }
  for (var k = 0; k < d; k++) {
    final out = List<double>.filled(x.length, 0);
    var acc = 0.0;
    for (var t = 0; t < x.length; t++) {
      acc += x[t];
      out[t] = acc;
    }
    x = out;
  }
  return [for (final v in x) v + mean];
}
