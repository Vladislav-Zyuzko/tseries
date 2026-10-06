/// Root-location checks for ARMA lag polynomials without finding the roots
/// (pure Dart).
///
/// A lag polynomial is written `a(z) = 1 − a_1 z − … − a_m z^m`. For an AR part
/// `φ(B) = 1 − φ_1 B − …` the coefficients are `a_j = φ_j`; for an MA part in
/// the `1 + θ_1 B + …` convention they are `a_j = −θ_j` (in ctsa's
/// `1 − θ_1 B − …` convention, used by this package's fit results, they are
/// `a_j = θ_j`).
///
/// Method (spec §5, rule 5): all roots of a(z) lie outside the circle
/// `|z| = ρ` iff all roots of `1 − Σ a_j ρ^j w^j` lie outside the unit circle,
/// which the step-down (reverse Levinson–Durbin) recursion decides: the
/// polynomial is stable iff every reflection coefficient has modulus < 1
/// (the Schur–Cohn criterion). Hyndman & Khandakar (2008, §3.2) reject a
/// candidate whose AR or MA polynomial has a root of modulus below 1.001;
/// `autoArima` defaults to 1.01 (spec §5, revision 9).
library;

import 'dart:math' as math;

/// Whether every root of `a(z) = 1 − Σ_{j=1..m} a_j z^j` has modulus
/// strictly greater than [radius] (default 1: the polynomial is stable /
/// invertible). A polynomial of degree 0 (no coefficients, or only zeros)
/// has no roots and passes.
///
/// [coefficients] are `a_1..a_m`. Throws [ArgumentError] for a non-finite
/// coefficient or a non-positive [radius].
bool allRootsOutside(List<double> coefficients, {double radius = 1}) {
  if (!(radius > 0) || !radius.isFinite) {
    throw ArgumentError.value(radius, 'radius', 'must be positive and finite');
  }
  for (final c in coefficients) {
    if (!c.isFinite) {
      throw ArgumentError.value(coefficients, 'coefficients', 'must be finite');
    }
  }
  return _stepDownStable(_trim(coefficients), radius);
}

/// The smallest modulus of the roots of `a(z) = 1 − Σ a_j z^j`, or null when
/// the polynomial has degree 0 (no roots).
///
/// Exact closed forms for degree 1 and 2; for higher degrees a bisection on
/// the radius of [allRootsOutside], to a relative width of 1e-12 — an
/// approximation whose accuracy is limited by the conditioning of the
/// step-down recursion near the boundary.
double? minRootModulus(List<double> coefficients) {
  final a = _trim(coefficients);
  for (final c in a) {
    if (!c.isFinite) {
      throw ArgumentError.value(coefficients, 'coefficients', 'must be finite');
    }
  }
  final m = a.length;
  if (m == 0) return null;
  if (m == 1) return 1 / a[0].abs();
  if (m == 2) {
    // The inverse roots u = 1/z solve u² − a_1 u − a_2 = 0.
    final a1 = a[0];
    final a2 = a[1];
    final disc = a1 * a1 + 4 * a2;
    final maxInverse = disc >= 0
        ? (a1.abs() + math.sqrt(disc)) / 2
        : math.sqrt(-a2);
    return 1 / maxInverse;
  }
  // The product of the root moduli is 1/|a_m|, so the smallest is at most
  // |a_m|^(−1/m); doubling that bracket guarantees an unstable upper end.
  var hi = 2 * math.pow(a[m - 1].abs(), -1 / m).toDouble();
  while (_stepDownStable(a, hi)) {
    hi *= 2;
  }
  var lo = 0.0;
  while (hi - lo > 1e-12 * hi) {
    final mid = 0.5 * (lo + hi);
    if (mid <= lo || mid >= hi) break;
    if (_stepDownStable(a, mid)) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return 0.5 * (lo + hi);
}

List<double> _trim(List<double> a) {
  var m = a.length;
  while (m > 0 && a[m - 1] == 0) {
    m--;
  }
  return a.sublist(0, m);
}

/// Step-down recursion on `a_j · ρ^j`: true iff every root lies outside |z|=ρ.
///
/// Carried out in double-double arithmetic (≈ 32 significant digits): near a
/// cluster of roots at the boundary the reflection coefficients approach ±1
/// and the update `a_j + κ·a_{k−j}` cancels catastrophically, so plain double
/// precision can flip the verdict for roots ~1e-6 away from ρ (found by the
/// gate-A comparison against high-precision roots; see
/// test/auto_arima_components_test.dart).
bool _stepDownStable(List<double> a, double rho) {
  final m = a.length;
  if (m == 0) return true;
  final r = _DD(rho, 0);
  var power = _DD.one;
  var cur = List<_DD>.generate(m, (j) {
    power = power * r;
    return power.scale(a[j]);
  });
  for (var k = m; k >= 1; k--) {
    final kappa = cur[k - 1];
    if (!(kappa.abs() < _DD.one)) return false;
    // 1 − κ² as (1 − κ)(1 + κ): both factors are exact-ish near |κ| = 1.
    final denom = (_DD.one - kappa) * (_DD.one + kappa);
    final next = List<_DD>.generate(
      k - 1,
      (i) => (cur[i] + kappa * cur[k - 2 - i]) / denom,
    );
    cur = next;
  }
  return true;
}

/// A minimal double-double number hi + lo (|lo| ≤ ulp(hi)/2), with the
/// classical error-free transformations (Dekker 1971; Knuth's TwoSum).
class _DD {
  const _DD(this.hi, this.lo);

  final double hi;
  final double lo;

  static const one = _DD(1, 0);
  static const double _split = 134217729; // 2^27 + 1

  static _DD _quickTwoSum(double a, double b) {
    final s = a + b;
    return _DD(s, b - (s - a));
  }

  static _DD _twoSum(double a, double b) {
    final s = a + b;
    final bb = s - a;
    return _DD(s, (a - (s - bb)) + (b - bb));
  }

  static _DD _twoProd(double a, double b) {
    final p = a * b;
    var t = _split * a;
    final ah = t - (t - a);
    final al = a - ah;
    t = _split * b;
    final bh = t - (t - b);
    final bl = b - bh;
    return _DD(p, ((ah * bh - p) + ah * bl + al * bh) + al * bl);
  }

  _DD operator +(_DD o) {
    final s = _twoSum(hi, o.hi);
    final t = _twoSum(lo, o.lo);
    var u = _quickTwoSum(s.hi, s.lo + t.hi);
    u = _quickTwoSum(u.hi, u.lo + t.lo);
    return u;
  }

  _DD operator -() => _DD(-hi, -lo);

  _DD operator -(_DD o) => this + (-o);

  _DD operator *(_DD o) {
    final p = _twoProd(hi, o.hi);
    return _quickTwoSum(p.hi, p.lo + (hi * o.lo + lo * o.hi));
  }

  _DD scale(double b) => this * _DD(b, 0);

  _DD operator /(_DD o) {
    final q1 = hi / o.hi;
    var r = this - o.scale(q1);
    final q2 = r.hi / o.hi;
    r = r - o.scale(q2);
    final q3 = r.hi / o.hi;
    final q = _quickTwoSum(q1, q2);
    return q + _DD(q3, 0);
  }

  _DD abs() => hi < 0 || (hi == 0 && lo < 0) ? -this : this;

  bool operator <(_DD o) => hi < o.hi || (hi == o.hi && lo < o.lo);
}
