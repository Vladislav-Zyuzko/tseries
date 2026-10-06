// Comparability of the two native fit paths and the definition of their
// log-likelihood — the preconditions for comparing models by an information
// criterion computed from `logLikelihood` (order selection must put every
// candidate through ONE path, with and without a constant, and compare like
// with like).
//
// Facts pinned here were established on 2026-10-05 against statsmodels 0.15.0
// (exact Gaussian likelihood, stationary initialisation, sigma2 concentrated)
// and closed forms. The pure-AR(1) case was wrong in ctsa until the fix in
// PROVENANCE.md "Local modifications" 10; see ar1_likelihood_test.dart.
import 'dart:io';
import 'dart:math' as math;

import 'package:test/test.dart';
import 'package:tseries/tseries.dart';

List<double> _series(String name) =>
    File('test/fixtures/auto_arima/series_$name.csv')
        .readAsLinesSync()
        .where((l) => l.trim().isNotEmpty)
        .map(double.parse)
        .toList();

/// Exact Gaussian log-likelihood of i.i.d. N(0, σ²) residuals, σ² profiled.
double _iidLogLik(List<double> e) {
  final t = e.length;
  final s = e.fold<double>(0, (a, v) => a + v * v);
  return -t / 2 * (math.log(2 * math.pi) + math.log(s / t) + 1);
}

List<double> _diff(List<double> y) => [
  for (var i = 1; i < y.length; i++) y[i] - y[i - 1],
];

void main() {
  final air = _series('air_passengers_log');
  final nile = _series('nile');
  final lynx = _series('lynx');
  final huron = _series('lake_huron');

  group('fitSarima and fitSarimax (no regressors) are the same fit', () {
    // fitSarima estimates a mean whenever d + D == 0 and has no switch for it,
    // so the like-for-like SARIMAX call is includeMean: true there and
    // includeMean: false (nothing to estimate) when the series is differenced.
    final cases = <(String, List<double>, ArimaOrder)>[
      (
        'airline on log(AirPassengers)',
        air,
        const ArimaOrder(
          p: 0,
          d: 1,
          q: 1,
          seasonalD: 1,
          seasonalQ: 1,
          seasonalPeriod: 12,
        ),
      ),
      ('Nile ARIMA(1,1,1)', nile, const ArimaOrder(p: 1, d: 1, q: 1)),
      ('lynx AR(2) with mean', lynx, const ArimaOrder(p: 2, d: 0, q: 0)),
    ];
    for (final (name, y, order) in cases) {
      test(name, () {
        final a = fitSarima(y, order: order);
        final b = fitSarimax(
          y,
          const {},
          order: order,
          includeMean: order.d + order.seasonalD == 0,
        );
        expect(a.status, ArimaFitStatus.probableSuccess);
        expect(b.status, a.status);
        expect(b.logLikelihood, isNotNull);
        // Same estimate, and both evaluated at it (PROVENANCE.md, Local
        // modifications 13): identical, not merely close.
        expect(b.logLikelihood, a.logLikelihood);
        expect(b.ar, a.ar);
        expect(b.ma, a.ma);
        expect(b.seasonalMa, a.seasonalMa);
        expect(b.sigma2, a.sigma2);
        if (order.d + order.seasonalD == 0) {
          expect(b.intercept, a.mean);
        } else {
          expect(b.intercept, isNull);
        }
      });
    }
  });

  group('logLikelihood is the Gaussian likelihood of the differenced series', () {
    test('random walk: closed form, with -T/2·ln 2π', () {
      final dy = _diff(huron);
      final f = fitSarimax(
        huron,
        const {},
        order: const ArimaOrder(p: 0, d: 1, q: 0),
        includeMean: false,
      );
      expect(f.logLikelihood!, closeTo(_iidLogLik(dy), 1e-9));
      // sigma2 is the profiled MLE S/T, not a degrees-of-freedom corrected one.
      final s = dy.fold<double>(0, (a, v) => a + v * v);
      expect(f.sigma2, closeTo(s / dy.length, 1e-12));
    });

    test(
      'drift is a regression on time with ARIMA errors, same definition',
      () {
        final dy = _diff(huron);
        final f = fitSarimax(
          huron,
          const {},
          order: const ArimaOrder(p: 0, d: 1, q: 0),
          includeDrift: true,
        );
        final mean = dy.reduce((a, b) => a + b) / dy.length;
        final drift = f.coefficient(sarimaxDriftColumn)!;
        expect(drift, closeTo(mean, 1e-12));
        expect(
          f.logLikelihood!,
          // At the returned drift (PROVENANCE.md, Local modifications 13);
          // until then the value came from a Hessian step away and needed 1e-7.
          closeTo(_iidLogLik([for (final v in dy) v - drift]), 1e-10),
        );
      },
    );

    test('pure AR(1): the first observation has variance σ²/(1−φ²)', () {
      // Until the fix recorded as PROVENANCE.md "Local modifications" 10,
      // ctsa's ip == 1 && iq == 0 special case left iupd = 0, so karma's first
      // prediction step overwrote the stationary P0 with V: the first
      // observation got variance σ² and ½·ln(1−φ²) was lost (Nile, this
      // model: 0.09 in log-likelihood). The full matrix is in
      // ar1_likelihood_test.dart; this pins the identity on the path used for
      // order selection.
      final f = fitSarimax(
        nile,
        const {},
        order: const ArimaOrder(p: 1, d: 0, q: 0),
      );
      final phi = f.ar.single;
      final mu = f.intercept!;
      final e = [for (final v in nile) v - mu];
      final exact = <double>[
        math.sqrt(1 - phi * phi) * e.first,
        for (var i = 1; i < e.length; i++) e[i] - phi * e[i - 1],
      ];
      final n = e.length;
      final exactLl = _iidLogLik(exact) + 0.5 * math.log(1 - phi * phi);
      expect(f.logLikelihood!, closeTo(exactLl, 1e-12 * exactLl.abs()));
      final s = exact.fold<double>(0, (a, v) => a + v * v);
      expect(f.sigma2, closeTo(s / n, 1e-12 * s / n));
    });
  });

  test('ModelTooLargeError is raised identically by both paths', () {
    const big = ArimaOrder(p: 0, d: 0, q: 0, seasonalQ: 1, seasonalPeriod: 288);
    final y = List<double>.generate(2000, (i) => math.sin(i / 7) + i % 5);
    expect(() => fitSarima(y, order: big), throwsA(isA<ModelTooLargeError>()));
    expect(
      () => fitSarimax(y, const {}, order: big),
      throwsA(isA<ModelTooLargeError>()),
    );
  });
}
