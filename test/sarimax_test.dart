/// Behavioural coverage of the SARIMAX path: equivalence with the SARIMA path,
/// regressor screening under both policies, the native guards (degrees of
/// freedom, infeasible seasonality, non-finite input), the trust gates, and
/// the forecast ladder. Numerical agreement with statsmodels is in
/// `sarimax_golden_test.dart`.
library;

import 'dart:math' as math;

import 'package:test/test.dart';
import 'package:tseries/src/ffi/bindings/tseries_ctsa_bindings.g.dart' as b;
import 'package:tseries/tseries.dart';

import 'fixtures/air_passengers.dart';
import 'fixtures/sarimax_reference.dart';

/// Deterministic standard-normal stream (Box-Muller over dart:math Random).
class _Gauss {
  _Gauss(int seed) : _r = math.Random(seed);
  final math.Random _r;
  double next() {
    final u1 = _r.nextDouble() + 1e-12;
    final u2 = _r.nextDouble();
    return math.sqrt(-2 * math.log(u1)) * math.cos(2 * math.pi * u2);
  }
}

/// ARIMA(1,1,1)-error series with two informative regressors (beta 0.8, -0.5).
({List<double> y, List<double> x1, List<double> x2}) _data({
  int n = 120,
  int seed = 3,
}) {
  final g = _Gauss(seed);
  final x1 = [for (var t = 0; t < n; t++) math.sin(t / 4) * 3 + g.next()];
  final x2 = [for (var t = 0; t < n; t++) math.cos(t / 9) * 2 + g.next()];
  final y = <double>[];
  var level = 5.0, w = 0.0, ePrev = 0.0;
  for (var t = 0; t < n; t++) {
    final e = 0.3 * g.next();
    w = 0.5 * w + e + 0.3 * ePrev;
    ePrev = e;
    level += w;
    y.add(level + 0.8 * x1[t] - 0.5 * x2[t]);
  }
  return (y: y, x1: x1, x2: x2);
}

const _arima111 = ArimaOrder(p: 1, d: 1, q: 1);

void main() {
  group('equivalence with fitSarima when there are no regressors', () {
    test('airline model on log(AirPassengers): identical fit and forecast', () {
      final logAir = airPassengers.map((v) => math.log(v)).toList();
      const order = ArimaOrder(
        p: 0,
        d: 1,
        q: 1,
        seasonalP: 0,
        seasonalD: 1,
        seasonalQ: 1,
        seasonalPeriod: 12,
      );
      final a = fitSarima(logAir, order: order, horizon: 12);
      final x = fitSarimax(logAir, const {}, order: order, horizon: 12);
      expect(x.ma, a.ma);
      expect(x.seasonalMa, a.seasonalMa);
      expect(x.sigma2, a.sigma2);
      expect(x.logLikelihood, a.logLikelihood);
      expect(x.aic, a.aic);
      expect(x.forecast, a.forecast);
      expect(x.standardErrors, a.standardErrors);
      expect(x.retval, a.retval);
      expect(x.regressors, isEmpty);
    });

    test('ARMA(2,0,1) with a mean (d = 0): identical, intercept == mean', () {
      final y = sarimaxRefStationaryWithMean.y;
      const order = ArimaOrder(p: 2, d: 0, q: 1);
      final a = fitSarima(y, order: order, horizon: 6);
      final x = fitSarimax(y, const {}, order: order, horizon: 6);
      expect(x.ar, a.ar);
      expect(x.ma, a.ma);
      expect(x.intercept, a.mean);
      expect(x.sigma2, a.sigma2);
      expect(x.forecast, a.forecast);
      expect(x.standardErrors, a.standardErrors);
    });
  });

  group('intercept', () {
    test('includeMean: false with d = 0 estimates no intercept', () {
      final c = sarimaxRefStationaryWithMean;
      final f = fitSarimax(
        c.y,
        c.exog,
        order: const ArimaOrder(p: 1, d: 0, q: 1),
        includeMean: false,
      );
      expect(f.intercept, isNull);
      expect(f.interceptStandardError, isNull);
      expect(f.covariance!.parameters, isNot(contains('intercept')));
    });

    test('with differencing the intercept is never estimated', () {
      final d = _data();
      final f = fitSarimax(d.y, {'x1': d.x1}, order: _arima111);
      expect(f.intercept, isNull);
    });
  });

  group('regressor screening — strict policy (default)', () {
    final d = _data();

    ExogenousColumnError strictError(
      Map<String, List<num>> exog, {
      ArimaOrder order = _arima111,
      bool includeMean = true,
    }) {
      try {
        fitSarimax(d.y, exog, order: order, includeMean: includeMean);
      } on ExogenousColumnError catch (e) {
        return e;
      }
      fail('expected ExogenousColumnError');
    }

    test('an all-zero column is named', () {
      final e = strictError({
        'x1': d.x1,
        'zeros': List.filled(d.y.length, 0),
        'x2': d.x2,
      });
      expect(e.column, 1);
      expect(e.columnName, 'zeros');
      expect(e.defect, ExogColumnDefect.allZero);
      expect(e.defects, [null, ExogColumnDefect.allZero, null]);
      expect(e, isA<ArgumentError>());
    });

    test('a constant column under d = 1 differences to zero', () {
      final e = strictError({'x1': d.x1, 'const': List.filled(d.y.length, 3)});
      expect(e.columnName, 'const');
      expect(e.defect, ExogColumnDefect.constantAfterDifferencing);
    });

    test('a linear trend under d = 2 differences to zero', () {
      final e = strictError({
        'trend': [for (var t = 0; t < d.y.length; t++) 0.1 * t + 2],
      }, order: const ArimaOrder(p: 1, d: 2, q: 1));
      expect(e.defect, ExogColumnDefect.constantAfterDifferencing);
    });

    test('a period-12 pattern under D = 1 differences to zero', () {
      final e = strictError(
        {
          'season': [for (var t = 0; t < d.y.length; t++) (t % 12) * 1.5],
        },
        order: const ArimaOrder(
          p: 0,
          d: 0,
          q: 1,
          seasonalD: 1,
          seasonalPeriod: 12,
        ),
      );
      expect(e.defect, ExogColumnDefect.constantAfterDifferencing);
    });

    test('a constant column with d = 0 duplicates the intercept', () {
      final e = strictError({
        'x1': d.x1,
        'ones': List.filled(d.y.length, 1),
      }, order: const ArimaOrder(p: 1, d: 0, q: 0));
      expect(e.defect, ExogColumnDefect.duplicatesIntercept);
    });

    test('...but is a legitimate regressor without an intercept', () {
      final f = fitSarimax(
        sarimaxRefStationaryWithMean.y,
        {'ones': List.filled(sarimaxRefStationaryWithMean.y.length, 1)},
        order: const ArimaOrder(p: 1, d: 0, q: 1),
        includeMean: false,
      );
      expect(f.regressors.single.isUsed, isTrue);
      expect(f.coefficient('ones'), closeTo(10, 1.5));
    });

    test('an exact linear combination flags the LATER column', () {
      final combo = [
        for (var t = 0; t < d.y.length; t++) 2 * d.x1[t] - d.x2[t],
      ];
      final e = strictError({'x1': d.x1, 'x2': d.x2, 'combo': combo});
      expect(e.columnName, 'combo');
      expect(e.defect, ExogColumnDefect.collinear);
    });

    test('a duplicated column is collinear', () {
      final e = strictError({'x1': d.x1, 'copy': List.of(d.x1)});
      expect(e.columnName, 'copy');
      expect(e.defect, ExogColumnDefect.collinear);
    });

    test('a nearly-but-not-exactly collinear column is accepted', () {
      final g = _Gauss(11);
      final near = [for (final v in d.x1) v + 1e-3 * g.next()];
      final f = fitSarimax(d.y, {'x1': d.x1, 'near': near}, order: _arima111);
      expect(f.regressors.every((r) => r.isUsed), isTrue);
    });
  });

  group('regressor screening — dropDegenerate policy', () {
    final d = _data();

    test('drops defective columns and fits exactly the remaining ones', () {
      final clean = fitSarimax(
        d.y,
        {'x1': d.x1, 'x2': d.x2},
        order: _arima111,
        horizon: 4,
        futureExog: {
          'x1': [1, 2, 3, 4],
          'x2': [0, 0, 1, 1],
        },
      );
      final dropped = fitSarimax(
        d.y,
        {
          'x1': d.x1,
          'zeros': List.filled(d.y.length, 0),
          'x2': d.x2,
          'copy': List.of(d.x2),
        },
        order: _arima111,
        horizon: 4,
        futureExog: {
          'x1': [1, 2, 3, 4],
          'zeros': [0, 0, 0, 0],
          'x2': [0, 0, 1, 1],
          'copy': [9, 9, 9, 9],
        },
        exogPolicy: ExogPolicy.dropDegenerate,
      );
      expect(dropped.regressors.map((r) => r.name), [
        'x1',
        'zeros',
        'x2',
        'copy',
      ]);
      expect(dropped.regressors.map((r) => r.defect), [
        null,
        ExogColumnDefect.allZero,
        null,
        ExogColumnDefect.collinear,
      ]);
      expect(dropped.usedRegressors.map((r) => r.name), ['x1', 'x2']);
      expect(dropped.coefficient('zeros'), isNull);
      expect(dropped.regressors[3].standardError, isNull);
      // The native fit saw literally the same design: identical numbers.
      expect(dropped.coefficient('x1'), clean.coefficient('x1'));
      expect(dropped.coefficient('x2'), clean.coefficient('x2'));
      expect(dropped.ar, clean.ar);
      expect(dropped.forecast, clean.forecast);
      expect(dropped.aic, clean.aic);
      expect(dropped.covariance!.parameters, ['ar.L1', 'ma.L1', 'x1', 'x2']);
    });

    test('dropping every column leaves a plain SARIMA fit', () {
      final f = fitSarimax(
        d.y,
        {'zeros': List.filled(d.y.length, 0)},
        order: _arima111,
        exogPolicy: ExogPolicy.dropDegenerate,
      );
      final s = fitSarima(d.y, order: _arima111);
      expect(f.usedRegressors, isEmpty);
      expect(f.ar, s.ar);
      expect(f.ma, s.ma);
    });
  });

  test(
    'regressors sampled at half rate are not falsely rejected (ctsa rank fix)',
    () {
      // Two full-rank regressors whose values are held for 2 steps. ctsa's
      // unpatched collinearity check read the column-major design as row-major
      // and rejected exactly this as collinear (retval 7); see PROVENANCE.md.
      final g = _Gauss(5);
      const n = 144;
      final x1 = [for (var t = 0; t < n; t++) ((t ~/ 2) % 5) + 0.1 * (t ~/ 2)];
      final x2 = [for (var t = 0; t < n; t++) math.sin((t ~/ 2) * 0.7) * 3];
      final y = <double>[];
      var u = 0.0;
      for (var t = 0; t < n; t++) {
        u = 0.5 * u + 0.3 * g.next();
        y.add(4 + 0.8 * x1[t] - 0.5 * x2[t] + u);
      }
      final f = fitSarimax(y, {
        'x1': x1,
        'x2': x2,
      }, order: const ArimaOrder(p: 1, d: 0, q: 0));
      expect(f.status, ArimaFitStatus.probableSuccess);
      expect(f.coefficient('x1'), closeTo(0.8, 0.1));
      expect(f.coefficient('x2'), closeTo(-0.5, 0.1));
    },
  );

  group('native guards', () {
    final d = _data();

    test(
      'degrees of freedom: refused before ctsa could reach exit() in tinv',
      () {
        // (1,1,1) + 2 regressors: Nused = n - 1 must exceed 4 coefficients by
        // >= 2, so n = 7 is the smallest admissible length.
        final ok = <Object>[];
        try {
          ok.add(
            fitSarimax(d.y.sublist(0, 7), {
              'x1': d.x1.sublist(0, 7),
              'x2': d.x2.sublist(0, 7),
            }, order: _arima111),
          );
        } on TseriesNumericException catch (e) {
          ok.add(e); // a bad fit is fine; the process must just survive.
        }
        expect(ok, hasLength(1));
        expect(
          () => fitSarimax(d.y.sublist(0, 6), {
            'x1': d.x1.sublist(0, 6),
            'x2': d.x2.sublist(0, 6),
          }, order: _arima111),
          throwsA(
            isA<ArgumentError>().having(
              (e) => '${e.message}',
              'message',
              contains('degrees of freedom'),
            ),
          ),
        );
      },
    );

    test('many regressors on a short series: ArgumentError, not exit()', () {
      final g = _Gauss(9);
      final exog = {
        for (var j = 0; j < 10; j++)
          'x$j': [for (var t = 0; t < 12; t++) g.next()],
      };
      expect(
        () => fitSarimax(
          d.y.sublist(0, 12),
          exog,
          order: const ArimaOrder(p: 0, d: 1, q: 0),
        ),
        throwsArgumentError,
      );
    });

    test('a seasonal MA at s = 288 is refused as too large (SARIMAX)', () {
      final g = _Gauss(1);
      final y = [for (var t = 0; t < 700; t++) g.next()];
      const order = ArimaOrder(
        p: 0,
        d: 0,
        q: 0,
        seasonalD: 1,
        seasonalQ: 1,
        seasonalPeriod: 288,
      );
      expect(
        () => fitSarimax(y, const {}, order: order),
        throwsA(
          isA<ModelTooLargeError>()
              .having(
                (e) => e.estimatedBytes,
                'estimated',
                greaterThan(1 << 32),
              )
              .having((e) => e.limitBytes, 'limit', 64 * 1024 * 1024),
        ),
      );
    });

    test('...and on the SARIMA path, through the same guard', () {
      final g = _Gauss(2);
      final y = [for (var t = 0; t < 700; t++) g.next()];
      expect(
        () => fitSarima(
          y,
          order: const ArimaOrder(
            p: 0,
            d: 0,
            q: 0,
            seasonalD: 1,
            seasonalQ: 1,
            seasonalPeriod: 288,
          ),
        ),
        throwsA(isA<ModelTooLargeError>()),
      );
    });

    test('...while D = 1 at s = 288 without P/Q stays feasible', () {
      final g = _Gauss(4);
      var level = 0.0;
      final y = [for (var t = 0; t < 700; t++) level += 0.1 * g.next()];
      final f = fitSarima(
        [for (var t = 0; t < 700; t++) y[t] + math.sin(2 * math.pi * t / 288)],
        order: const ArimaOrder(
          p: 1,
          d: 0,
          q: 1,
          seasonalD: 1,
          seasonalPeriod: 288,
        ),
        horizon: 3,
      );
      expect(f.forecast, hasLength(3));
    });

    test('arimaForecast skips an infeasible rung and falls back', () {
      final g = _Gauss(6);
      final y = [for (var t = 0; t < 700; t++) g.next()];
      final f = arimaForecast(
        y,
        orders: const [
          ArimaOrder(
            p: 0,
            d: 0,
            q: 0,
            seasonalD: 1,
            seasonalQ: 1,
            seasonalPeriod: 288,
          ),
          ArimaOrder(p: 1, d: 0, q: 0),
        ],
        horizon: 2,
      );
      expect(f.order, const ArimaOrder(p: 1, d: 0, q: 0));
      expect(f.attempts.first.accepted, isFalse);
      expect(f.attempts.first.retval, isNull);
    });

    test('non-finite regressor values are rejected with the column named', () {
      final bad = List<double>.of(d.x1)..[7] = double.nan;
      expect(
        () => fitSarimax(d.y, {'x1': bad}, order: _arima111),
        throwsA(
          isA<ArgumentError>().having(
            (e) => '${e.name}',
            'name',
            contains('exog["x1"][7]'),
          ),
        ),
      );
    });

    test('a column that overflows when differenced is rejected', () {
      final huge = [
        for (var t = 0; t < d.y.length; t++) t.isEven ? 9e149 : -9e149,
      ];
      expect(
        () => fitSarimax(d.y, {'huge': huge}, order: _arima111),
        throwsA(
          isA<ArgumentError>().having(
            (e) => '${e.message}',
            'message',
            contains('when differenced'),
          ),
        ),
      );
    });

    test('values beyond ±1e150 are refused (ctsa squares them unscaled)', () {
      final big = List<double>.of(d.x1)..[3] = 2e150;
      expect(
        () => fitSarimax(d.y, {'big': big}, order: _arima111),
        throwsA(
          isA<ArgumentError>().having(
            (e) => '${e.message}',
            'message',
            contains('"big"'),
          ),
        ),
      );
    });

    test('regressor units do not change the fit (column scaling)', () {
      // Unscaled, ctsa optimises beta in raw units and is not scale-invariant:
      // at x * 1e6 it returned status 1 at an optimum 52 log-lik units worse.
      // The shim scales each column by an exact power of two and maps beta
      // and its covariance back, so any unit gives the same model.
      final c = sarimaxRefNonSeasonal;
      final base = fitSarimax(c.y, c.exog, order: _arima111);
      final b0 = base.coefficient('x1')!;
      final se0 = base.regressors.first.standardError!;
      for (final k in [1e-100, 1e-6, 1e6, 1e100, 1e140]) {
        final f = fitSarimax(c.y, {
          'x1': [for (final v in c.exog['x1']!) v * k],
          'x2': c.exog['x2']!,
        }, order: _arima111);
        expect(f.status, ArimaFitStatus.probableSuccess, reason: '$k');
        expect(f.coefficient('x1')! * k, closeTo(b0, 1e-5 * b0.abs()));
        expect(f.regressors.first.standardError! * k, closeTo(se0, 1e-5 * se0));
        expect(f.logLikelihood, closeTo(base.logLikelihood!, 1e-6));
        expect(f.coefficient('x2'), closeTo(base.coefficient('x2')!, 1e-5));
      }
    });

    test('drift with d + D = 2 is refused', () {
      expect(
        () => fitSarimax(
          d.y,
          const {},
          order: const ArimaOrder(p: 0, d: 2, q: 1),
          includeDrift: true,
        ),
        throwsArgumentError,
      );
    });

    test('future regressors are required and must match', () {
      expect(
        () => fitSarimax(d.y, {'x1': d.x1}, order: _arima111, horizon: 3),
        throwsArgumentError,
      );
      expect(
        () => fitSarimax(
          d.y,
          {'x1': d.x1},
          order: _arima111,
          horizon: 3,
          futureExog: {
            'x1': [1, 2],
          },
        ),
        throwsArgumentError,
      );
      expect(
        () => fitSarimax(
          d.y,
          {'x1': d.x1},
          order: _arima111,
          horizon: 2,
          futureExog: {
            'x1': [1, 2],
            'other': [1, 2],
          },
        ),
        throwsArgumentError,
      );
    });

    test('reserved and mismatched regressor names are refused', () {
      for (final name in ['intercept', 'ar.L1', 'ma.S.L12', '']) {
        expect(
          () => fitSarimax(d.y, {name: d.x1}, order: _arima111),
          throwsArgumentError,
          reason: name,
        );
      }
      expect(
        () => fitSarimax(d.y, {'x1': d.x1.sublist(1)}, order: _arima111),
        throwsArgumentError,
      );
    });
  });

  group('trust gates by method', () {
    final d = _data();
    final exog = {'x1': d.x1, 'x2': d.x2};

    test('method codes match the native shim constants', () {
      expect(SarimaxMethod.cssMle.code, b.TSERIES_SARIMAX_METHOD_CSS_MLE);
      expect(SarimaxMethod.mle.code, b.TSERIES_SARIMAX_METHOD_MLE);
      expect(SarimaxMethod.css.code, b.TSERIES_SARIMAX_METHOD_CSS);
    });

    test('pure MLE agrees with CSS-then-MLE on a well-posed problem', () {
      final a = fitSarimax(d.y, exog, order: _arima111);
      final m = fitSarimax(
        d.y,
        exog,
        order: _arima111,
        method: SarimaxMethod.mle,
      );
      expect(m.status, ArimaFitStatus.probableSuccess);
      expect(m.coefficient('x1'), closeTo(a.coefficient('x1')!, 1e-3));
      expect(m.logLikelihood, closeTo(a.logLikelihood!, 1e-3));
      expect(m.aic, isNotNull);
      expect(m.covariance, isNotNull);
    });

    test('CSS: a (CSS) log-likelihood, but no AIC and no covariance', () {
      final f = fitSarimax(
        d.y,
        exog,
        order: _arima111,
        method: SarimaxMethod.css,
      );
      expect(f.status, ArimaFitStatus.probableSuccess);
      expect(f.logLikelihood, isNotNull);
      expect(f.aic, isNull);
      expect(f.covariance, isNull);
      expect(f.regressors.every((r) => r.standardError == null), isTrue);
      expect(f.coefficient('x1'), closeTo(0.8, 0.15));
    });

    test(
      'a non-converged CSS start (status 4): no AIC, log-lik or covariance',
      () {
        // Explosive AR(1) with a regressor: ctsa's CSS pre-estimate exhausts
        // its iterations and as154x returns BEFORE the MLE step (as it does on
        // 10/12), leaving CSS-only coefficients. Every diagnostic that would
        // describe an MLE fit must therefore be withheld.
        final rnd = math.Random(42);
        final y = <double>[10];
        for (var i = 1; i < 80; i++) {
          y.add(y[i - 1] * 1.15 + rnd.nextDouble() - 0.5);
        }
        final x = [for (var i = 0; i < 80; i++) math.sin(i / 3)];
        final f = fitSarimax(y, {
          'x': x,
        }, order: const ArimaOrder(p: 1, d: 0, q: 0));
        expect(f.status, isNot(ArimaFitStatus.probableSuccess));
        expect(f.status, ArimaFitStatus.maxIterations);
        expect(f.aic, isNull);
        expect(f.logLikelihood, isNull);
        expect(f.covariance, isNull);
        expect(f.regressors.single.coefficient, isNotNull);
        expect(f.regressors.single.standardError, isNull);
      },
    );
  });

  group('sarimaxForecast', () {
    final d = _data(n: 150);
    final hist = {'x1': d.x1.sublist(0, 140), 'x2': d.x2.sublist(0, 140)};
    final fut = {'x1': d.x1.sublist(140), 'x2': d.x2.sublist(140)};
    final y = d.y.sublist(0, 140);

    test('band is point ± z·se and the accepted rung is recorded', () {
      final f = sarimaxForecast(
        y,
        hist,
        orders: const [_arima111],
        horizon: 10,
        futureExog: fut,
        confidenceLevel: 0.9,
      );
      const z = 1.6448536269514722;
      for (var i = 0; i < 10; i++) {
        expect(
          f.ciLower[i],
          closeTo(f.point[i] - z * f.standardErrors[i], 1e-8),
        );
        expect(
          f.ciUpper[i],
          closeTo(f.point[i] + z * f.standardErrors[i], 1e-8),
        );
      }
      expect(f.order, _arima111);
      expect(f.attempts.single.accepted, isTrue);
      expect(f.fit.coefficient('x1'), closeTo(0.8, 0.1));
      // The regressors carry the signal: the forecast tracks y's true future.
      for (var i = 0; i < 10; i++) {
        expect(
          (f.point[i] - d.y[140 + i]).abs(),
          lessThan(4 * f.standardErrors[i] + 0.5),
        );
      }
    });

    test('rungs below the floor (counting regressors) are skipped', () {
      final f = sarimaxForecast(
        y.sublist(0, 45),
        {'x1': hist['x1']!.sublist(0, 45), 'x2': hist['x2']!.sublist(0, 45)},
        // (2,1,2): 4 ARMA + 2 regressors -> floor 61 > 45.
        orders: const [
          ArimaOrder(p: 2, d: 1, q: 2),
          ArimaOrder(p: 0, d: 1, q: 1),
        ],
        horizon: 3,
        futureExog: {
          'x1': fut['x1']!.sublist(0, 3),
          'x2': fut['x2']!.sublist(0, 3),
        },
      );
      expect(f.order, const ArimaOrder(p: 0, d: 1, q: 1));
      expect(f.attempts.first.note, contains('floor'));
    });

    test('strict regressor errors are rethrown, not laddered around', () {
      expect(
        () => sarimaxForecast(
          y,
          {...hist, 'const': List.filled(140, 2)},
          orders: const [_arima111, ArimaOrder(p: 0, d: 1, q: 1)],
          horizon: 2,
          futureExog: {
            'x1': fut['x1']!.sublist(0, 2),
            'x2': fut['x2']!.sublist(0, 2),
            'const': [2, 2],
          },
        ),
        throwsA(isA<ExogenousColumnError>()),
      );
    });

    test('dropDegenerate lets each rung drop what its differencing kills', () {
      final f = sarimaxForecast(
        y,
        {...hist, 'const': List.filled(140, 2)},
        orders: const [_arima111],
        horizon: 2,
        futureExog: {
          'x1': fut['x1']!.sublist(0, 2),
          'x2': fut['x2']!.sublist(0, 2),
          'const': [2, 2],
        },
        exogPolicy: ExogPolicy.dropDegenerate,
      );
      expect(
        f.fit.regressors.last.defect,
        ExogColumnDefect.constantAfterDifferencing,
      );
    });

    test('an infeasible seasonal rung is skipped like any other', () {
      final g = _Gauss(8);
      final yy = [for (var t = 0; t < 700; t++) g.next()];
      final xx = [for (var t = 0; t < 702; t++) math.sin(t / 5)];
      final f = sarimaxForecast(
        yy,
        {'x': xx.sublist(0, 700)},
        orders: const [
          ArimaOrder(
            p: 0,
            d: 0,
            q: 0,
            seasonalD: 1,
            seasonalQ: 1,
            seasonalPeriod: 288,
          ),
          ArimaOrder(p: 1, d: 0, q: 0),
        ],
        horizon: 2,
        futureExog: {'x': xx.sublist(700)},
      );
      expect(f.order, const ArimaOrder(p: 1, d: 0, q: 0));
      expect(f.attempts.first.note, contains('MiB'));
    });
  });
}
