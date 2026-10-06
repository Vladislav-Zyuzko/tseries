/// Contract tests for the robust forecast API [arimaForecast]. These do NOT
/// re-validate the native numeric core (that is done in
/// `arima_golden_test.dart`); they validate the wrapper contract: that the
/// forecast-oriented method does not distort the numbers coming out of
/// [fitSarima], builds the confidence band correctly, degrades along the order
/// ladder on a hard/short series, and rejects bad input cleanly.
library;

import 'dart:math' as math;

import 'package:test/test.dart';
import 'package:tseries/tseries.dart';

import 'fixtures/air_passengers.dart';

void main() {
  final logAir = airPassengers.map((v) => math.log(v)).toList();
  const airline = ArimaOrder(
    p: 0,
    d: 1,
    q: 1,
    seasonalP: 0,
    seasonalD: 1,
    seasonalQ: 1,
    seasonalPeriod: 12,
  );

  group('arimaForecast — numeric passthrough (does not distort fitSarima)', () {
    test('point and stdErr match a direct fitSarima of the same order', () {
      const horizon = 12;
      final direct = fitSarima(logAir, order: airline, horizon: horizon);
      final forecast = arimaForecast(
        logAir,
        orders: const [airline],
        horizon: horizon,
      );

      // The wrapper must hand back exactly the native numbers, untouched.
      expect(forecast.point, orderedEquals(direct.forecast));
      expect(forecast.standardErrors, orderedEquals(direct.standardErrors));
      // Diagnostics/coefficients are preserved too.
      expect(forecast.sigma2, direct.sigma2);
      expect(forecast.aic, direct.aic);
      expect(forecast.logLikelihood, direct.logLikelihood);
      expect(forecast.ma, orderedEquals(direct.ma));
      expect(forecast.seasonalMa, orderedEquals(direct.seasonalMa));
      // The order that fitted is the (only) candidate supplied.
      expect(forecast.order.toString(), airline.toString());
    });
  });

  group('arimaForecast — shape and confidence band', () {
    late final ArimaForecast forecast;
    const horizon = 9;

    setUpAll(() {
      forecast = arimaForecast(
        logAir,
        orders: const [airline],
        horizon: horizon,
      );
    });

    test('every output list has length == horizon', () {
      expect(forecast.point, hasLength(horizon));
      expect(forecast.standardErrors, hasLength(horizon));
      expect(forecast.ciLower, hasLength(horizon));
      expect(forecast.ciUpper, hasLength(horizon));
    });

    test('CI == point ± 1.96·stdErr at the default 95% level', () {
      expect(forecast.confidenceLevel, 0.95);
      for (var i = 0; i < horizon; i++) {
        final halfWidth = 1.96 * forecast.standardErrors[i];
        expect(
          forecast.ciUpper[i] - forecast.point[i],
          closeTo(halfWidth, 1e-4),
        );
        expect(
          forecast.point[i] - forecast.ciLower[i],
          closeTo(halfWidth, 1e-4),
        );
        // Band is symmetric about the point forecast to full precision.
        expect(
          forecast.ciUpper[i] - forecast.point[i],
          closeTo(forecast.point[i] - forecast.ciLower[i], 1e-12),
        );
        // Non-degenerate band (positive standard errors on real data).
        expect(forecast.ciUpper[i], greaterThan(forecast.ciLower[i]));
      }
    });

    test('a wider confidence level yields a wider band', () {
      final wide = arimaForecast(
        logAir,
        orders: const [airline],
        horizon: horizon,
        confidenceLevel: 0.99,
      );
      // z(0.99) ≈ 2.576 > z(0.95) ≈ 1.96, so each half-width must grow.
      for (var i = 0; i < horizon; i++) {
        final w95 = forecast.ciUpper[i] - forecast.point[i];
        final w99 = wide.ciUpper[i] - wide.point[i];
        expect(w99, greaterThan(w95));
      }
    });
  });

  group('arimaForecast — degradation ladder', () {
    test(
      'falls back to a simpler order when the first is too short to fit',
      () {
        // 25 points: below the numeric floor of (2,1,1) (=31) but above that of
        // (1,1,1) (=21) and (0,1,1) (=11). The ladder must skip (2,1,1) and use
        // the next order that fits, reporting which one actually worked.
        final shortSeries = airPassengers.sublist(0, 25);
        final result = arimaForecast(
          shortSeries,
          orders: const [
            ArimaOrder(p: 2, d: 1, q: 1),
            ArimaOrder(p: 1, d: 1, q: 1),
            ArimaOrder(p: 0, d: 1, q: 1),
          ],
          horizon: 4,
        );

        // It degraded away from the first (richest) candidate.
        expect(result.order.p, lessThan(2));
        expect(result.point, hasLength(4));
        for (final v in result.point) {
          expect(v.isFinite, isTrue);
        }
        for (final v in result.standardErrors) {
          expect(v.isFinite, isTrue);
        }
      },
    );

    test('a single-order ladder behaves like a plain fixed-order forecast', () {
      final result = arimaForecast(
        logAir,
        orders: const [ArimaOrder(p: 1, d: 1, q: 1)],
        horizon: 3,
      );
      expect(result.order.toString(), 'ARIMA(1,1,1)');
      expect(result.point, hasLength(3));
    });
  });

  group('arimaForecast — input guards', () {
    test('series shorter than the numeric floor of every order throws', () {
      // Smallest candidate is (0,1,1) with floor 11; 8 points is below it.
      expect(
        () => arimaForecast(
          List<double>.filled(8, 1),
          orders: const [
            ArimaOrder(p: 2, d: 1, q: 1),
            ArimaOrder(p: 0, d: 1, q: 1),
          ],
          horizon: 2,
        ),
        throwsArgumentError,
      );
    });

    test('non-finite input (NaN/Inf) throws ArgumentError', () {
      final withNan = List<double>.filled(40, 1)..[10] = double.nan;
      final withInf = List<double>.filled(40, 1)..[10] = double.infinity;
      expect(
        () => arimaForecast(
          withNan,
          orders: const [ArimaOrder(p: 1, d: 1, q: 1)],
          horizon: 2,
        ),
        throwsArgumentError,
      );
      expect(
        () => arimaForecast(
          withInf,
          orders: const [ArimaOrder(p: 1, d: 1, q: 1)],
          horizon: 2,
        ),
        throwsArgumentError,
      );
    });

    test('empty order ladder throws ArgumentError', () {
      expect(
        () => arimaForecast(logAir, orders: const [], horizon: 4),
        throwsArgumentError,
      );
    });

    test('non-positive horizon throws ArgumentError', () {
      expect(
        () => arimaForecast(logAir, orders: const [airline], horizon: 0),
        throwsArgumentError,
      );
    });

    test('confidence level outside (0, 1) throws ArgumentError', () {
      for (final bad in [0.0, 1.0, -0.1, 1.5]) {
        expect(
          () => arimaForecast(
            logAir,
            orders: const [airline],
            horizon: 4,
            confidenceLevel: bad,
          ),
          throwsArgumentError,
        );
      }
    });

    test('a malformed order (seasonal term without a period) throws', () {
      expect(
        () => arimaForecast(
          logAir,
          orders: const [ArimaOrder(p: 0, d: 1, q: 1, seasonalQ: 1)],
          horizon: 4,
        ),
        throwsArgumentError,
      );
    });
  });
}
