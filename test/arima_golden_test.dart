/// Golden-value validation of the native ctsa core against published
/// reference fits. This is the whole point of the POC: proving ctsa is
/// numerically correct, not merely that it "compiles and returns a number".
///
/// Reference: the classic multiplicative "airline model"
/// SARIMA(0,1,1)(0,1,1)[12] fitted to `log(AirPassengers)` (Box & Jenkins
/// Series G). The R `stats::arima` / `forecast::Arima` fit is one of the most
/// widely published results in time-series analysis:
///
///   ma1 (θ)   = -0.4018
///   sma1 (Θ)  = -0.5569
///   sigma^2   =  0.001348
///   logLik    =  244.7
///   aic       = -483.4
///
/// Sources: Hyndman & Athanasopoulos, "Forecasting: Principles and Practice"
/// (otexts.com/fpp2/seasonal-arima.html); R `stats`/`forecast` documentation.
///
/// SIGN CONVENTION: ctsa parameterises the MA polynomial as θ(B) = 1 − θ₁B,
/// whereas R/statsmodels use 1 + θ₁B. ctsa's MA estimates therefore come out
/// with the opposite sign to R's; we compare magnitudes and document the flip.
library;

import 'dart:math' as math;

import 'package:test/test.dart';
import 'package:tseries/tseries.dart';

import 'fixtures/air_passengers.dart';

void main() {
  group('airline model — SARIMA(0,1,1)(0,1,1)[12] on log(AirPassengers)', () {
    late final SarimaFitResult fit;

    setUpAll(() {
      final logAir = airPassengers.map((v) => math.log(v)).toList();
      fit = fitSarima(
        logAir,
        order: const ArimaOrder(
          p: 0,
          d: 1,
          q: 1,
          seasonalP: 0,
          seasonalD: 1,
          seasonalQ: 1,
          seasonalPeriod: 12,
        ),
        horizon: 12,
      );
    });

    test('MA coefficient matches R |ma1| = 0.4018', () {
      // ctsa ~ 0.3991 vs R 0.4018 -> gap 0.0027. Tolerance 0.01 is tight.
      expect(fit.ma, hasLength(1));
      expect(fit.ma.single.abs(), closeTo(0.4018, 0.01));
      // ctsa's sign convention is opposite R's (1 - θB), so it reports it > 0.
      expect(fit.ma.single, greaterThan(0));
    });

    test('seasonal MA coefficient matches R |sma1| = 0.5569', () {
      // ctsa ~ 0.5544 vs R 0.5569 -> gap 0.0025.
      expect(fit.seasonalMa, hasLength(1));
      expect(fit.seasonalMa.single.abs(), closeTo(0.5569, 0.01));
    });

    test('residual variance matches R sigma^2 = 0.001348', () {
      // ctsa ~ 0.0013479 -> agrees to ~4 significant figures.
      expect(fit.sigma2, closeTo(0.001348, 5e-5));
    });

    test('log-likelihood matches R logLik = 244.7', () {
      expect(fit.logLikelihood, closeTo(244.7, 0.5));
    });

    test('AIC matches R aic = -483.4', () {
      expect(fit.aic, closeTo(-483.4, 1.0));
    });

    test('12-step forecast is finite and in a plausible level range', () {
      final level = fit.forecast.map(math.exp).toList();
      expect(level, hasLength(12));
      for (final v in level) {
        expect(v.isFinite, isTrue);
        // 1961 passenger levels continue the 1960 range (~ 350–700k).
        expect(v, inInclusiveRange(340, 720));
      }
      // The model must reproduce the summer peak / winter trough seasonality.
      expect(level.reduce(math.max), greaterThan(600));
      expect(level.reduce(math.min), lessThan(480));
    });
  });

  group('input guards (safe wrapper rejects bad input, no native crash)', () {
    test('non-finite input throws ArgumentError', () {
      expect(
        () => fitSarima([
          1,
          2,
          double.nan,
          4,
          5,
        ], order: const ArimaOrder(p: 1, d: 0, q: 0)),
        throwsArgumentError,
      );
      expect(
        () => fitSarima([
          1,
          2,
          double.infinity,
          4,
          5,
        ], order: const ArimaOrder(p: 1, d: 0, q: 0)),
        throwsArgumentError,
      );
    });

    test('too-short series throws ArgumentError', () {
      expect(
        () => fitSarima([1, 2], order: const ArimaOrder(p: 0, d: 0, q: 0)),
        throwsArgumentError,
      );
    });
  });
}
