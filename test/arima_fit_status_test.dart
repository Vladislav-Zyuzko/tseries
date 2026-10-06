/// Coverage for ctsa's fit-status (`retval`) surface and the trustworthiness of
/// the diagnostics derived from it.
///
/// This suite deliberately exercises the ctsa branches the golden suite CANNOT
/// reach. The golden fit is the airline model (0,1,1)(0,1,1)[12], which has
/// p == 0 and P == 0; `checkroots_cerr` only inspects the AR roots when p > 0
/// (code 10) or P > 0 (code 12), so the golden test is *structurally* incapable
/// of entering the 10/12 branch. Everything validated against R therefore
/// speaks only for the ordinary MLE path — these tests cover the rest.
library;

import 'dart:math' as math;

import 'package:test/test.dart';
import 'package:tseries/tseries.dart';

/// An explosive AR(1) (`x[t] = 1.15*x[t-1] + noise`). Fitted WITHOUT
/// differencing its AR root is far outside the unit circle, so ctsa's CSS
/// pre-estimate fails `archeck` and it returns code 10.
List<double> _explosiveSeries({int n = 80}) {
  final rnd = math.Random(42);
  final out = <double>[10];
  for (var i = 1; i < n; i++) {
    out.add(out[i - 1] * 1.15 + rnd.nextDouble() - 0.5);
  }
  return out;
}

/// An explosive SEASONAL AR (`x[t] = 1.2*x[t-12] + noise`), which drives the
/// seasonal-AR root check into code 12.
List<double> _explosiveSeasonalSeries({int n = 120, int s = 12}) {
  final rnd = math.Random(7);
  final out = <double>[for (var i = 0; i < s; i++) 10 + rnd.nextDouble()];
  for (var i = s; i < n; i++) {
    out.add(out[i - s] * 1.2 + rnd.nextDouble() - 0.5);
  }
  return out;
}

/// A well-behaved stationary AR(1) that converges normally.
List<double> _benignSeries({int n = 120}) {
  final rnd = math.Random(3);
  final out = <double>[5];
  for (var i = 1; i < n; i++) {
    out.add(5 + 0.6 * (out[i - 1] - 5) + (rnd.nextDouble() - 0.5) * 0.4);
  }
  return out;
}

void main() {
  group('ArimaFitStatus', () {
    test('maps every documented ctsa code', () {
      expect(ArimaFitStatus.fromNative(0), ArimaFitStatus.notRun);
      expect(ArimaFitStatus.fromNative(1), ArimaFitStatus.probableSuccess);
      expect(ArimaFitStatus.fromNative(4), ArimaFitStatus.maxIterations);
      expect(ArimaFitStatus.fromNative(7), ArimaFitStatus.collinearExogenous);
      expect(ArimaFitStatus.fromNative(10), ArimaFitStatus.nonStationaryAr);
      expect(
        ArimaFitStatus.fromNative(12),
        ArimaFitStatus.nonStationarySeasonalAr,
      );
      expect(ArimaFitStatus.fromNative(15), ArimaFitStatus.nonFinite);
    });

    test('maps an undocumented code to unknown rather than throwing', () {
      // Never trust the native side to stay inside a known set: an unexpected
      // int must degrade to `unknown`, not crash the caller.
      expect(ArimaFitStatus.fromNative(99), ArimaFitStatus.unknown);
      expect(ArimaFitStatus.fromNative(-3), ArimaFitStatus.unknown);
    });
  });

  group('retval is reported, never enforced', () {
    test('a benign series fits with probableSuccess', () {
      final fit = fitSarima(
        _benignSeries(),
        order: const ArimaOrder(p: 1, d: 0, q: 0),
        horizon: 3,
      );
      expect(fit.retval, 1);
      expect(fit.status, ArimaFitStatus.probableSuccess);
    });

    test('code 10 (non-stationary AR) is reachable and still returns a fit', () {
      // The point of Step A: we surface the code without acting on it. If this
      // ever starts throwing, the "report, do not reject" contract has been
      // broken silently.
      final fit = fitSarima(
        _explosiveSeries(),
        order: const ArimaOrder(p: 1, d: 0, q: 0),
        horizon: 3,
      );
      expect(fit.retval, 10);
      expect(fit.status, ArimaFitStatus.nonStationaryAr);
      expect(fit.forecast, hasLength(3));
    });

    test('code 12 (non-stationary seasonal AR) is reachable', () {
      final fit = fitSarima(
        _explosiveSeasonalSeries(),
        order: const ArimaOrder(
          p: 0,
          d: 0,
          q: 0,
          seasonalP: 1,
          seasonalPeriod: 12,
        ),
        horizon: 3,
      );
      expect(fit.retval, 12);
      expect(fit.status, ArimaFitStatus.nonStationarySeasonalAr);
    });
  });

  group('aic/logLikelihood are null unless genuinely computed', () {
    test('MLE + probableSuccess reports both', () {
      final fit = fitSarima(
        _benignSeries(),
        order: const ArimaOrder(p: 1, d: 0, q: 0),
        horizon: 1,
      );
      expect(fit.retval, 1);
      expect(fit.aic, isNotNull);
      expect(fit.logLikelihood, isNotNull);
    });

    test(
      'CSS reports a log-likelihood but NO aic (ctsa never computes one)',
      () {
        final fit = fitSarima(
          _benignSeries(),
          order: const ArimaOrder(p: 1, d: 0, q: 0),
          horizon: 1,
          method: ArimaMethod.css,
        );
        expect(fit.logLikelihood, isNotNull);
        // Before this was gated, the field held malloc garbage from sarima_init.
        expect(fit.aic, isNull);
      },
    );

    test('Box-Jenkins reports neither (ctsa assigns neither field)', () {
      final fit = fitSarima(
        _benignSeries(),
        order: const ArimaOrder(p: 1, d: 0, q: 0),
        horizon: 1,
        method: ArimaMethod.boxJenkins,
      );
      expect(fit.logLikelihood, isNull);
      expect(fit.aic, isNull);
    });

    test('code 10 reports NO aic: it would be a CSS-based, incomparable AIC', () {
      // On 10/12 ctsa returns before the MLE step, so `loglik` still holds a
      // CSS log-likelihood while sarima_exec computes `aic` from it anyway.
      // Reporting that as an MLE AIC would silently corrupt any AIC comparison.
      final fit = fitSarima(
        _explosiveSeries(),
        order: const ArimaOrder(p: 1, d: 0, q: 0),
        horizon: 1,
      );
      expect(fit.retval, 10);
      expect(fit.aic, isNull);
      expect(fit.logLikelihood, isNull);
    });
  });

  group('degenerate fits are rejected on a field we trust', () {
    test('a constant series has sigma2 == 0 and is refused', () {
      // Such a model fits the input exactly and reports zero forecast
      // uncertainty — "confident nonsense". It must not be handed back.
      expect(
        () => fitSarima(
          List<double>.filled(60, 5.5),
          order: const ArimaOrder(p: 2, d: 1, q: 1),
          horizon: 3,
        ),
        throwsA(isA<TseriesNumericException>()),
      );
    });

    test('the ladder surfaces the refusal rather than a zero-width band', () {
      expect(
        () => arimaForecast(
          List<double>.filled(60, 5.5),
          orders: const [
            ArimaOrder(p: 2, d: 1, q: 1),
            ArimaOrder(p: 0, d: 1, q: 1),
          ],
          horizon: 3,
        ),
        throwsA(isA<TseriesNumericException>()),
      );
    });
  });

  group('ladder attempts are measurable', () {
    test('every rung is recorded, in order, with the accepted one last', () {
      final forecast = arimaForecast(
        _benignSeries(),
        orders: const [
          ArimaOrder(p: 2, d: 1, q: 1),
          ArimaOrder(p: 1, d: 1, q: 1),
          ArimaOrder(p: 0, d: 1, q: 1),
        ],
        horizon: 3,
      );
      expect(forecast.attempts, isNotEmpty);
      expect(forecast.attempts.last.accepted, isTrue);
      expect(forecast.attempts.last.order, forecast.order);
      expect(forecast.attempts.last.retval, forecast.retval);
      // Only the final rung is ever the accepted one.
      expect(forecast.attempts.where((a) => a.accepted), hasLength(1));
    });

    test('a rung skipped below the length floor has no ctsa code', () {
      // 20 points is below the floor of (2,1,1)=31 but above (0,1,1)=11, so the
      // ladder must skip the first rung without ever entering C.
      final short = _benignSeries(n: 20);
      final forecast = arimaForecast(
        short,
        orders: const [
          ArimaOrder(p: 2, d: 1, q: 1),
          ArimaOrder(p: 0, d: 1, q: 1),
        ],
        horizon: 2,
      );
      expect(forecast.attempts, hasLength(2));
      final skipped = forecast.attempts.first;
      expect(skipped.accepted, isFalse);
      expect(skipped.retval, isNull, reason: 'native never ran for this rung');
      expect(skipped.status, isNull);
      expect(skipped.note, contains('floor'));
      expect(forecast.attempts.last.accepted, isTrue);
      expect(forecast.order, const ArimaOrder(p: 0, d: 1, q: 1));
    });

    test(
      'attempts can be aggregated by order (ArimaOrder has value equality)',
      () {
        const a = ArimaOrder(p: 1, d: 1, q: 1);
        const b = ArimaOrder(p: 1, d: 1, q: 1);
        expect(a, b);
        expect(a.hashCode, b.hashCode);
        expect(a, isNot(const ArimaOrder(p: 2, d: 1, q: 1)));

        // Exercise the actual downstream use: the feature layer keys a
        // histogram on the order, which only collapses correctly if equal
        // orders hash and compare equal.
        const other = ArimaOrder(p: 2, d: 1, q: 1);
        final histogram = <ArimaOrder, int>{};
        for (final order in <ArimaOrder>[a, b, other]) {
          histogram[order] = (histogram[order] ?? 0) + 1;
        }
        expect(histogram, hasLength(2));
        expect(histogram[a], 2);
        expect(histogram[other], 1);
      },
    );

    test('(0,1,1) has p == 0 so the ladder is always reachable to its end', () {
      // checkroots_cerr only inspects AR roots when p > 0 / P > 0, so this rung
      // structurally cannot return 10 or 12. This is what makes it a sound
      // last resort for a degradation ladder.
      final fit = fitSarima(
        _explosiveSeries(),
        order: const ArimaOrder(p: 0, d: 1, q: 1),
        horizon: 3,
      );
      expect(fit.retval, isNot(10));
      expect(fit.retval, isNot(12));
    });
  });
}
