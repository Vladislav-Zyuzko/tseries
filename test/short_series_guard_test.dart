/// Structural length guards of `fitSarima` / `arimaForecast` (PROVENANCE.md,
/// "Local modifications" 15): inputs on which ctsa reads uninitialised memory —
/// and so returns run-to-run different results — are refused with a typed
/// [SeriesTooShortError] before the native call, and the forecast ladder skips
/// such an order with a note instead of failing.
library;

import 'dart:math' as math;

import 'package:test/test.dart';
import 'package:tseries/tseries.dart';

List<double> _series(int n) => [
  for (var i = 0; i < n; i++)
    5 + 0.1 * i + math.sin(i * 0.9) + (i % 12) * 0.3 + math.cos(i * 2.1) * 0.5,
];

Matcher _tooShort({required int length, required int minimum}) => throwsA(
  isA<SeriesTooShortError>()
      .having((e) => e.length, 'length', length)
      .having((e) => e.minimumLength, 'minimumLength', minimum),
);

/// Runs [fit] and fails only on a [SeriesTooShortError]: at the boundary the
/// fit is allowed to be numerically poor, just not refused for its length.
void _notRefused(void Function() fit) {
  try {
    fit();
  } on SeriesTooShortError catch (e) {
    fail('unexpectedly refused: $e');
  } on TseriesNumericException {
    // A legitimate outcome for a tiny series; not a length refusal.
  }
}

void main() {
  group('MLE/CSS: the differenced series must be longer than p + s·P', () {
    // p + s·P = 13: ctsa's CSS step conditions on 13 observations.
    const order = ArimaOrder(
      p: 1,
      d: 0,
      q: 0,
      seasonalP: 1,
      seasonalPeriod: 12,
    );
    for (final m in [ArimaMethod.mle, ArimaMethod.css]) {
      test('${m.name}: refused at Nd = p + s·P, accepted one later', () {
        expect(
          () => fitSarima(_series(13), order: order, method: m),
          _tooShort(length: 13, minimum: 14),
        );
        _notRefused(
          () => fitSarima(_series(14), order: order, method: m, horizon: 2),
        );
      });
    }

    test('differencing is counted (Nd = n − d − s·D)', () {
      const o = ArimaOrder(
        p: 0,
        d: 1,
        q: 0,
        seasonalP: 1,
        seasonalD: 1,
        seasonalPeriod: 4,
      );
      // d + s·D = 5 observations consumed, p + s·P = 4 -> n ≥ 5 + 4 + 1 = 10.
      expect(
        () => fitSarima(_series(9), order: o),
        _tooShort(length: 9, minimum: 10),
      );
      _notRefused(() => fitSarima(_series(10), order: o));
    });

    test('seasonal MA terms alone do not raise the minimum', () {
      // Nd = 8 ≤ q + s·Q = 12 is fine: the CSS sum is not empty, no
      // uninitialised read, identical results run to run.
      const o = ArimaOrder(p: 0, d: 0, q: 0, seasonalQ: 1, seasonalPeriod: 12);
      _notRefused(() => fitSarima(_series(8), order: o));
      final a = fitSarima(_series(8), order: o, horizon: 3);
      final b = fitSarima(_series(8), order: o, horizon: 3);
      expect(b.seasonalMa, a.seasonalMa);
      expect(b.forecast, a.forecast);
    });
  });

  group('Box-Jenkins: seasonal fits need (P + Q + 1)·s + 2 observations', () {
    // ctsa's seasonal start values read autocovariances up to lag (P + Q)·s,
    // computed only up to lag Nd − 2: uninitialised for Nd ≤ (P + Q)·s + 1.
    // The guard adds one season: Nd ≥ (P + Q + 1)·s + 2.
    const o = ArimaOrder(p: 0, d: 0, q: 0, seasonalP: 1, seasonalPeriod: 4);
    test('refused below, accepted at the minimum', () {
      expect(
        () => fitSarima(_series(9), order: o, method: ArimaMethod.boxJenkins),
        _tooShort(length: 9, minimum: 10),
      );
      _notRefused(
        () => fitSarima(_series(10), order: o, method: ArimaMethod.boxJenkins),
      );
    });

    test('non-seasonal Box-Jenkins is unaffected', () {
      _notRefused(
        () => fitSarima(
          _series(8),
          order: const ArimaOrder(p: 1, d: 0, q: 1),
          method: ArimaMethod.boxJenkins,
        ),
      );
    });
  });

  group('arimaForecast', () {
    test('skips an order below its structural minimum, with a note', () {
      // (0,0,0)(1,0,0)[24]: numeric floor 20, structural minimum 25. At 22
      // points the first rung passes the numeric floor but not the
      // structural minimum; the ladder moves on to (0,1,1).
      final result = arimaForecast(
        _series(22),
        orders: const [
          ArimaOrder(p: 0, d: 0, q: 0, seasonalP: 1, seasonalPeriod: 24),
          ArimaOrder(p: 0, d: 1, q: 1),
        ],
        horizon: 3,
      );
      expect(result.order, const ArimaOrder(p: 0, d: 1, q: 1));
      final first = result.attempts.first;
      expect(first.accepted, isFalse);
      expect(first.retval, isNull, reason: 'refused before reaching C');
      expect(first.note, contains('structural minimum 25'));
    });

    test('nothing long enough: ArgumentError, not a ladder exhaustion', () {
      expect(
        () => arimaForecast(
          _series(22),
          orders: const [
            ArimaOrder(p: 0, d: 0, q: 0, seasonalP: 1, seasonalPeriod: 24),
          ],
          horizon: 3,
        ),
        throwsArgumentError,
      );
    });

    test('Sweet Limit ladder is never refused at its working lengths', () {
      // The app's backtest scores an origin only after 60 in-segment steps,
      // so the ARIMA window is 61..144 points; the ladder's structural
      // minima (4, 3, 2) are far below its numeric floors (31, 21, 11).
      const ladder = [
        ArimaOrder(p: 2, d: 1, q: 1),
        ArimaOrder(p: 1, d: 1, q: 1),
        ArimaOrder(p: 0, d: 1, q: 1),
      ];
      for (final n in [11, 21, 31, 61, 62, 90, 144]) {
        try {
          final r = arimaForecast(_series(n), orders: ladder, horizon: 6);
          for (final a in r.attempts) {
            expect(a.note ?? '', isNot(contains('structural')), reason: '$n');
          }
        } on TseriesNumericException catch (e) {
          expect('$e', isNot(contains('structural')), reason: '$n');
        }
      }
    });
  });

  group('fitSarimax: regressors under d > 0 and D > 0', () {
    // ctsa's regression set-up writes past a heap buffer once d × (regressors
    // used) exceeds the differenced length (PROVENANCE.md, Local
    // modifications 16): the shim refuses it.
    const order = ArimaOrder(
      p: 0,
      d: 2,
      q: 0,
      seasonalD: 1,
      seasonalPeriod: 4,
    );
    Map<String, List<double>> exog(int n) => {
      for (var j = 0; j < 3; j++)
        'x$j': [
          for (var i = 0; i < n; i++)
            math.sin((i + 1) * (0.7 + j)) * 3 + 0.01 * j * i * i,
        ],
    };

    test('d·r = 6 > Nd = 5 (n = 11) is refused', () {
      for (final method in SarimaxMethod.values) {
        expect(
          () => fitSarimax(_series(11), exog(11), order: order, method: method),
          throwsA(
            isA<ArgumentError>().having(
              (e) => e.message,
              'message',
              contains('must not exceed the differenced length 5'),
            ),
          ),
        );
      }
    });

    test('d·r = 6 = Nd (n = 12) is fitted', () {
      for (final method in SarimaxMethod.values) {
        try {
          fitSarimax(_series(12), exog(12), order: order, method: method);
        } on TseriesNumericException {
          // A legitimate outcome on 6 differenced points; not a refusal.
        }
      }
    });
  });
}
