/// autoArima: API contract, journal invariants, and gate F (robustness) of
/// the specification (§9). Fast; part of the default `dart test` run.
library;

import 'dart:math' as math;
import 'dart:typed_data';

import 'package:test/test.dart';
import 'package:tseries/src/domain/internal/auto_arima_candidate.dart';
import 'package:tseries/tseries.dart';

import 'support/arima_simulator.dart';
import 'support/auto_arima_fixtures.dart';

/// Journal invariants that must hold for every search.
void _checkJournal(
  List<CandidateRecord> trace,
  int fits, {
  required int n,
  required DifferencingReport differencing,
}) {
  var calls = 0;
  for (var i = 0; i < trace.length; i++) {
    final r = trace[i];
    expect(r.index, i);
    // No length-rule defect may ever surface (spec §5, ред. 3).
    expect(r.detail ?? '', isNot(contains('DEFECT')), reason: '$r');
    expect(r.verdict, isNot(CandidateVerdict.rejectedFitError), reason: '$r');
    expect(r.order.d, differencing.d);
    expect(r.order.seasonalD, differencing.seasonalD);
    expect(
      r.effectiveLength,
      n -
          differencing.d -
          differencing.seasonalD * math.max(differencing.period, 1),
    );
    if (r.verdict == CandidateVerdict.cachedDuplicate) {
      final orig = trace[r.duplicateOf!];
      expect(orig.index, lessThan(i));
      expect(orig.order, r.order);
      expect(orig.constant, r.constant);
      expect(orig.verdict, isNot(CandidateVerdict.cachedDuplicate));
      continue;
    }
    expect(r.duplicateOf, isNull);
    if (r.verdict != CandidateVerdict.rejectedTooShort) calls++;
    if (r.verdict == CandidateVerdict.accepted) {
      expect(r.criterion, isNotNull);
      expect(r.ctsaStatus, 1);
    }
  }
  expect(calls, fits);
}

void main() {
  final nile = loadReferenceSeries('nile');
  final air = loadReferenceSeries('air_passengers_log');

  group('result contract', () {
    test('Nile stepwise: order, criterion, fit identity, journal', () {
      final r = autoArima(nile, horizon: 5);
      expect(r.order, const ArimaOrder(p: 1, d: 1, q: 1));
      expect(r.constant, isFalse);
      expect(r.search, ArimaSearch.stepwise);
      expect(r.criterion, InformationCriterion.aicc);
      expect(r.truncated, isFalse);
      _checkJournal(
        r.trace,
        r.fitsEvaluated,
        n: 100,
        differencing: r.differencing,
      );

      // The selected record is the last one that became current, and the
      // result's fit is that record's fit (no refit).
      final current = r.trace.lastWhere((t) => t.becameCurrent);
      final selected = current.duplicateOf == null
          ? current
          : r.trace[current.duplicateOf!];
      expect(identical(selected.fit, r.fit), isTrue);
      expect(r.criterionValue, selected.criterion);

      // Criterion from logLik, K, T (spec §4).
      expect(
        r.criterionValue,
        InformationCriterion.aicc.evaluate(
          logLikelihood: r.fit.logLikelihood!,
          parameters: 3,
          observations: 99,
        ),
      );

      // Forecast and band.
      expect(r.forecast, hasLength(5));
      final z = 1.959963984540054;
      for (var i = 0; i < 5; i++) {
        expect(r.lower[i], lessThan(r.forecast[i]));
        expect(
          r.upper[i] - r.forecast[i],
          closeTo(z * r.fit.standardErrors[i], 1e-6),
        );
      }

      // The same fit through the public fixed-order path is bit-identical.
      final direct = fitSarimax(
        nile,
        const {},
        order: r.order,
        horizon: 5,
        includeMean: false,
      );
      expect(direct.logLikelihood, r.fit.logLikelihood);
      expect(direct.ar, r.fit.ar);
      expect(direct.forecast, r.forecast);
    });

    test('differencing report (Nile): KPSS steps and sources', () {
      final dr = autoArima(nile).differencing;
      expect(dr.d, 1);
      expect(dr.dSource, DifferencingSource.test);
      expect(dr.kpss, hasLength(2));
      expect(dr.kpss[0].rejectsStationarity, isTrue);
      expect(dr.kpss[1].rejectsStationarity, isFalse);
      expect(dr.kpss[1].length, 99);
      expect(dr.seasonalD, 0);
      expect(dr.seasonalDSource, DifferencingSource.notSeasonal);
      expect(dr.seasonalStrength, isNull);
    });

    test('user-given d and D: no tests run', () {
      final r = autoArima(air, period: 12, d: 1, seasonalD: 1);
      expect(r.differencing.dSource, DifferencingSource.user);
      expect(r.differencing.seasonalDSource, DifferencingSource.user);
      expect(r.differencing.kpss, isEmpty);
      expect(r.differencing.seasonalStrength, isNull);
      // d + D = 2: no constant anywhere in the journal.
      expect(r.trace.every((t) => !t.constant), isTrue);
      expect(
        r.order,
        const ArimaOrder(
          p: 0,
          d: 1,
          q: 1,
          seasonalP: 0,
          seasonalD: 1,
          seasonalQ: 1,
          seasonalPeriod: 12,
        ),
      );
    });

    test('seasonal D: stepped rule by full cycles (spec §15.1)', () {
      expect(autoArimaSeasonalStrengthThreshold, 0.35);
      expect(autoArimaSeasonalStrengthThresholdFewCycles, 0.70);
      expect(autoArimaSeasonalMinCycles, 5);
      expect(autoArimaSeasonalManyCycles, 10);
      // log AirPassengers: k = 12 → τ = 0.35, F_S ≈ 0.933.
      final a = autoArima(air, period: 12, maxP: 1, maxQ: 1).differencing;
      expect(a.seasonalCycles, 12);
      expect(a.seasonalThreshold, 0.35);
      expect(a.seasonalD, 1);
      expect(a.seasonalDSource, DifferencingSource.test);
      expect(a.seasonalStrength, closeTo(0.9333231826436825, 1e-12));
      // USAccDeaths: k = 6 → τ_hi = 0.70, F_S ≈ 0.936.
      final usAcc = loadReferenceSeries('us_acc_deaths');
      final u = autoArima(usAcc, period: 12, maxP: 1, maxQ: 1).differencing;
      expect(u.seasonalCycles, 6);
      expect(u.seasonalThreshold, 0.70);
      expect(u.seasonalD, 1);
      // 4 cycles (59 points): no automatic D, F_S not computed.
      final w = autoArima(
        usAcc.sublist(0, 59),
        period: 12,
        maxP: 1,
        maxQ: 1,
      ).differencing;
      expect(w.seasonalCycles, 4);
      expect(w.seasonalD, 0);
      expect(w.seasonalDSource, DifferencingSource.skippedShort);
      expect(w.seasonalStrength, isNull);
      expect(w.seasonalThreshold, isNull);
      // Known seasonality: the caller passes seasonalD: 1.
      final f = autoArima(
        usAcc.sublist(0, 59),
        period: 12,
        seasonalD: 1,
        maxP: 1,
        maxQ: 1,
      ).differencing;
      expect(f.seasonalD, 1);
      expect(f.seasonalDSource, DifferencingSource.user);
      expect(f.seasonalCycles, 4);
      // Non-seasonal: no cycles.
      expect(autoArima(nile, maxP: 1).differencing.seasonalCycles, isNull);
    });

    test('seasonal D at 5–9 cycles needs F_S ≥ 0.70', () {
      // A weak deterministic pattern over noise, 6 cycles of s = 12.
      final w = simulateArima(n: 72, seed: 5, ar: [0.3]);
      final y = [
        for (var t = 0; t < 72; t++)
          10 + w[t] + 0.8 * math.sin(2 * math.pi * t / 12),
      ];
      final dr = autoArima(y, period: 12, maxP: 1, maxQ: 1).differencing;
      final fs = dr.seasonalStrength!;
      expect(dr.seasonalThreshold, 0.70);
      expect(dr.seasonalD, fs >= 0.70 ? 1 : 0);
    });

    test('allowConstant: false', () {
      final r = autoArima(loadReferenceSeries('lynx'), allowConstant: false);
      expect(r.trace.every((t) => !t.constant), isTrue);
    });

    test('maxD: the test at maxD is run and reported (spec §2.1, ред. 5)', () {
      // log AirPassengers level: η ≫ crit, so maxD = 0 is limitedByMax.
      final r = autoArima(air, period: 12, maxD: 0, maxSeasonalD: 0, maxP: 1);
      expect(r.differencing.d, 0);
      expect(r.differencing.kpss, hasLength(1));
      expect(r.differencing.kpss.single.rejectsStationarity, isTrue);
      expect(r.differencing.dSource, DifferencingSource.limitedByMax);
      expect(r.differencing.seasonalDSource, DifferencingSource.limitedByMax);
      // lynx: stationary, so the test at maxD = 0 decides: source test.
      final l = autoArima(loadReferenceSeries('lynx'), maxD: 0, maxP: 1);
      expect(l.differencing.dSource, DifferencingSource.test);
      expect(l.differencing.kpss, hasLength(1));
      // Nile with maxD = 1: η at d = 1 is computed and does not reject.
      final n1 = autoArima(nile, maxD: 1, maxP: 1);
      expect(n1.differencing.d, 1);
      expect(n1.differencing.kpss, hasLength(2));
      expect(n1.differencing.dSource, DifferencingSource.test);
    });

    test('KPSS lag rule: re-applied to T of each step and reported', () {
      final dr = autoArima(nile, maxP: 1).differencing; // T = 100, then 99
      expect(dr.kpssLag, KpssLag.short);
      expect([for (final k in dr.kpss) k.lag], [2, 2]);
      expect([for (final k in dr.kpss) k.length], [100, 99]);
      final l4 = autoArima(nile, maxP: 1, kpssLag: KpssLag.l4).differencing;
      expect([for (final k in l4.kpss) k.lag], [4, 3]);
      expect(l4.kpssLag, KpssLag.l4);
    });

    test('period 0 and period s give the same likelihood (spec §3.4)', () {
      // A candidate without seasonal terms at D = 0 is fitted with period 0;
      // the same model declared with period s must be identical.
      final lynx = loadReferenceSeries('lynx');
      final cases = <(String, List<double>, int, int, bool)>[
        ('air_log d=1', air, 12, 1, false),
        ('air_log d=1 drift', air, 12, 1, true),
        ('us_acc d=0 mean', loadReferenceSeries('us_acc_deaths'), 12, 0, true),
        ('us_acc d=0', loadReferenceSeries('us_acc_deaths'), 12, 0, false),
        ('lynx d=0 mean', lynx, 10, 0, true),
        ('nile d=1', nile, 4, 1, false),
      ];
      for (final (label, y, s, d, c) in cases) {
        for (final (p, q) in [(1, 1), (2, 0), (0, 2)]) {
          SarimaxFitResult fit(int period) => fitSarimax(
            y,
            const {},
            order: ArimaOrder(p: p, d: d, q: q, seasonalPeriod: period),
            includeMean: c && d == 0,
            includeDrift: c && d == 1,
          );
          final a = fit(0).logLikelihood!;
          final b = fit(s).logLikelihood!;
          expect(
            (a - b).abs(),
            lessThanOrEqualTo(1e-10 * math.max(1, a.abs())),
            reason: '$label ($p,$d,$q): $a vs $b',
          );
        }
      }
    });

    test('BIC and AIC', () {
      for (final c in [InformationCriterion.aic, InformationCriterion.bic]) {
        final r = autoArima(nile, criterion: c);
        expect(r.criterion, c);
        expect(
          r.criterionValue,
          c.evaluate(
            logLikelihood: r.fit.logLikelihood!,
            parameters: r.trace
                .firstWhere((t) => identical(t.fit, r.fit))
                .parameterCount,
            observations: 99,
          ),
        );
      }
    });

    test('maxModels truncates and says so', () {
      final r = autoArima(nile, maxModels: 3);
      expect(r.truncated, isTrue);
      expect(r.fitsEvaluated, 3);
      _checkJournal(r.trace, 3, n: 100, differencing: r.differencing);
    });

    test('exhaustive grid: lexicographic, sum bound, both constants', () {
      final r = autoArima(
        nile,
        search: ArimaSearch.exhaustive,
        maxP: 2,
        maxQ: 2,
        maxOrderSum: 3,
      );
      final keys = [
        for (final t in r.trace) (t.order.p, t.order.q, t.constant),
      ];
      expect(keys, [
        for (var p = 0; p <= 2; p++)
          for (var q = 0; q <= 2; q++)
            if (p + q <= 3) ...[(p, q, false), (p, q, true)],
      ]);
      expect(r.trace.where((t) => t.becameCurrent), hasLength(1));
      _checkJournal(
        r.trace,
        r.fitsEvaluated,
        n: 100,
        differencing: r.differencing,
      );
    });

    test('stepwise neighbourhood: 17 neighbours in the order of rev. 9', () {
      List<(int, int, int, int, bool)> expected(
        (int, int, int, int, bool) k, {
        required bool seasonal,
        required bool constantAllowed,
      }) {
        final (p, q, sp, sq, c) = k;
        return [
          if (seasonal) ...[
            (p, q, sp - 1, sq, c),
            (p, q, sp, sq - 1, c),
            (p, q, sp + 1, sq, c),
            (p, q, sp, sq + 1, c),
            (p, q, sp - 1, sq - 1, c),
            (p, q, sp - 1, sq + 1, c),
            (p, q, sp + 1, sq - 1, c),
            (p, q, sp + 1, sq + 1, c),
          ],
          (p - 1, q, sp, sq, c),
          (p, q - 1, sp, sq, c),
          (p + 1, q, sp, sq, c),
          (p, q + 1, sp, sq, c),
          (p - 1, q - 1, sp, sq, c),
          (p - 1, q + 1, sp, sq, c),
          (p + 1, q - 1, sp, sq, c),
          (p + 1, q + 1, sp, sq, c),
          if (constantAllowed) (p, q, sp, sq, !c),
        ].where((n) {
          final (p, q, sp, sq, _) = n;
          return p >= 0 &&
              p <= 5 &&
              q >= 0 &&
              q <= 5 &&
              sp >= 0 &&
              sp <= (seasonal ? 2 : 0) &&
              sq >= 0 &&
              sq <= (seasonal ? 2 : 0);
        }).toList();
      }

      (int, int, int, int, bool) key(CandidateRecord r) => (
        r.order.p,
        r.order.q,
        r.order.seasonalP,
        r.order.seasonalQ,
        r.constant,
      );

      for (final (y, s) in [
        (nile, 1),
        (loadReferenceSeries('lynx'), 1),
        (loadReferenceSeries('lake_huron'), 1),
        (air, 12),
      ]) {
        final r = autoArima(y, period: s);
        final constantAllowed = r.differencing.d + r.differencing.seasonalD < 2;
        var current = key(r.trace.firstWhere((t) => t.becameCurrent));
        final maxStep = r.trace.map((t) => t.step).reduce(math.max);
        for (var step = 1; step <= maxStep; step++) {
          final visited = [
            for (final t in r.trace)
              if (t.stage == CandidateStage.step && t.step == step) key(t),
          ];
          final want = expected(
            current,
            seasonal: s > 1,
            constantAllowed: constantAllowed,
          );
          expect(visited, want.sublist(0, visited.length));
          final moved = r.trace.where(
            (t) =>
                t.stage == CandidateStage.step &&
                t.step == step &&
                t.becameCurrent,
          );
          if (moved.isNotEmpty) {
            current = key(moved.single);
          } else {
            expect(visited, want, reason: 'last scan is complete');
          }
        }
      }
    });

    test('in-process determinism: identical journals', () {
      final a = autoArima(air, period: 12).toJournal();
      final b = autoArima(air, period: 12).toJournal();
      expect(a, b);
    });
  });

  group('arguments', () {
    test('invalid arguments are ArgumentErrors', () {
      expect(() => autoArima(const []), throwsArgumentError);
      expect(() => autoArima(nile, period: 0), throwsArgumentError);
      expect(() => autoArima(nile, seasonalD: 1), throwsArgumentError);
      expect(() => autoArima(nile, d: -1), throwsArgumentError);
      expect(() => autoArima(nile, maxSeasonalD: 2), throwsArgumentError);
      expect(() => autoArima(nile, kpssAlpha: 0.2), throwsArgumentError);
      expect(() => autoArima(nile, level: 1), throwsArgumentError);
      expect(() => autoArima(nile, rootMargin: 0.99), throwsArgumentError);
      expect(() => autoArima(nile, maxModels: 0), throwsArgumentError);
      expect(() => autoArima(nile, horizon: -1), throwsArgumentError);
    });
  });

  group('gate F: robustness', () {
    test('NaN in the input is an ArgumentError', () {
      expect(
        () => autoArima([...nile.take(50), double.nan, ...nile.skip(51)]),
        throwsArgumentError,
      );
    });

    Matcher constantSeriesRefusal({required int d}) => throwsA(
      isA<AutoArimaNoModelException>()
          .having(
            (e) => e.reason,
            'reason',
            AutoArimaNoModelReason.constantSeries,
          )
          .having((e) => e.trace, 'trace', isEmpty)
          .having((e) => e.fitsEvaluated, 'fits', 0)
          .having((e) => e.differencing.d, 'd', d)
          .having(
            (e) => e.differencing.dSource,
            'dSource',
            DifferencingSource.constantSeries,
          ),
    );

    test('constant series: refused before the search (spec §2.3)', () {
      expect(
        () => autoArima(List<double>.filled(50, 5)),
        constantSeriesRefusal(d: 0),
      );
      expect(
        () => autoArima(List<double>.filled(50, 0)),
        constantSeriesRefusal(d: 0),
      );
      // Constant up to rounding relative to the input scale.
      expect(
        () => autoArima([
          for (var t = 0; t < 50; t++) 1e6 + (t.isEven ? 1e-6 : 0),
        ]),
        constantSeriesRefusal(d: 0),
      );
    });

    test('linear ramp: one difference, then constant → refused', () {
      expect(
        () => autoArima([for (var t = 0; t < 50; t++) 2 + 0.5 * t]),
        constantSeriesRefusal(d: 1),
      );
      // Same with rounding noise from non-representable steps.
      expect(
        () => autoArima([for (var t = 0; t < 60; t++) 0.1 * t + 1000.3]),
        constantSeriesRefusal(d: 1),
      );
    });

    test('given d: constancy of the working series is still checked', () {
      expect(
        () => autoArima([for (var t = 0; t < 50; t++) 2 + 0.5 * t], d: 1),
        throwsA(
          isA<AutoArimaNoModelException>().having(
            (e) => e.reason,
            'reason',
            AutoArimaNoModelReason.constantSeries,
          ),
        ),
      );
    });

    test('N < 12: KPSS skipped, short candidates refused before fitting', () {
      final y = [for (var t = 0; t < 10; t++) math.sin(t * 1.3) + t * 0.1];
      final r = autoArima(y);
      expect(r.differencing.dSource, DifferencingSource.skippedShort);
      expect(
        r.trace.any((t) => t.verdict == CandidateVerdict.rejectedTooShort),
        isTrue,
      );
      _checkJournal(
        r.trace,
        r.fitsEvaluated,
        n: 10,
        differencing: r.differencing,
      );
    });

    test('too short for any candidate: exception carries the journal', () {
      final y = [1.0, 3.0, 2.0, 5.0, 4.0, 6.0];
      expect(
        () => autoArima(y),
        throwsA(
          isA<AutoArimaNoModelException>()
              .having((e) => e.trace, 'trace', isNotEmpty)
              .having(
                (e) => e.trace.every(
                  (t) =>
                      t.verdict == CandidateVerdict.rejectedTooShort ||
                      t.verdict == CandidateVerdict.cachedDuplicate,
                ),
                'all too short',
                isTrue,
              )
              .having((e) => e.fitsEvaluated, 'fits', 0)
              .having(
                (e) => e.differencing.dSource,
                'dSource',
                DifferencingSource.skippedShort,
              ),
        ),
      );
    });

    test('s = 288: seasonal ARMA terms rejected as too large', () {
      final w = simulateArima(n: 600, seed: 1, ar: [0.5]);
      final y = [
        for (var t = 0; t < 600; t++)
          w[t] + 5 * math.sin(2 * math.pi * t / 288),
      ];
      final r = autoArima(y, period: 288);
      expect(
        r.trace.where((t) => t.verdict == CandidateVerdict.rejectedTooLarge),
        isNotEmpty,
      );
      expect(r.order.seasonalP + r.order.seasonalQ, 0);
      _checkJournal(
        r.trace,
        r.fitsEvaluated,
        n: 600,
        differencing: r.differencing,
      );
    });

    // Short seasonal series around the length rule of spec §5: no
    // ArgumentError of the fit path may surface, ever.
    for (final s in [4, 12]) {
      for (final sD in [0, 1]) {
        test(
          'short seasonal series, s = $s, D = $sD, lengths around the rule',
          () {
            final lengths = s == 12
                ? [for (var n = 18; n <= 40; n++) n]
                : [for (var n = 8; n <= 30; n++) n];
            for (final n in lengths) {
              final y = simulateArima(
                n: n,
                seed: 100 + n,
                seasonalMa: [0.5],
                ar: [0.3],
                period: s,
                seasonalD: sD,
              );
              for (final mode in ArimaSearch.values) {
                try {
                  final r = autoArima(
                    y,
                    period: s,
                    seasonalD: sD,
                    search: mode,
                    maxP: 2,
                    maxQ: 2,
                  );
                  _checkJournal(
                    r.trace,
                    r.fitsEvaluated,
                    n: n,
                    differencing: r.differencing,
                  );
                } on AutoArimaNoModelException catch (e) {
                  _checkJournal(
                    e.trace,
                    e.fitsEvaluated,
                    n: n,
                    differencing: e.differencing,
                  );
                }
              }
            }
          },
        );
      }
    }

    test('length rule at threshold − 1, threshold, threshold + 1', () {
      // (p, q, P, Q, c) at s = 12: rule binding values computed by hand.
      // (0,0,0,1,c=0): K = 2; T ≥ 10, T − 3 ≥ 1, T > max(0, 12) + 2 = 14.
      const k1 = (p: 0, q: 0, sp: 0, sq: 1, c: false);
      expect(lengthRuleViolation(k1, t: 14, period: 12), isNotNull);
      expect(lengthRuleViolation(k1, t: 15, period: 12), isNull);
      expect(lengthRuleViolation(k1, t: 16, period: 12), isNull);
      // (1,1,0,0,c=1): K = 4; T ≥ 20 binds.
      const k2 = (p: 1, q: 1, sp: 0, sq: 0, c: true);
      expect(lengthRuleViolation(k2, t: 19, period: 1), isNotNull);
      expect(lengthRuleViolation(k2, t: 20, period: 1), isNull);
      // (0,0,0,0,c=0): K = 1; T ≥ 8 binds.
      const k3 = (p: 0, q: 0, sp: 0, sq: 0, c: false);
      expect(lengthRuleViolation(k3, t: 7, period: 1), isNotNull);
      expect(lengthRuleViolation(k3, t: 8, period: 1), isNull);
      // (2,0,2,0,c=0) at s = 4: K = 5; T ≥ 40 binds before T > 10 + 5.
      const k4 = (p: 2, q: 0, sp: 2, sq: 0, c: false);
      expect(lengthRuleViolation(k4, t: 39, period: 4), isNotNull);
      expect(lengthRuleViolation(k4, t: 40, period: 4), isNull);

      // N = 22, s = 12, D = 1 → T = 10: seasonal AR (T ≤ 12 + K) is refused,
      // the small non-seasonal forms pass.
      final y = simulateArima(
        n: 22,
        seed: 7,
        seasonalMa: [0.5],
        period: 12,
        seasonalD: 1,
      );
      for (final key in [
        (p: 0, q: 0, sp: 0, sq: 0, c: false),
        (p: 0, q: 0, sp: 0, sq: 0, c: true),
        (p: 1, q: 0, sp: 0, sq: 0, c: false),
        (p: 0, q: 0, sp: 1, sq: 0, c: false),
      ]) {
        final out = evaluateCandidate(
          Float64List.fromList(y),
          key,
          d: 0,
          seasonalD: 1,
          period: 12,
          criterion: InformationCriterion.aicc,
          horizon: 2,
          rootMargin: 1.01,
        );
        final passes = key.sp == 0;
        expect(
          out.verdict == CandidateVerdict.rejectedTooShort,
          !passes,
          reason: '$key: ${out.verdict} ${out.detail}',
        );
        expect(out.detail ?? '', isNot(contains('DEFECT')));
      }
    });
  });
}
