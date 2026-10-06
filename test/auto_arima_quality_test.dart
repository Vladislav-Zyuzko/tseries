/// Gate B of the autoArima specification (§9, revision 2): quality of the
/// order search on the six reference series against the orders R's
/// auto.arima chose (external tool; only the orders are used).
///
/// Every AICc here is ours: a reference order o is evaluated by our own fit
/// at OUR (d*, D*) — AICc_ours(o). R's d/D agreement is reported, not gated.
library;

import 'dart:typed_data';

import 'package:test/test.dart';
import 'package:tseries/src/domain/internal/auto_arima_candidate.dart';
import 'package:tseries/tseries.dart';

import 'support/auto_arima_fixtures.dart';

const _eps = 1e-6;

void main() {
  final rOrders = loadROrders();
  final rows = <String>[];
  var stepWithin2 = 0;
  var stepCompared = 0;
  final stepOver4 = <String>[];
  final exhaustiveFailures = <String>[];

  for (final (name, period) in autoArimaReferenceSeries) {
    test(name, () {
      final y = loadReferenceSeries(name);
      final step = autoArima(y, period: period);
      final ex = autoArima(
        y,
        period: period,
        search: ArimaSearch.exhaustive,
        maxModels: 1000,
      );
      expect(ex.truncated, isFalse);
      expect(step.truncated, isFalse);
      final dd = ex.differencing.d;
      final sD = ex.differencing.seasonalD;
      expect(step.differencing.d, dd);
      expect(step.differencing.seasonalD, sD);

      /// AICc_ours of an R order at our (d*, D*), or null with the reason
      /// when it did not pass our rules.
      (double?, String) oursOf(ROrder o) {
        final c = o.constant != 'none' && dd + sD < 2;
        final key = (p: o.p, q: o.q, sp: o.sp, sq: o.sq, c: c);
        // The exhaustive journal already holds it when it is in the grid.
        for (final r in ex.trace) {
          if (r.duplicateOf != null) continue;
          final ro = r.order;
          if (ro.p == key.p &&
              ro.q == key.q &&
              ro.seasonalP == key.sp &&
              ro.seasonalQ == key.sq &&
              r.constant == key.c) {
            return r.verdict == CandidateVerdict.accepted
                ? (r.criterion, 'grid #${r.index}')
                : (null, '${r.verdict.name} (grid #${r.index})');
          }
        }
        final out = evaluateCandidate(
          Float64List.fromList(y),
          key,
          d: dd,
          seasonalD: sD,
          period: period,
          criterion: InformationCriterion.aicc,
          horizon: 0,
          rootMargin: 1.01,
        );
        return out.verdict == CandidateVerdict.accepted
            ? (out.criterion, 'off-grid fit')
            : (null, '${out.verdict.name} (off-grid)');
      }

      final rStep = rOrders['$name|stepwise']!;
      final rEx = rOrders['$name|exhaustive']!;
      final (aStep, whyStep) = oursOf(rStep);
      final (aEx, whyEx) = oursOf(rEx);

      // Exhaustive may not lose to a candidate of its own grid.
      for (final (a, label) in [(aEx, 'o_R^ex'), (aStep, 'o_R^step')]) {
        if (a != null && ex.criterionValue > a + _eps) {
          exhaustiveFailures.add(
            '$name: exhaustive ${ex.criterionValue} > AICc_ours($label) $a',
          );
        }
      }
      double? delta;
      if (aStep != null) {
        delta = step.criterionValue - aStep;
        stepCompared++;
        if (delta <= 2) stepWithin2++;
        if (delta > 4) stepOver4.add('$name Δ=$delta');
      }

      String fmt(double? v) => v == null ? '—' : v.toStringAsFixed(3);
      final dMatch =
          (rStep.d == dd ? 'd=' : 'd≠') + (rStep.sd == sD ? 'D=' : 'D≠');
      rows.add(
        '| $name | ($dd,$sD) vs R (${rStep.d},${rStep.sd}) $dMatch '
        '| ${step.order} c=${step.constant ? 1 : 0} ${fmt(step.criterionValue)} '
        '[${step.fitsEvaluated}] '
        '| ${ex.order} c=${ex.constant ? 1 : 0} ${fmt(ex.criterionValue)} '
        '[${ex.fitsEvaluated}] '
        '| ${describeROrder(rStep)} → ${fmt(aStep)} ($whyStep) '
        '| ${describeROrder(rEx)} → ${fmt(aEx)} ($whyEx) '
        '| ${fmt(delta)} | ${fmt(step.criterionValue - ex.criterionValue)} |',
      );
    }, timeout: const Timeout(Duration(minutes: 3)));
  }

  tearDownAll(() {
    // ignore: avoid_print
    print(
      [
        '| series | (d*,D*) vs R | ours stepwise AICc [fits] | '
            'ours exhaustive AICc [fits] | o_R^step → AICc_ours | '
            'o_R^ex → AICc_ours | Δ_step | step − ex |',
        ...rows,
        'stepwise: Δ ≤ 2 on $stepWithin2 of $stepCompared compared; '
            'Δ > 4: $stepOver4; exhaustive failures: $exhaustiveFailures',
      ].join('\n'),
    );
  });

  test('gates', () {
    expect(exhaustiveFailures, isEmpty, reason: 'exhaustive gate');
    expect(stepOver4, isEmpty, reason: 'stepwise: nowhere above +4');
    expect(
      stepWithin2,
      greaterThanOrEqualTo(stepCompared - 1),
      reason: 'stepwise: Δ ≤ +2 on all but at most one of 6',
    );
    expect(stepCompared, greaterThanOrEqualTo(5));
  });
}
