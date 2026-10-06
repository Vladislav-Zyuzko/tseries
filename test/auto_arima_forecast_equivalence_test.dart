/// Gate D of the autoArima specification (§9.D): forecast equivalence with
/// the R reference orders on the six reference series (tagged `heavy`).
///
/// Rolling origin, 10 origins (the last 10 that leave a full horizon), the
/// horizon h = s for a seasonal series and 10 otherwise. At each origin our
/// autoArima is re-run on the training part (same mode as the reference: our
/// stepwise against o_R^step, our exhaustive against o_R^ex), and o_R's
/// (p, q, P, Q, constant) is fitted by our fitSarimax at OUR (d*, D*) of that
/// origin, as in gate B. MASE (scale: in-sample (seasonal) naive) of each on
/// the next h values; gate: the median over origins of
/// MASE(ours) / MASE(o_R) is within [0.95, 1.05] for every series and mode.
@Tags(['heavy'])
library;

import 'package:test/test.dart';
import 'package:tseries/tseries.dart';

import 'support/auto_arima_fixtures.dart';
import 'support/auto_arima_synthetic.dart';

void main() {
  final rOrders = loadROrders();
  final rows = <String>[
    '| series | mode | h | origins | median ratio | min / max ratio | '
        'origins with the same order | our orders (count) |',
    '|---|---|---|---|---|---|---|---|',
  ];
  final failures = <String>[];

  for (final (name, period) in autoArimaReferenceSeries) {
    for (final mode in ArimaSearch.values) {
      test('$name ${mode.name}', () {
        final y = loadReferenceSeries(name);
        final h = period > 1 ? period : 10;
        final mm = period > 1 ? period : 1;
        final o =
            rOrders['$name|${mode == ArimaSearch.stepwise ? 'stepwise' : 'exhaustive'}']!;
        final ratios = <double>[];
        final orders = <String, int>{};
        var same = 0;
        for (var origin = y.length - h - 9; origin <= y.length - h; origin++) {
          final train = y.sublist(0, origin);
          final test = y.sublist(origin, origin + h);
          final AutoArimaResult ours;
          try {
            ours = autoArima(
              train,
              period: period,
              horizon: h,
              search: mode,
              maxModels: 1000,
            );
          } on AutoArimaNoModelException {
            continue;
          }
          final dd = ours.differencing.d;
          final sD = ours.differencing.seasonalD;
          final c = o.constant != 'none' && dd + sD < 2;
          final List<double> fr;
          try {
            fr = fitSarimax(
              train,
              const {},
              order: ArimaOrder(
                p: o.p,
                d: dd,
                q: o.q,
                seasonalP: o.sp,
                seasonalD: sD,
                seasonalQ: o.sq,
                seasonalPeriod: o.sp + o.sq + sD > 0 ? period : 0,
              ),
              horizon: h,
              includeMean: c && dd + sD == 0,
              includeDrift: c && dd + sD == 1,
            ).forecast;
          } on Object {
            continue;
          }
          final key = '${ours.order} c=${ours.constant ? 1 : 0}';
          orders[key] = (orders[key] ?? 0) + 1;
          final ro = ours.order;
          if (ro.p == o.p &&
              ro.q == o.q &&
              ro.seasonalP == o.sp &&
              ro.seasonalQ == o.sq &&
              ours.constant == c) {
            same++;
          }
          ratios.add(
            mase(train, test, ours.forecast, mm) / mase(train, test, fr, mm),
          );
        }
        final med = median(ratios);
        ratios.sort();
        rows.add(
          '| $name | ${mode.name} | $h | ${ratios.length} | '
          '${med.toStringAsFixed(4)} | ${ratios.first.toStringAsFixed(3)} / '
          '${ratios.last.toStringAsFixed(3)} | $same | '
          '${orders.entries.map((e) => '${e.key} (${e.value})').join('; ')} |',
        );
        if (!(med >= 0.95 && med <= 1.05)) {
          failures.add('$name ${mode.name}: median $med');
        }
        expect(ratios.length, greaterThanOrEqualTo(10));
      }, timeout: const Timeout(Duration(minutes: 20)));
    }
  }

  tearDownAll(() {
    // ignore: avoid_print
    print([...rows, '', 'failures: $failures'].join('\n'));
  });

  test('D: medians within [0.95, 1.05]', () {
    expect(failures, isEmpty);
  });
}
