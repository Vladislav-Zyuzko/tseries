// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// Gate D′ of the autoArima specification (§9.D, §15.2, §15.4): exports the
// rolling-origin training windows of the six reference series for the
// external R run, together with our own per-window results.
//
// Windows: the last 10 origins that leave a full horizon, h = s for a
// seasonal series and 10 otherwise; the window is y[0 .. origin), the
// holdout y[origin .. origin + h). series_id = `<series>_o<k>`, k = 1..10
// (k = 10 is the longest window). At each window our autoArima is run anew
// in both modes.
//
// Files (same formats as the synthetic export, doc/auto_arima_r_exchange.md;
// in series.csv `model` = the series name, `n` = the window length,
// `rep` = k, `seed` = 0):
//   index.csv, series.csv, holdout.csv, ours.csv, ours_kpss.csv, and
//   ours_windows.csv — per window: k, cycles, D source, F_S, threshold, our
//   (d, D), MASE of our stepwise and exhaustive choice, MASE of our stepwise
//   with d forced to 1 (wwwusage: the price of d = 0), MASE of our stepwise
//   with seasonalD: 1 (seasonal series), and the old gate-D comparison: MASE
//   of the full-series R order fitted at our (d*, D*) of the window.
//
// Usage (from packages/tseries):
//   dart run tool/auto_arima_rolling_windows.dart --export <dir>
import 'dart:io';

import 'package:tseries/tseries.dart';

import '../test/support/auto_arima_fixtures.dart';
import '../test/support/auto_arima_synthetic.dart';

double? _maseOf(
  List<double> train,
  List<double> test,
  int mm,
  AutoArimaResult Function() run,
) {
  try {
    return mase(train, test, run().forecast, mm);
  } on AutoArimaNoModelException {
    return null;
  }
}

void main(List<String> args) {
  final i = args.indexOf('--export');
  if (i < 0 || i + 1 >= args.length) {
    stderr.writeln('usage: --export <dir>');
    exit(64);
  }
  final dir = args[i + 1];
  String f(String name) => '$dir${Platform.pathSeparator}$name';
  final sIndex = LineSink(
    f('index.csv'),
    header:
        'series_id,model,n,rep,seed,period,d_star,D_star,d_true,D_true,'
        'holdout',
  );
  final sTrain = LineSink(
    f('series.csv'),
    header: 'series_id,model,n,rep,seed,period,d_star,D_star,t,value',
  );
  final sHold = LineSink(f('holdout.csv'), header: 'series_id,t,value');
  final sOurs = LineSink(
    f('ours.csv'),
    header:
        'series_id,d_star,D_star,step_p,step_q,step_P,step_Q,step_c,'
        'step_aicc,ex_p,ex_q,ex_P,ex_Q,ex_c,ex_aicc,mase_step,no_model',
  );
  final sKpss = LineSink(
    f('ours_kpss.csv'),
    header: 'series_id,step,T,lag,eta,rejects',
  );
  final sWin = LineSink(
    f('ours_windows.csv'),
    header:
        'series_id,series,k,window_n,cycles,D_source,F_S,threshold,d_star,'
        'D_star,mase_step,mase_ex,mase_step_d1,mase_step_sD1,'
        'mase_oldR_step,mase_oldR_ex',
  );
  final rOrders = loadROrders();
  String g(double? v) => v == null ? '' : g17(v);

  for (final (name, period) in autoArimaReferenceSeries) {
    final y = loadReferenceSeries(name);
    final h = period > 1 ? period : 10;
    final mm = period > 1 ? period : 1;
    var k = 0;
    for (var origin = y.length - h - 9; origin <= y.length - h; origin++) {
      k++;
      final id = '${name}_o$k';
      final train = y.sublist(0, origin);
      final test = y.sublist(origin, origin + h);

      AutoArimaResult? step, ex;
      DifferencingReport dr;
      try {
        step = autoArima(train, period: period, horizon: h);
        ex = autoArima(
          train,
          period: period,
          horizon: h,
          search: ArimaSearch.exhaustive,
          maxModels: 1000,
        );
        dr = step.differencing;
      } on AutoArimaNoModelException catch (e) {
        dr = e.differencing;
      }
      final dStar = dr.d, sDStar = dr.seasonalD;

      // Old gate D: the full-series R order at our (d*, D*) of the window.
      double? oldR(String mode) {
        final o = rOrders['$name|$mode']!;
        final c = o.constant != 'none' && dStar + sDStar < 2;
        try {
          final fc = fitSarimax(
            train,
            const {},
            order: ArimaOrder(
              p: o.p,
              d: dStar,
              q: o.q,
              seasonalP: o.sp,
              seasonalD: sDStar,
              seasonalQ: o.sq,
              seasonalPeriod: o.sp + o.sq + sDStar > 0 ? period : 0,
            ),
            horizon: h,
            includeMean: c && dStar + sDStar == 0,
            includeDrift: c && dStar + sDStar == 1,
          ).forecast;
          return mase(train, test, fc, mm);
        } on Object {
          return null;
        }
      }

      final d1 = _maseOf(
        train,
        test,
        mm,
        () => autoArima(train, period: period, horizon: h, d: 1),
      );
      final sD1 = period > 1
          ? _maseOf(
              train,
              test,
              mm,
              () => autoArima(train, period: period, horizon: h, seasonalD: 1),
            )
          : null;

      sIndex.writeln(
        '$id,$name,${train.length},$k,0,$period,$dStar,$sDStar,,,$h',
      );
      final pre = '$id,$name,${train.length},$k,0,$period,$dStar,$sDStar';
      for (var t = 0; t < train.length; t++) {
        sTrain.writeln('$pre,${t + 1},${g17(train[t])}');
      }
      for (var t = 0; t < test.length; t++) {
        sHold.writeln('$id,${origin + t + 1},${g17(test[t])}');
      }
      String o(AutoArimaResult? r) => r == null
          ? ',,,,,'
          : '${r.order.p},${r.order.q},${r.order.seasonalP},'
                '${r.order.seasonalQ},${r.constant ? 1 : 0},'
                '${g17(r.criterionValue)}';
      final ms = step == null ? null : mase(train, test, step.forecast, mm);
      final me = ex == null ? null : mase(train, test, ex.forecast, mm);
      sOurs.writeln(
        '$id,$dStar,$sDStar,${o(step)},${o(ex)},${g(ms)},'
        '${step == null ? 1 : 0}',
      );
      for (var j = 0; j < dr.kpss.length; j++) {
        final kr = dr.kpss[j];
        sKpss.writeln(
          '$id,${j + 1},${kr.length},${kr.lag},${g17(kr.statistic)},'
          '${kr.rejectsStationarity ? 1 : 0}',
        );
      }
      sWin.writeln(
        '$id,$name,$k,${train.length},${dr.seasonalCycles ?? ''},'
        '${dr.seasonalDSource.name},${g(dr.seasonalStrength)},'
        '${g(dr.seasonalThreshold)},$dStar,$sDStar,${g(ms)},${g(me)},'
        '${g(d1)},${g(sD1)},${g(oldR('stepwise'))},${g(oldR('exhaustive'))}',
      );
      stdout.writeln(
        '$id n=${train.length} cycles=${dr.seasonalCycles ?? '-'} '
        'D=${dr.seasonalDSource.name} (d,D)=($dStar,$sDStar) '
        'step=${step?.order} ex=${ex?.order}',
      );
    }
  }
  for (final s in [sIndex, sTrain, sHold, sOurs, sKpss, sWin]) {
    s.close();
  }
}
