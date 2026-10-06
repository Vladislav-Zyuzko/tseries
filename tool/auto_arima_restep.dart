// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// Re-runs only our STEPWISE search on an existing export (synthetic study or
// D′ windows) and writes a new export directory: index.csv, series.csv,
// holdout.csv and ours_kpss.csv are copied (differencing does not depend on
// the search), ours.csv gets the new stepwise columns (the exhaustive
// columns are copied: the grid does not depend on the neighbourhood), and
// ours_windows.csv (if present) gets the new mase_step. Also writes
// fits.csv — `series_id,fits,truncated` of the stepwise search at the
// default maxModels — for the cost comparison of spec §16.6 item 4.
//
// Usage (from packages/tseries):
//   dart run tool/auto_arima_restep.dart --from <export> --to <new export>
import 'dart:io';

import 'package:tseries/tseries.dart';

import '../test/support/auto_arima_synthetic.dart';

List<List<String>> _rows(String path) => File(path)
    .readAsLinesSync()
    .where((l) => l.isNotEmpty)
    .map((l) => l.split(','))
    .toList();

void main(List<String> args) {
  String arg(String name) => args[args.indexOf(name) + 1];
  final from = arg('--from');
  final to = arg('--to');
  String f(String dir, String name) => '$dir${Platform.pathSeparator}$name';
  Directory(to).createSync(recursive: true);
  for (final name in [
    'index.csv',
    'series.csv',
    'holdout.csv',
    'ours_kpss.csv',
  ]) {
    File(f(from, name)).copySync(f(to, name));
  }

  final index = _rows(f(from, 'index.csv'));
  final head = index.first;
  int col(String c) => head.indexOf(c);
  final train = <String, List<double>>{};
  for (final line in File(f(from, 'series.csv')).readAsLinesSync().skip(1)) {
    final id = line.substring(0, line.indexOf(','));
    (train[id] ??= []).add(
      double.parse(line.substring(line.lastIndexOf(',') + 1)),
    );
  }
  final test = <String, List<double>>{};
  for (final r in _rows(f(from, 'holdout.csv')).skip(1)) {
    (test[r[0]] ??= []).add(double.parse(r[2]));
  }
  final oldOurs = {for (final r in _rows(f(from, 'ours.csv')).skip(1)) r[0]: r};

  final ours = LineSink(
    f(to, 'ours.csv'),
    header:
        'series_id,d_star,D_star,step_p,step_q,step_P,step_Q,step_c,'
        'step_aicc,ex_p,ex_q,ex_P,ex_Q,ex_c,ex_aicc,mase_step,no_model',
  );
  final fits = LineSink(f(to, 'fits.csv'), header: 'series_id,fits,truncated');
  final newMase = <String, String>{};
  var done = 0;
  for (final r in index.skip(1)) {
    final id = r[0];
    final period = int.parse(r[col('period')]);
    final h = int.parse(r[col('holdout')]);
    final y = train[id]!;
    final old = oldOurs[id]!;
    final mm = period > 1 ? period : 1;
    try {
      final s = autoArima(y, period: period, horizon: h);
      final ms = g17(mase(y, test[id]!, s.forecast, mm));
      newMase[id] = ms;
      ours.writeln(
        '$id,${s.differencing.d},${s.differencing.seasonalD},'
        '${s.order.p},${s.order.q},${s.order.seasonalP},${s.order.seasonalQ},'
        '${s.constant ? 1 : 0},${g17(s.criterionValue)},'
        '${old.sublist(9, 15).join(',')},$ms,0',
      );
      fits.writeln('$id,${s.fitsEvaluated},${s.truncated ? 1 : 0}');
    } on AutoArimaNoModelException catch (e) {
      ours.writeln(old.join(','));
      fits.writeln('$id,${e.fitsEvaluated},${e.truncated ? 1 : 0}');
    }
    if (++done % 200 == 0) stdout.writeln('$done series');
  }
  ours.close();
  fits.close();

  final winFile = File(f(from, 'ours_windows.csv'));
  if (winFile.existsSync()) {
    final lines = winFile.readAsLinesSync();
    final h = lines.first.split(',');
    final c = h.indexOf('mase_step');
    final out = StringBuffer()..writeln(lines.first);
    for (final l in lines.skip(1).where((l) => l.isNotEmpty)) {
      final v = l.split(',');
      if (newMase[v[0]] != null) v[c] = newMase[v[0]]!;
      out.writeln(v.join(','));
    }
    File(f(to, 'ours_windows.csv')).writeAsStringSync(out.toString());
  }
  stdout.writeln('done: $done series');
}
