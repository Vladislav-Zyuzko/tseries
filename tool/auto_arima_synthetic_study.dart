/// Gate C of the autoArima specification (§9.C, revisions 5–6), the local
/// part: 100 seeded series per (model, n) cell (tagged `heavy`).
///
/// Computed here: the d̂ × d_true matrix at the default KPSS lag (short) and
/// at l4 (C-d, nominal gate ≥ 89 % on white noise, MA(1), (0,1,1),
/// (0,1,1)+drift); D for the seasonal model (d conditional on a correct D);
/// stepwise = exhaustive share and the hard "step < ex ⇒ winner outside the
/// grid" check; the MASE components "search" (stepwise vs the true order at
/// our d*, D*: median of per-series ratios ≤ 1.03, bootstrap 95 % upper bound
/// ≤ 1.08) and "differencing" (true order at d*, D* vs at d_true: report).
/// All shares carry 95 % Wilson intervals.
///
/// The parts that need R on the same series (C-d′, C-S′, MASE against the R
/// order) are computed by tool/auto_arima_r_compare.dart from the export
/// written here when AUTO_ARIMA_C_EXPORT names a directory (format:
/// doc/auto_arima_r_exchange.md).
///
/// Usage (from packages/tseries):
///   dart run tool/auto_arima_synthetic_study.dart [--reps 100]
///       [--models wn,ar1,...] [--out table.md] [--export <dir>]
/// The last output line starts with `gate notes:`.
library;

import 'dart:io';

import 'package:tseries/tseries.dart';

import '../test/support/auto_arima_synthetic.dart';

class _Cell {
  _Cell(this.m, this.n);
  final SyntheticModel m;
  final int n;
  var reps = 0, valid = 0, noModel = 0;
  final dShort = <int, int>{}, dL4 = <int, int>{};
  var dOkShort = 0, dOkL4 = 0, sDOk = 0, dGivenD = 0;
  var same = 0, stepBetter = 0, fitsStep = 0, fitsEx = 0;
  final stepBetterDefects = <String>[];
  final searchRatios = <double>[], diffRatios = <double>[];
}

ArimaOrder _trueOrder(SyntheticModel m, int d, int sD) => ArimaOrder(
  p: m.ar.length,
  d: d,
  q: m.ma.length,
  seasonalQ: m.sma.length,
  seasonalD: sD,
  seasonalPeriod: m.sma.isNotEmpty || sD > 0 ? m.period : 0,
);

List<double>? _trueForecast(
  SyntheticModel m,
  List<double> train,
  int d,
  int sD,
  int h,
) {
  try {
    final c = m.constantAt(d, sD);
    return fitSarimax(
      train,
      const {},
      order: _trueOrder(m, d, sD),
      horizon: h,
      includeMean: c && d + sD == 0,
      includeDrift: c && d + sD == 1,
    ).forecast;
  } on Object {
    return null;
  }
}

String _dist(Map<int, int> m) => [
  for (final d in [0, 1, 2]) '$d:${m[d] ?? 0}',
].join(' ');

String _ratio(List<double> xs) {
  if (xs.isEmpty) return '—';
  final b = bootstrapMedian(xs);
  return '${b.median.toStringAsFixed(3)} '
      '[${b.lo.toStringAsFixed(3)}, ${b.hi.toStringAsFixed(3)}]';
}

void main(List<String> args) {
  String? arg(String name) {
    final i = args.indexOf(name);
    return i >= 0 && i + 1 < args.length ? args[i + 1] : null;
  }

  final reps = int.tryParse(arg('--reps') ?? '') ?? 100;
  final only = arg('--models')?.split(',').toSet();
  final exportDir = arg('--export');
  final outPath = arg('--out');
  final cells = <_Cell>[];
  LineSink? sIndex, sTrain, sHoldout, sOurs, sKpss;

  if (exportDir != null) {
    String f(String name) => '$exportDir${Platform.pathSeparator}$name';
    sIndex = LineSink(
      f('index.csv'),
      header:
          'series_id,model,n,rep,seed,period,d_star,D_star,d_true,D_true,'
          'holdout',
    );
    sTrain = LineSink(
      f('series.csv'),
      header: 'series_id,model,n,rep,seed,period,d_star,D_star,t,value',
    );
    sHoldout = LineSink(f('holdout.csv'), header: 'series_id,t,value');
    sOurs = LineSink(
      f('ours.csv'),
      header:
          'series_id,d_star,D_star,step_p,step_q,step_P,step_Q,step_c,'
          'step_aicc,ex_p,ex_q,ex_P,ex_Q,ex_c,ex_aicc,mase_step,no_model',
    );
    sKpss = LineSink(
      f('ours_kpss.csv'),
      header: 'series_id,step,T,lag,eta,rejects',
    );
  }

  for (final m in syntheticModels) {
    if (only != null && !only.contains(m.id)) continue;
    for (final n in syntheticLengths) {
      {
        final sw = Stopwatch()..start();
        final cell = _Cell(m, n);
        cells.add(cell);
        final h = m.holdout;
        final mm = m.period > 1 ? m.period : 1;
        for (var rep = 0; rep < reps; rep++) {
          cell.reps++;
          final id = syntheticId(m, n, rep);
          final seed = syntheticSeed(m, n, rep);
          final (:train, :test) = syntheticSeries(m, n, rep);

          // d̂ at l4 (differencing only: no ARMA terms searched).
          DifferencingReport l4;
          try {
            l4 = autoArima(
              train,
              period: m.period,
              kpssLag: KpssLag.l4,
              maxP: 0,
              maxQ: 0,
              maxSeasonalP: 0,
              maxSeasonalQ: 0,
            ).differencing;
          } on AutoArimaNoModelException catch (e) {
            l4 = e.differencing;
          }
          cell.dL4[l4.d] = (cell.dL4[l4.d] ?? 0) + 1;
          if (l4.d == m.d) cell.dOkL4++;

          AutoArimaResult? step;
          AutoArimaResult? ex;
          DifferencingReport dr;
          try {
            step = autoArima(train, period: m.period, horizon: h);
            ex = autoArima(
              train,
              period: m.period,
              horizon: h,
              search: ArimaSearch.exhaustive,
              maxModels: 1000,
            );
            dr = step.differencing;
          } on AutoArimaNoModelException catch (e) {
            dr = e.differencing;
            cell.noModel++;
          }
          final dStar = dr.d, sDStar = dr.seasonalD;
          cell.dShort[dStar] = (cell.dShort[dStar] ?? 0) + 1;
          if (dStar == m.d) cell.dOkShort++;
          if (sDStar == m.sD) {
            cell.sDOk++;
            if (dStar == m.d) cell.dGivenD++;
          }

          double? masStep;
          if (step != null && ex != null) {
            cell.valid++;
            cell.fitsStep += step.fitsEvaluated;
            cell.fitsEx += ex.fitsEvaluated;
            final delta = step.criterionValue - ex.criterionValue;
            if (delta.abs() <= 1e-6) cell.same++;
            if (delta < -1e-6) {
              cell.stepBetter++;
              final o = step.order;
              if (o.p + o.q + o.seasonalP + o.seasonalQ <=
                  autoArimaDefaultMaxOrderSum) {
                cell.stepBetterDefects.add('$id ${step.order} Δ=$delta');
              }
            }
            masStep = mase(train, test, step.forecast, mm);
            final fTrueStar = _trueForecast(m, train, dStar, sDStar, h);
            final fTrue = _trueForecast(m, train, m.d, m.sD, h);
            if (fTrueStar != null) {
              final mt = mase(train, test, fTrueStar, mm);
              cell.searchRatios.add(masStep / mt);
              if (fTrue != null) {
                cell.diffRatios.add(mt / mase(train, test, fTrue, mm));
              }
            }
          }

          if (exportDir != null) {
            sIndex!.writeln(
              '$id,${m.id},$n,$rep,$seed,${m.period},$dStar,$sDStar,'
              '${m.d},${m.sD},$h',
            );
            final pre = '$id,${m.id},$n,$rep,$seed,${m.period},$dStar,$sDStar';
            for (var t = 0; t < train.length; t++) {
              sTrain!.writeln('$pre,${t + 1},${g17(train[t])}');
            }
            for (var t = 0; t < test.length; t++) {
              sHoldout!.writeln('$id,${n + t + 1},${g17(test[t])}');
            }
            String o(AutoArimaResult? r) => r == null
                ? ',,,,,'
                : '${r.order.p},${r.order.q},${r.order.seasonalP},'
                      '${r.order.seasonalQ},${r.constant ? 1 : 0},'
                      '${g17(r.criterionValue)}';
            sOurs!.writeln(
              '$id,$dStar,$sDStar,${o(step)},${o(ex)},'
              '${masStep == null ? '' : g17(masStep)},'
              '${step == null ? 1 : 0}',
            );
            for (var i = 0; i < dr.kpss.length; i++) {
              final k = dr.kpss[i];
              sKpss!.writeln(
                '$id,${i + 1},${k.length},${k.lag},${g17(k.statistic)},'
                '${k.rejectsStationarity ? 1 : 0}',
              );
            }
          }
        }
        stdout.writeln(
          '${m.id} n=$n: d short ${_dist(cell.dShort)} | l4 ${_dist(cell.dL4)}'
          ' | step=ex ${cell.same}/${cell.valid} | search '
          '${_ratio(cell.searchRatios)}',
        );
        stdout.writeln('  (${sw.elapsed.inSeconds} s)');
      }
    }
  }

  {
    for (final s in [sIndex, sTrain, sHoldout, sOurs, sKpss]) {
      s?.close();
    }
    final notes = <String>[];
    final rows = <String>[
      '| model | n | reps | d̂ short (0:1:2) | d ok short % [W] | '
          'd̂ l4 | d ok l4 % [W] | D ok % [W] | d ok given D % [W] | '
          'step=ex % [W] | step<ex (defects) | MASE search med [CI] | '
          'MASE diff med [CI] | fits step/ex | no model |',
      '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|',
    ];
    for (final c in cells) {
      final m = c.m;
      rows.add(
        '| ${m.label} | ${c.n} | ${c.reps} | ${_dist(c.dShort)} | '
        '${pctCi(c.dOkShort, c.reps)} | ${_dist(c.dL4)} | '
        '${pctCi(c.dOkL4, c.reps)} | '
        '${m.period > 1 ? pctCi(c.sDOk, c.reps) : '—'} | '
        '${m.period > 1 ? pctCi(c.dGivenD, c.sDOk) : '—'} | '
        '${pctCi(c.same, c.valid)} | ${c.stepBetter} '
        '(${c.stepBetterDefects.length}) | ${_ratio(c.searchRatios)} | '
        '${_ratio(c.diffRatios)} | '
        '${c.valid == 0 ? '—' : (c.fitsStep / c.valid).toStringAsFixed(1)}/'
        '${c.valid == 0 ? '—' : (c.fitsEx / c.valid).toStringAsFixed(1)} | '
        '${c.noModel} |',
      );
      if (m.nominalGate && c.dOkShort < 0.89 * c.reps) {
        notes.add('C-d nominal: ${m.id} n=${c.n} d ok ${c.dOkShort}/${c.reps}');
      }
      if (c.searchRatios.isNotEmpty) {
        final b = bootstrapMedian(c.searchRatios);
        if (b.median > 1.03 || b.hi > 1.08) {
          notes.add(
            'MASE search: ${m.id} n=${c.n} median '
            '${b.median.toStringAsFixed(3)} CI hi ${b.hi.toStringAsFixed(3)}',
          );
        }
      }
      for (final d in c.stepBetterDefects) {
        notes.add('step<ex inside grid: $d');
      }
    }
    final text = [
      ...rows,
      '',
      'gate notes: ${notes.isEmpty ? 'none' : notes.join('; ')}',
    ].join('\n');
    stdout.writeln(text);
    if (outPath != null) File(outPath).writeAsStringSync(text);
  }
}
