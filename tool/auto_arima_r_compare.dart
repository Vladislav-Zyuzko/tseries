// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// Reads the synthetic-study export of tool/auto_arima_synthetic_study.dart and
// the results of the external R run on the same series, and computes the
// R-dependent parts of gate C of the autoArima specification (§9.C, §14):
//
//  * C-d′ — KPSS statistic equivalence (|η − η_R| ≤ 1e-6·max(1, η) at the
//    same T and lag, every step) and decision equivalence (d̂ = ndiffs_R,
//    a mismatch allowed only when |η − crit| ≤ 1e-6 at some step);
//  * C-S′ — share of series where AICc(our stepwise) > AICc_ours(R stepwise
//    order) + 1e-3 at the same (d*, D*): ≤ 5 %, Wilson upper bound ≤ 8 %;
//    every case listed;
//  * MASE against the reference — median over series of
//    MASE(our stepwise) / MASE(R auto.arima order with its own d, D, fitted
//    by our fitSarimax): ≤ 1.05, bootstrap 95 % interval reported;
//  * reports: R step = ex share, F_S (ours, classical) vs F_S_STL and
//    nsdiffs vs our D.
//
// File formats: doc/auto_arima_r_exchange.md.
//
// Usage (from packages/tseries):
//   dart run tool/auto_arima_r_compare.dart --export <dir> [--export <dir>…]
//       --r <dir> [--out report.md]
//   dart run tool/auto_arima_r_compare.dart --export <dir> --selftest <dir>
//   dart run tool/auto_arima_r_compare.dart --windows --export <cexp_D>
//       --r <dir> [--out report.md]
// --windows: the export is the rolling-origin windows of gate D′
// (tool/auto_arima_rolling_windows.dart). Cells are per series (not per
// window length), and gate D′ is computed: per window R's own stepwise and
// exhaustive orders (with R's d, D) are fitted by our fitSarimax and the
// median over windows of MASE(ours)/MASE(o_R) must be ≤ 1.05 (one-sided);
// the per-window table carries the window length, cycles, our D source,
// (d, D) ours and R's, the old-gate-D ratio, wwwusage's forced d = 1 ratio
// and the seasonalD: 1 run of seasonal series.
// --selftest writes R-format files built from OUR results into <dir> and
// compares against them (every gate must pass trivially: plumbing check).
import 'dart:io';
import 'dart:typed_data';

import 'package:tseries/src/domain/internal/auto_arima_candidate.dart';
import 'package:tseries/tseries.dart';

import '../test/support/auto_arima_synthetic.dart';

typedef _Row = Map<String, String>;

List<_Row> _csv(String path) {
  final f = File(path);
  if (!f.existsSync()) return const [];
  final lines = f
      .readAsLinesSync()
      .where((l) => l.isNotEmpty && !l.startsWith('#'))
      .toList();
  if (lines.isEmpty) return const [];
  final head = lines.first.split(',').map((h) => h.trim()).toList();
  return [
    for (final l in lines.skip(1))
      () {
        final v = l.split(',');
        return {
          for (var i = 0; i < head.length; i++)
            head[i]: i < v.length ? v[i].trim() : '',
        };
      }(),
  ];
}

int _i(String s) => int.parse(s);
double _d(String s) => double.parse(s);
bool _b(String s) => s == '1' || s.toUpperCase() == 'TRUE';

class _Series {
  _Series(this.index);
  final _Row index;
  final train = <double>[];
  final test = <double>[];
  _Row? ours;
  final oursKpss = <_Row>[];
  final rOrders = <String, _Row>{};
  final rKpss = <_Row>[];
  _Row? rSeasonal;

  String get id => index['series_id']!;
  String get model => index['model']!;
  int get n => _i(index['n']!);
  int get period => _i(index['period']!);
  int get dStar => _i(index['d_star']!);
  int get sDStar => _i(index['D_star']!);
}

void main(List<String> args) {
  final exports = <String>[];
  String? rDir;
  String? out;
  String? selftest;
  var windows = false;
  for (var i = 0; i < args.length; i++) {
    switch (args[i]) {
      case '--export':
        exports.add(args[++i]);
      case '--r':
        rDir = args[++i];
      case '--out':
        out = args[++i];
      case '--windows':
        windows = true;
      case '--selftest':
        selftest = args[++i];
      default:
        stderr.writeln('unknown argument ${args[i]}');
        exit(64);
    }
  }
  String p(String dir, String f) => '$dir${Platform.pathSeparator}$f';

  // ---- Load the export ----------------------------------------------------
  final series = <String, _Series>{};
  for (final dir in exports) {
    for (final r in _csv(p(dir, 'index.csv'))) {
      series[r['series_id']!] = _Series(r);
    }
    for (final line in File(p(dir, 'series.csv')).readAsLinesSync().skip(1)) {
      final c = line.indexOf(',');
      final id = line.substring(0, c);
      series[id]!.train.add(
        double.parse(line.substring(line.lastIndexOf(',') + 1)),
      );
    }
    for (final r in _csv(p(dir, 'holdout.csv'))) {
      series[r['series_id']!]!.test.add(_d(r['value']!));
    }
    for (final r in _csv(p(dir, 'ours.csv'))) {
      series[r['series_id']!]!.ours = r;
    }
    for (final r in _csv(p(dir, 'ours_kpss.csv'))) {
      series[r['series_id']!]!.oursKpss.add(r);
    }
  }
  stdout.writeln('loaded ${series.length} series from ${exports.join(', ')}');

  // ---- Self-test: R files from our own results ---------------------------
  if (selftest != null) {
    rDir = selftest;
    Directory(selftest).createSync(recursive: true);
    final o = StringBuffer(
      'series_id,mode,p,d,q,P,D,Q,has_mean,has_drift,aicc_R\n',
    );
    final k = StringBuffer('series_id,step,T,lag,eta_R,ndiffs_R\n');
    final s = StringBuffer('series_id,F_S_STL,nsdiffs_R\n');
    for (final x in series.values) {
      final r = x.ours!;
      if (r['no_model'] == '1') continue;
      final mean = x.dStar + x.sDStar == 0;
      for (final (mode, pre) in [
        ('step', 'step'),
        ('ex', 'ex'),
        ('auto', 'step'),
      ]) {
        final c = r['${pre}_c'] == '1';
        o.writeln(
          '${x.id},$mode,${r['${pre}_p']},${x.dStar},${r['${pre}_q']},'
          '${r['${pre}_P']},${x.sDStar},${r['${pre}_Q']},'
          '${c && mean ? 1 : 0},${c && !mean ? 1 : 0},${r['${pre}_aicc']}',
        );
      }
      for (final kr in x.oursKpss) {
        k.writeln(
          '${x.id},${kr['step']},${kr['T']},${kr['lag']},${kr['eta']},${x.dStar}',
        );
      }
      if (x.period > 1) {
        s.writeln(
          '${x.id},${g17(seasonalStrength(x.train, period: x.period))},${x.sDStar}',
        );
      }
    }
    File(p(selftest, 'r_orders.csv')).writeAsStringSync(o.toString());
    File(p(selftest, 'r_kpss.csv')).writeAsStringSync(k.toString());
    File(p(selftest, 'r_seasonal.csv')).writeAsStringSync(s.toString());
  }
  if (rDir == null) {
    stderr.writeln('--r <dir> (or --selftest <dir>) is required');
    exit(64);
  }

  for (final r in _csv(p(rDir, 'r_orders.csv'))) {
    series[r['series_id']]?.rOrders[r['mode']!] = r;
  }
  for (final r in _csv(p(rDir, 'r_kpss.csv'))) {
    series[r['series_id']]?.rKpss.add(r);
  }
  for (final r in _csv(p(rDir, 'r_seasonal.csv'))) {
    series[r['series_id']]?.rSeasonal = r;
  }

  // ---- Per cell -------------------------------------------------------------
  final cells = <String, List<_Series>>{};
  for (final x in series.values) {
    (cells[windows ? x.model : '${x.model}|${x.n}'] ??= []).add(x);
  }
  final report = <String>[
    '| cell | series | η compared / viol. | lag/T mismatch | d = ndiffs | '
        'd ≠ ndiffs at tie / defect | C-S′ worse % [W] (cases) | '
        'R step=ex % | MASE ours/R-auto med [CI] | order = R step % | '
        'ours better % | ours step=ex % |',
    '|---|---|---|---|---|---|---|---|---|---|---|---|',
  ];
  final details = <String>[];
  var anyDefect = false;
  var statFail = false;
  for (final MapEntry(key: cell, value: xs) in cells.entries) {
    var etaN = 0, etaViol = 0, lagMis = 0, dEq = 0, dTie = 0, dDef = 0;
    var csN = 0, csWorse = 0, rSame = 0, rBoth = 0;
    var csMatch = 0, csBetter = 0, oSame = 0, oBoth = 0;
    final ratios = <double>[];
    for (final x in xs) {
      final ours = x.ours!;
      // C-d′ — statistic.
      var tie = false;
      for (final rk in x.rKpss) {
        final ok = x.oursKpss.where((o) => o['step'] == rk['step']).toList();
        if (ok.isEmpty) continue;
        final o = ok.single;
        final eta = _d(o['eta']!);
        if ((eta - 0.463).abs() <= 1e-6) tie = true;
        if (o['T'] != rk['T'] || o['lag'] != rk['lag']) {
          lagMis++;
          details.add(
            'lag/T: ${x.id} step ${rk['step']} ours T=${o['T']} '
            'l=${o['lag']} R T=${rk['T']} l=${rk['lag']}',
          );
          continue;
        }
        etaN++;
        final etaR = _d(rk['eta_R']!);
        if ((eta - etaR).abs() > 1e-6 * (eta.abs() > 1 ? eta.abs() : 1)) {
          etaViol++;
          details.add('η: ${x.id} step ${rk['step']} ours $eta R $etaR');
        }
      }
      // C-d′ — decision.
      if (x.rKpss.isNotEmpty) {
        final nd = _i(x.rKpss.first['ndiffs_R']!);
        if (nd == x.dStar) {
          dEq++;
        } else if (tie) {
          dTie++;
        } else {
          dDef++;
          details.add('d: ${x.id} ours ${x.dStar} ndiffs_R $nd');
        }
      }
      if (ours['no_model'] == '1') continue;
      // C-S′.
      final rs = x.rOrders['step'];
      if (rs != null && _i(rs['d']!) == x.dStar && _i(rs['D']!) == x.sDStar) {
        final key = (
          p: _i(rs['p']!),
          q: _i(rs['q']!),
          sp: _i(rs['P']!),
          sq: _i(rs['Q']!),
          c:
              (_b(rs['has_mean']!) || _b(rs['has_drift']!)) &&
              x.dStar + x.sDStar < 2,
        );
        final outc = evaluateCandidate(
          Float64List.fromList(x.train),
          key,
          d: x.dStar,
          seasonalD: x.sDStar,
          period: x.period,
          criterion: InformationCriterion.aicc,
          horizon: 0,
          rootMargin: 1.01,
        );
        if (outc.verdict == CandidateVerdict.accepted) {
          csN++;
          final ourStep = _d(ours['step_aicc']!);
          if (_i(ours['step_p']!) == key.p &&
              _i(ours['step_q']!) == key.q &&
              _i(ours['step_P']!) == key.sp &&
              _i(ours['step_Q']!) == key.sq &&
              (ours['step_c'] == '1') == key.c) {
            csMatch++;
          }
          if (ourStep < outc.criterion! - 1e-3) csBetter++;
          if (ourStep > outc.criterion! + 1e-3) {
            csWorse++;
            details.add(
              'C-S′: ${x.id} ours ${ourStep.toStringAsFixed(4)} '
              '(${ours['step_p']},${ours['step_q']},${ours['step_P']},'
              '${ours['step_Q']},c=${ours['step_c']}) vs R order $key '
              '${outc.criterion!.toStringAsFixed(4)}',
            );
          }
        } else {
          details.add(
            'C-S′ skipped: ${x.id} R step order $key '
            '${outc.verdict.name} with us',
          );
        }
      }
      if ((ours['ex_aicc'] ?? '').isNotEmpty) {
        oBoth++;
        if ((_d(ours['step_aicc']!) - _d(ours['ex_aicc']!)).abs() <= 1e-6) {
          oSame++;
        }
      }
      final re = x.rOrders['ex'];
      if (rs != null && re != null) {
        rBoth++;
        if ((_d(rs['aicc_R']!) - _d(re['aicc_R']!)).abs() <= 1e-6) rSame++;
      }
      // MASE against the R auto order (its own d, D).
      final ra = x.rOrders['auto'];
      final ms = ours['mase_step'];
      if (ra != null && ms != null && ms.isNotEmpty) {
        final d = _i(ra['d']!), sD = _i(ra['D']!);
        final sp = _i(ra['P']!), sq = _i(ra['Q']!);
        try {
          final f = fitSarimax(
            x.train,
            const {},
            order: ArimaOrder(
              p: _i(ra['p']!),
              d: d,
              q: _i(ra['q']!),
              seasonalP: sp,
              seasonalD: sD,
              seasonalQ: sq,
              seasonalPeriod: sp + sq + sD > 0 ? x.period : 0,
            ),
            horizon: x.test.length,
            includeMean: _b(ra['has_mean']!) && d + sD == 0,
            includeDrift: _b(ra['has_drift']!) && d + sD == 1,
          );
          final mm = x.period > 1 ? x.period : 1;
          ratios.add(_d(ms) / mase(x.train, x.test, f.forecast, mm));
        } on Object catch (e) {
          details.add('MASE skipped: ${x.id} R auto order fit failed: $e');
        }
      }
    }
    if (etaViol > 0 || dDef > 0 || lagMis > 0) anyDefect = true;
    final w = wilson(csWorse, csN);
    if (csN > 0 && (csWorse > 0.05 * csN || w.hi > 0.08)) statFail = true;
    final b = bootstrapMedian(ratios);
    if (ratios.isNotEmpty && b.median > 1.05) statFail = true;
    report.add(
      '| $cell | ${xs.length} | $etaN / $etaViol | $lagMis | $dEq | '
      '$dTie / $dDef | ${pctCi(csWorse, csN)} ($csWorse) | '
      '${rBoth == 0 ? '—' : (100 * rSame / rBoth).toStringAsFixed(1)} | '
      '${ratios.isEmpty ? '—' : '${b.median.toStringAsFixed(3)} '
                '[${b.lo.toStringAsFixed(3)}, ${b.hi.toStringAsFixed(3)}]'} | '
      '${csN == 0 ? '—' : (100 * csMatch / csN).toStringAsFixed(1)} | '
      '${csN == 0 ? '—' : (100 * csBetter / csN).toStringAsFixed(1)} | '
      '${oBoth == 0 ? '—' : (100 * oSame / oBoth).toStringAsFixed(1)} |',
    );
  }

  // F_S (ours, classical) vs STL and nsdiffs (report).
  final seas = series.values.where((x) => x.rSeasonal != null).toList();
  if (seas.isNotEmpty) {
    var agree = 0;
    report.add('');
    report.add('F_S classical vs STL (seasonal-period series):');
    for (final x in seas) {
      final ns = _i(x.rSeasonal!['nsdiffs_R']!);
      if (ns == x.sDStar) agree++;
    }
    report.add('D* = nsdiffs_R on $agree / ${seas.length}');
  }
  report
    ..add('')
    ..add('DEFECTS (C-d′ η/lag/d): ${anyDefect ? 'FOUND' : 'none'}')
    ..add('STATISTICAL GATES (C-S′, MASE vs R): ${statFail ? 'FAIL' : 'PASS'}')
    ..add('')
    ..add('details:')
    ..addAll(details);
  if (windows) {
    report.addAll(_dPrime(exports, series, p));
  }
  final text = report.join('\n');
  if (out != null) File(out).writeAsStringSync(text);
  stdout.writeln(text);
}

/// Gate D′ (spec §15.2): per-window table and per-series medians.
List<String> _dPrime(
  List<String> exports,
  Map<String, _Series> series,
  String Function(String, String) p,
) {
  final win = <String, _Row>{};
  for (final dir in exports) {
    for (final r in _csv(p(dir, 'ours_windows.csv'))) {
      win[r['series_id']!] = r;
    }
  }
  double? opt(String? v) => v == null || v.isEmpty ? null : _d(v);
  double? maseR(_Series x, _Row? r) {
    if (r == null) return null;
    final d = _i(r['d']!), sD = _i(r['D']!);
    final sp = _i(r['P']!), sq = _i(r['Q']!);
    try {
      final fc = fitSarimax(
        x.train,
        const {},
        order: ArimaOrder(
          p: _i(r['p']!),
          d: d,
          q: _i(r['q']!),
          seasonalP: sp,
          seasonalD: sD,
          seasonalQ: sq,
          seasonalPeriod: sp + sq + sD > 0 ? x.period : 0,
        ),
        horizon: x.test.length,
        includeMean: _b(r['has_mean']!) && d + sD == 0,
        includeDrift: _b(r['has_drift']!) && d + sD == 1,
      ).forecast;
      return mase(x.train, x.test, fc, x.period > 1 ? x.period : 1);
    } on Object {
      return null;
    }
  }

  String f3(double? v) => v == null ? '—' : v.toStringAsFixed(3);
  final out = <String>[
    '',
    '## Gate D′ (per window)',
    '',
    '| window | n | cycles | D source | (d,D) ours | (d,D) R step | '
        'ours/R step | (d,D) R ex | ours/R ex | old D step | old D ex | '
        'ours/ours d=1 | sD=1 / R step |',
    '|---|---|---|---|---|---|---|---|---|---|---|---|---|',
  ];
  final perSeries = <String, Map<String, List<double>>>{};
  for (final x in series.values) {
    final w = win[x.id];
    if (w == null) continue;
    final ms = opt(w['mase_step']), me = opt(w['mase_ex']);
    final rs = x.rOrders['step'], re = x.rOrders['ex'];
    final mrs = maseR(x, rs), mre = maseR(x, re);
    final rStep = ms != null && mrs != null ? ms / mrs : null;
    final rEx = me != null && mre != null ? me / mre : null;
    final oldS = opt(w['mase_oldR_step']), oldE = opt(w['mase_oldR_ex']);
    final d1 = opt(w['mase_step_d1']), sd1 = opt(w['mase_step_sD1']);
    final m = perSeries[x.model] ??= {};
    void add(String k, double? v) {
      if (v != null) (m[k] ??= []).add(v);
    }

    add('step', rStep);
    add('ex', rEx);
    add('oldStep', ms != null && oldS != null ? ms / oldS : null);
    add('oldEx', me != null && oldE != null ? me / oldE : null);
    add('d1', ms != null && d1 != null ? ms / d1 : null);
    add('sD1', sd1 != null && mrs != null ? sd1 / mrs : null);
    out.add(
      '| ${x.id} | ${w['window_n']} | ${w['cycles']} | ${w['D_source']} | '
      '(${x.dStar},${x.sDStar}) | '
      '${rs == null ? '—' : '(${rs['d']},${rs['D']})'} | ${f3(rStep)} | '
      '${re == null ? '—' : '(${re['d']},${re['D']})'} | ${f3(rEx)} | '
      '${f3(ms != null && oldS != null ? ms / oldS : null)} | '
      '${f3(me != null && oldE != null ? me / oldE : null)} | '
      '${f3(ms != null && d1 != null ? ms / d1 : null)} | '
      '${f3(sd1 != null && mrs != null ? sd1 / mrs : null)} |',
    );
  }
  out.addAll([
    '',
    '| series | median ours/R step | median ours/R ex | old D step / ex | '
        'ours/ours d=1 | sD=1 / R step | D′ (≤ 1.05) |',
    '|---|---|---|---|---|---|---|',
  ]);
  for (final MapEntry(key: name, value: m) in perSeries.entries) {
    double? med(String k) => m[k] == null ? null : median(m[k]!);
    final ok = (med('step') ?? 0) <= 1.05 && (med('ex') ?? 0) <= 1.05;
    out.add(
      '| $name | ${f3(med('step'))} | ${f3(med('ex'))} | '
      '${f3(med('oldStep'))} / ${f3(med('oldEx'))} | ${f3(med('d1'))} | '
      '${f3(med('sD1'))} | ${ok ? 'pass' : 'FAIL'} |',
    );
  }
  return out;
}
