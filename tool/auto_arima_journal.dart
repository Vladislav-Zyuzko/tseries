// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// Runs `autoArima` on the six reference series of
// test/fixtures/auto_arima/ in both search modes and writes:
//
//  * --journal <file>: the full journals (differencing report + every
//    candidate record), as stable text. Two fresh processes of the same build
//    must write byte-identical files (gate E of the autoArima spec, §9).
//  * --g0 <file>: one line per fitted candidate (no cache hits) with the
//    order, constant, estimated parameters and logLik — the export for the
//    external likelihood check against statsmodels (gate G0, spec §9). The
//    gate itself is run by someone other than the search implementer.
//
// Usage (from packages/tseries):
//   dart run tool/auto_arima_journal.dart --journal out.txt [--g0 g0.txt]
//       [--reverse] [--timing]
//
// --reverse processes the series in reverse order (different heap history).
// --timing prints fits and wall time per series and mode to stdout.
// --boundary also appends the journals of short seasonal series around the
// length rule of spec §5 (s in {4, 12}, D in {0, 1}; gate F/E).
import 'dart:io';

import 'package:tseries/tseries.dart';

import '../test/support/arima_simulator.dart';

const _series = <(String, int)>[
  ('air_passengers_log', 12),
  ('wwwusage', 1),
  ('lynx', 1),
  ('lake_huron', 1),
  ('nile', 1),
  ('us_acc_deaths', 12),
];

/// Budget large enough that neither mode is truncated on these series.
const _maxModels = 1000;

List<double> _load(String name) =>
    File('test/fixtures/auto_arima/series_$name.csv')
        .readAsLinesSync()
        .where((l) => l.trim().isNotEmpty)
        .map(double.parse)
        .toList();

String _list(List<double> xs) => xs.join(',');

void main(List<String> args) {
  String? journalPath;
  String? g0Path;
  var reverse = false;
  var timing = false;
  var boundary = false;
  for (var i = 0; i < args.length; i++) {
    switch (args[i]) {
      case '--journal':
        journalPath = args[++i];
      case '--g0':
        g0Path = args[++i];
      case '--reverse':
        reverse = true;
      case '--timing':
        timing = true;
      case '--boundary':
        boundary = true;
      default:
        stderr.writeln('unknown argument ${args[i]}');
        exit(64);
    }
  }

  final journals = <String, String>{};
  final g0 = StringBuffer()
    ..writeln(
      '# autoArima candidate fits for the external logLik check (gate G0).',
    )
    ..writeln(
      '# MA/SMA coefficients are in ctsa\'s 1 - theta*B convention '
      '(opposite sign to statsmodels). c = mean when d+D = 0, drift when '
      'd+D = 1. sigma2 = profile MLE S/T.',
    )
    ..writeln(
      '# series|mode|index|p|d|q|P|D|Q|s|c|verdict|status|loglik|sigma2|'
      'ar|ma|sar|sma|intercept|drift',
    );

  final order = reverse ? _series.reversed.toList() : _series;
  for (final (name, period) in order) {
    final y = _load(name);
    for (final mode in ArimaSearch.values) {
      final sw = Stopwatch()..start();
      String text;
      List<CandidateRecord> trace;
      int fits;
      try {
        final r = autoArima(
          y,
          period: period,
          search: mode,
          maxModels: _maxModels,
        );
        text = r.toJournal();
        trace = r.trace;
        fits = r.fitsEvaluated;
      } on AutoArimaNoModelException catch (e) {
        text =
            '${e.differencing.toJournal()}'
            '${e.trace.map((r) => r.toJournalLine()).join('\n')}\n'
            'NO MODEL: ${e.message}\n';
        trace = e.trace;
        fits = e.fitsEvaluated;
      }
      sw.stop();
      if (timing) {
        stdout.writeln(
          '$name ${mode.name}: fits=$fits time_ms=${sw.elapsedMilliseconds}',
        );
      }
      journals['$name|${mode.name}'] = text;
      for (final r in trace) {
        final f = r.fit;
        if (r.duplicateOf != null || f == null) continue;
        final o = r.order;
        final drift = [
          for (final x in f.regressors)
            if (x.name == sarimaxDriftColumn) x.coefficient,
        ];
        g0.writeln(
          [
            name,
            mode.name,
            r.index,
            o.p,
            o.d,
            o.q,
            o.seasonalP,
            o.seasonalD,
            o.seasonalQ,
            o.seasonalPeriod,
            r.constant ? 1 : 0,
            r.verdict.name,
            r.ctsaStatus ?? '-',
            r.logLikelihood ?? '-',
            f.sigma2,
            _list(f.ar),
            _list(f.ma),
            _list(f.seasonalAr),
            _list(f.seasonalMa),
            f.intercept ?? '-',
            drift.isEmpty ? '-' : drift.single,
          ].join('|'),
        );
      }
    }
  }

  // Always write in the canonical series order, whatever the run order.
  final out = StringBuffer();
  for (final (name, _) in _series) {
    for (final mode in ArimaSearch.values) {
      out
        ..writeln('== $name ${mode.name}')
        ..write(journals['$name|${mode.name}']);
    }
  }
  if (boundary) {
    for (final s in [4, 12]) {
      for (final sD in [0, 1]) {
        final lengths = s == 12
            ? [21, 22, 23, 26, 27, 28]
            : [9, 10, 11, 14, 15];
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
            out.writeln('== boundary s=$s D=$sD n=$n ${mode.name}');
            try {
              out.write(
                autoArima(
                  y,
                  period: s,
                  seasonalD: sD,
                  search: mode,
                  maxP: 2,
                  maxQ: 2,
                ).toJournal(),
              );
            } on AutoArimaNoModelException catch (e) {
              out
                ..write(e.differencing.toJournal())
                ..writeln(e.trace.map((r) => r.toJournalLine()).join('\n'))
                ..writeln('NO MODEL: ${e.message}');
            }
          }
        }
      }
    }
  }
  if (journalPath != null) File(journalPath).writeAsStringSync(out.toString());
  if (g0Path != null) File(g0Path).writeAsStringSync(g0.toString());
  if (journalPath == null && g0Path == null && !timing) stdout.write(out);
}
