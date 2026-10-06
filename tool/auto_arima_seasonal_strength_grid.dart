// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// "C-D" of the autoArima specification (§9.C, revision 5): the strength of
// seasonality F_S (classical decomposition) on seeded series, for the
// calibration of the seasonal-difference threshold τ_S. Only F_S is
// computed — no fits.
//
// Grid: s ∈ {4, 7, 12, 24} × k ∈ {3, 5, 10, 20} cycles (N = k·s), 100
// repetitions per DGP:
//  * null family (D = 0 is right): white noise, AR(1) 0.7, AR(1) 0.98,
//    AR(2) 0.5/0.3, ARIMA(0,1,1) θ=0.4, (0,1,1)+drift 0.3, (1,1,0) φ=0.5;
//  * seasonal family (D = 1 is right): (0,1,1)(0,1,1)_s, θ = −0.4,
//    Θ ∈ {−0.3, −0.6, −0.9}; (1,0,0)(0,1,1)_s, φ = 0.5, Θ = −0.6;
//  * deterministic seasonality over AR(1) 0.5: amplitude·(sin(2πt/s) +
//    ½·sin(4πt/s)), amplitude ∈ {0.5, 1, 2} (σ = 1) — no "right" D, report.
//
// Usage (from packages/tseries):
//   dart run tool/auto_arima_seasonal_strength_grid.dart --values <file>
//       [--out <table.md>] [--reps 100] [--cycles 3,5,10,20]
//
// Each cell is judged at the threshold of the stepped rule of spec §15.1
// for its k (0.70 for k < 10 — also reported for k < 5, where autoArima
// takes no automatic D at all —, 0.35 for k ≥ 10). Gate (§15.1): false
// D = 1 ≤ 2 % (point) in every (s, null DGP) cell with k ≥ 5; power on
// Θ ∈ {−0.3, −0.6} with k ≥ 10: Wilson upper bound ≥ 90 %.
// --values: one line per series `dgp|family|s|k|N|rep|seed|F_S`.
import 'dart:io';
import 'dart:math' as math;

import 'package:tseries/tseries.dart';

import '../test/support/arima_simulator.dart';
import '../test/support/auto_arima_synthetic.dart';

typedef _Dgp = ({
  String id,
  String family,
  List<double> Function(int n, int s, int seed) make,
});

List<double> _sim(
  int n,
  int seed, {
  List<double> ar = const [],
  List<double> ma = const [],
  List<double> sma = const [],
  int d = 0,
  int sD = 0,
  int s = 1,
  double drift = 0,
  double mean = 0,
}) => simulateArima(
  n: n,
  seed: seed,
  ar: ar,
  ma: ma,
  seasonalMa: sma,
  d: d,
  seasonalD: sD,
  period: s,
  constant: drift,
  mean: mean,
);

final List<_Dgp> _dgps = [
  (id: 'wn', family: 'null', make: (n, s, seed) => _sim(n, seed, mean: 10)),
  (
    id: 'ar1_0.7',
    family: 'null',
    make: (n, s, seed) => _sim(n, seed, ar: [0.7], mean: 10),
  ),
  (
    id: 'ar1_0.98',
    family: 'null',
    make: (n, s, seed) => _sim(n, seed, ar: [0.98], mean: 10),
  ),
  (
    id: 'ar2',
    family: 'null',
    make: (n, s, seed) => _sim(n, seed, ar: [0.5, 0.3], mean: 10),
  ),
  (
    id: 'ima',
    family: 'null',
    make: (n, s, seed) => _sim(n, seed, ma: [0.4], d: 1),
  ),
  (
    id: 'ima_drift',
    family: 'null',
    make: (n, s, seed) => _sim(n, seed, ma: [0.4], d: 1, drift: 0.3),
  ),
  (
    id: 'ari',
    family: 'null',
    make: (n, s, seed) => _sim(n, seed, ar: [0.5], d: 1),
  ),
  for (final big in [-0.3, -0.6, -0.9])
    (
      id: 'airline_Θ$big',
      family: 'seasonal',
      make: (n, s, seed) =>
          _sim(n, seed, ma: [-0.4], sma: [big], d: 1, sD: 1, s: s),
    ),
  (
    id: 'ar1_sma',
    family: 'seasonal',
    make: (n, s, seed) => _sim(n, seed, ar: [0.5], sma: [-0.6], sD: 1, s: s),
  ),
  for (final amp in [0.5, 1.0, 2.0])
    (
      id: 'determ_$amp',
      family: 'deterministic',
      make: (n, s, seed) {
        final w = _sim(n, seed, ar: [0.5]);
        return [
          for (var t = 0; t < n; t++)
            w[t] +
                amp *
                    (math.sin(2 * math.pi * t / s) +
                        0.5 * math.sin(4 * math.pi * t / s)),
        ];
      },
    ),
];

void main(List<String> args) {
  String? arg(String name) {
    final i = args.indexOf(name);
    return i >= 0 && i + 1 < args.length ? args[i + 1] : null;
  }

  final reps = int.tryParse(arg('--reps') ?? '') ?? 100;
  final cycles = (arg('--cycles') ?? '3,5,10,20')
      .split(',')
      .map(int.parse)
      .toList();
  final values = LineSink(
    arg('--values') ?? 'fs_values.txt',
    header: '# dgp|family|s|k|N|rep|seed|F_S',
  );
  // (s, k) -> dgp -> F_S list
  final byCell = <(int, int), Map<String, List<double>>>{};
  for (final s in [4, 7, 12, 24]) {
    for (final k in cycles) {
      final n = s * k;
      final cell = byCell[(s, k)] = {};
      for (var g = 0; g < _dgps.length; g++) {
        final dgp = _dgps[g];
        final list = cell[dgp.id] = [];
        for (var rep = 0; rep < reps; rep++) {
          final seed = 7000003 * (g + 1) + 10007 * s + 101 * k + rep;
          final y = dgp.make(n, s, seed);
          final f = seasonalStrength(y, period: s);
          list.add(f);
          values.writeln(
            '${dgp.id}|${dgp.family}|$s|$k|$n|$rep|$seed|${g17(f)}',
          );
        }
      }
    }
  }
  values.close();

  double tauOf(int k) => k >= autoArimaSeasonalManyCycles
      ? autoArimaSeasonalStrengthThreshold
      : autoArimaSeasonalStrengthThresholdFewCycles;
  String fam(String id) => _dgps.firstWhere((d) => d.id == id).family;
  final rows = <String>[
    '### At the §15.1 threshold of each k: false D = 1 on the null family '
        '(worst DGP) and power on the seasonal family, % [Wilson 95 %]',
    '',
    '| s | k | N | τ | worst null DGP | false D % | null max F_S (DGP) | '
        'airline Θ=−0.3 | Θ=−0.6 | Θ=−0.9 | (1,0,0)(0,1,1) | '
        'seasonal min F_S (DGP) | determ 0.5 / 1 / 2 med F_S |',
    '|---|---|---|---|---|---|---|---|---|---|---|---|---|',
  ];
  final notes = <String>[];
  for (final MapEntry(key: (s, k), value: cell) in byCell.entries) {
    var worst = '';
    var worstCount = -1;
    var nullMax = -1.0;
    var nullMaxId = '';
    var seasMin = 2.0;
    var seasMinId = '';
    final power = <String, String>{};
    final tau = tauOf(k);
    for (final MapEntry(key: id, value: fs) in cell.entries) {
      final above = fs.where((v) => v >= tau).length;
      switch (fam(id)) {
        case 'null':
          if (above > worstCount) {
            worstCount = above;
            worst = id;
          }
          for (final v in fs) {
            if (v > nullMax) {
              nullMax = v;
              nullMaxId = id;
            }
          }
          if (k >= autoArimaSeasonalMinCycles && above > 0.02 * fs.length) {
            notes.add('false D: s=$s k=$k $id $above/${fs.length}');
          }
        case 'seasonal':
          power[id] = pctCi(above, fs.length);
          for (final v in fs) {
            if (v < seasMin) {
              seasMin = v;
              seasMinId = id;
            }
          }
          if (k >= autoArimaSeasonalManyCycles &&
              (id == 'airline_Θ-0.3' || id == 'airline_Θ-0.6') &&
              wilson(above, fs.length).hi < 0.9) {
            notes.add('power: s=$s k=$k $id $above/${fs.length}');
          }
        default:
      }
    }
    String med(String id) => median(cell[id]!).toStringAsFixed(3);
    rows.add(
      '| $s | $k | ${s * k} | ${tauOf(k)} | $worst | ${pctCi(worstCount, reps)} | '
      '${nullMax.toStringAsFixed(3)} ($nullMaxId) | '
      '${power['airline_Θ-0.3']} | ${power['airline_Θ-0.6']} | '
      '${power['airline_Θ-0.9']} | ${power['ar1_sma']} | '
      '${seasMin.toStringAsFixed(3)} ($seasMinId) | '
      '${med('determ_0.5')} / ${med('determ_1.0')} / ${med('determ_2.0')} |',
    );
  }

  // Threshold sweep: worst false-D share (any null DGP, any cell with
  // k ≥ kMin) and worst power (Θ ∈ {−0.3, −0.6}, k ≥ 5) per candidate τ.
  rows.addAll([
    '',
    '### Threshold sweep (cells with k ≥ kMin for false D; power over k ≥ 5)',
    '',
    '| τ_S | worst false D %, k≥3 | k≥4 (k≥5) | worst power Θ∈{−0.3,−0.6}, k≥5 |',
    '|---|---|---|---|',
  ]);
  for (final t in [0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.64]) {
    double worstFalse(int kMin) {
      var w = 0.0;
      for (final MapEntry(key: (_, k), value: cell) in byCell.entries) {
        if (k < kMin) continue;
        for (final MapEntry(key: id, value: fs) in cell.entries) {
          if (fam(id) != 'null') continue;
          w = math.max(w, fs.where((v) => v >= t).length / fs.length);
        }
      }
      return w;
    }

    var worstPower = 1.0;
    for (final MapEntry(key: (_, k), value: cell) in byCell.entries) {
      if (k < 5) continue;
      for (final id in ['airline_Θ-0.3', 'airline_Θ-0.6']) {
        final fs = cell[id]!;
        worstPower = math.min(
          worstPower,
          fs.where((v) => v >= t).length / fs.length,
        );
      }
    }
    rows.add(
      '| $t | ${(100 * worstFalse(3)).toStringAsFixed(1)} | '
      '${(100 * worstFalse(5)).toStringAsFixed(1)} | '
      '${(100 * worstPower).toStringAsFixed(1)} |',
    );
  }
  rows.add('');
  rows.add(
    'gate notes (§15.1 rule): ${notes.isEmpty ? 'none' : notes.join('; ')}',
  );
  final text = rows.join('\n');
  final out = arg('--out');
  if (out != null) File(out).writeAsStringSync(text);
  stdout.writeln(text);
}
