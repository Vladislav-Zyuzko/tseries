/// The synthetic study of the autoArima specification (§9.C): data-generating
/// processes, seeded series, the statistics of the study (Wilson intervals,
/// bootstrap of a median, MASE) and the CSV exchange format with the external
/// R run. Shared by test/auto_arima_synthetic_test.dart and
/// tool/auto_arima_r_compare.dart.
library;

import 'dart:io';
import 'dart:math' as math;

import 'arima_simulator.dart';

/// One data-generating process of the C study.
class SyntheticModel {
  const SyntheticModel(
    this.id,
    this.label, {
    this.ar = const [],
    this.ma = const [],
    this.sma = const [],
    this.d = 0,
    this.sD = 0,
    this.period = 1,
    this.drift = 0,
    this.mean = 0,
    this.nominalGate = false,
    this.persistent = false,
    this.nearUnit = false,
  });

  /// Short id used in series ids (`<id>_n<n>_r<rep>`).
  final String id;
  final String label;
  final List<double> ar;

  /// MA in the `1 + θB` convention.
  final List<double> ma;
  final List<double> sma;
  final int d;
  final int sD;
  final int period;
  final double drift;
  final double mean;

  /// C-d nominal gate applies (white noise, MA(1), (0,1,1), (0,1,1)+drift).
  final bool nominalGate;

  /// C-d′ applies instead (persistent models).
  final bool persistent;

  /// φ = 0.98: report only.
  final bool nearUnit;

  /// Holdout length: 2s for a seasonal model, 10 otherwise.
  int get holdout => period > 1 ? 2 * period : 10;

  /// Whether the true model has a constant at differencing (dd, sD): a mean
  /// when nothing is differenced (every DGP has a non-zero level then), a
  /// drift when d + D = 1 and the DGP has one.
  bool constantAt(int dd, int sDD) =>
      dd + sDD == 0 ? true : (dd + sDD == 1 ? drift != 0 : false);
}

/// The 11 models of §9.C.
const syntheticModels = <SyntheticModel>[
  SyntheticModel('wn', 'white noise', mean: 10, nominalGate: true),
  SyntheticModel('ar1', 'AR(1) φ=0.7', ar: [0.7], mean: 10, persistent: true),
  SyntheticModel('ma1', 'MA(1) θ=0.5', ma: [0.5], mean: 10, nominalGate: true),
  SyntheticModel(
    'arma11',
    'ARMA(1,1) 0.5/0.4',
    ar: [0.5],
    ma: [0.4],
    mean: 10,
    persistent: true,
  ),
  SyntheticModel(
    'ar2',
    'AR(2) 0.5/0.3',
    ar: [0.5, 0.3],
    mean: 10,
    persistent: true,
  ),
  SyntheticModel(
    'ima',
    'ARIMA(0,1,1) θ=0.4',
    ma: [0.4],
    d: 1,
    nominalGate: true,
  ),
  SyntheticModel(
    'ari',
    'ARIMA(1,1,0) φ=0.5',
    ar: [0.5],
    d: 1,
    persistent: true,
  ),
  SyntheticModel(
    'arima111',
    'ARIMA(1,1,1) 0.5/0.3',
    ar: [0.5],
    ma: [0.3],
    d: 1,
    persistent: true,
  ),
  SyntheticModel(
    'imadrift',
    'ARIMA(0,1,1)+drift 0.3',
    ma: [0.4],
    d: 1,
    drift: 0.3,
    nominalGate: true,
  ),
  SyntheticModel(
    'airline',
    '(0,1,1)(0,1,1)[12] −0.4/−0.6',
    ma: [-0.4],
    sma: [-0.6],
    d: 1,
    sD: 1,
    period: 12,
  ),
  SyntheticModel('ar098', 'AR(1) φ=0.98', ar: [0.98], mean: 10, nearUnit: true),
];

const syntheticLengths = [200, 500];

/// Seed of a series (fixed formula; identical to the first C run).
int syntheticSeed(SyntheticModel m, int n, int rep) =>
    1000003 * (syntheticModels.indexOf(m) + 1) + 7919 * n + rep;

String syntheticId(SyntheticModel m, int n, int rep) => '${m.id}_n${n}_r$rep';

/// Training part (n values) and holdout of one series.
({List<double> train, List<double> test}) syntheticSeries(
  SyntheticModel m,
  int n,
  int rep,
) {
  final h = m.holdout;
  final y = simulateArima(
    n: n + h,
    seed: syntheticSeed(m, n, rep),
    ar: m.ar,
    ma: m.ma,
    seasonalMa: m.sma,
    d: m.d,
    seasonalD: m.sD,
    period: m.period,
    constant: m.drift,
    mean: m.mean,
  );
  return (train: y.sublist(0, n), test: y.sublist(n));
}

/// MASE of [fc] on [test]; scale = in-sample mean absolute (seasonal) naive
/// error of [train] at lag [m].
double mase(List<double> train, List<double> test, List<double> fc, int m) {
  var scale = 0.0;
  for (var t = m; t < train.length; t++) {
    scale += (train[t] - train[t - m]).abs();
  }
  scale /= train.length - m;
  var err = 0.0;
  for (var i = 0; i < test.length; i++) {
    err += (test[i] - fc[i]).abs();
  }
  return err / test.length / scale;
}

/// 95 % Wilson score interval of k successes in n trials.
({double lo, double hi}) wilson(int k, int n) {
  if (n == 0) return (lo: double.nan, hi: double.nan);
  const z = 1.959963984540054;
  final p = k / n;
  final den = 1 + z * z / n;
  final centre = (p + z * z / (2 * n)) / den;
  final half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den;
  return (lo: math.max(0, centre - half), hi: math.min(1, centre + half));
}

String pctCi(int k, int n) {
  if (n == 0) return '—';
  final w = wilson(k, n);
  return '${(100 * k / n).toStringAsFixed(1)} '
      '[${(100 * w.lo).toStringAsFixed(1)}, ${(100 * w.hi).toStringAsFixed(1)}]';
}

double median(List<double> xs) {
  if (xs.isEmpty) return double.nan;
  final s = [...xs]..sort();
  final m = s.length ~/ 2;
  return s.length.isOdd ? s[m] : 0.5 * (s[m - 1] + s[m]);
}

/// Median of [xs] with a 95 % percentile-bootstrap interval (2000
/// resamples, fixed seed).
({double median, double lo, double hi}) bootstrapMedian(List<double> xs) {
  if (xs.isEmpty) return (median: double.nan, lo: double.nan, hi: double.nan);
  final rng = math.Random(20261006);
  final meds = <double>[];
  for (var b = 0; b < 2000; b++) {
    final sample = [
      for (var i = 0; i < xs.length; i++) xs[rng.nextInt(xs.length)],
    ];
    meds.add(median(sample));
  }
  meds.sort();
  return (
    median: median(xs),
    lo: meds[(0.025 * (meds.length - 1)).round()],
    hi: meds[(0.975 * (meds.length - 1)).round()],
  );
}

/// 17 significant digits: round-trips every double.
String g17(double v) => v.toStringAsPrecision(17);

/// Appends lines to a file (creating its directory).
class LineSink {
  LineSink(String path, {String? header})
    : _f = (File(
        path,
      )..parent.createSync(recursive: true)).openSync(mode: FileMode.write) {
    if (header != null) writeln(header);
  }

  final RandomAccessFile _f;
  final _buf = StringBuffer();

  void writeln(String line) {
    _buf.writeln(line);
    if (_buf.length > 1 << 20) flush();
  }

  void flush() {
    _f.writeStringSync(_buf.toString());
    _buf.clear();
  }

  void close() {
    flush();
    _f.closeSync();
  }
}
