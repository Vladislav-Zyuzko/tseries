/// Gate A of the autoArima specification (§9): the pure-Dart components
/// against independent references. statsmodels/numpy were run as external
/// tools by tool/auto_arima_components_reference.py; only their numbers are
/// in test/fixtures/auto_arima_components/.
library;

import 'dart:io';
import 'dart:math' as math;

import 'package:test/test.dart';
import 'package:tseries/tseries.dart';

const _dir = 'test/fixtures/auto_arima_components';

List<double> _series(String name) =>
    File('test/fixtures/auto_arima/series_$name.csv')
        .readAsLinesSync()
        .where((l) => l.trim().isNotEmpty)
        .map(double.parse)
        .toList();

List<List<String>> _rows(String file) => File('$_dir/$file')
    .readAsLinesSync()
    .where((l) => l.isNotEmpty && !l.startsWith('#'))
    .map((l) => l.split('|'))
    .toList();

List<double> _diff(List<double> x) => [
  for (var i = 1; i < x.length; i++) x[i] - x[i - 1],
];

List<double> _sdiff(List<double> x, int s) => [
  for (var i = s; i < x.length; i++) x[i] - x[i - s],
];

double _rel(double a, double b) => (a - b).abs() / math.max(b.abs(), 1e-300);

void main() {
  group('KPSS η against statsmodels kpss (same lag)', () {
    final rows = _rows('kpss_reference.txt');
    test('fixture has the 6 series and their differences', () {
      expect(rows.length, 22);
    });
    for (final r in rows) {
      final [name, transform, tText, lagText, etaText] = r;
      test('$name $transform', () {
        final y = _series(name);
        final x = switch (transform) {
          'level' => y,
          'diff1' => _diff(y),
          'diff2' => _diff(_diff(y)),
          'sdiff' => _sdiff(y, 12),
          'sdiff_diff1' => _diff(_sdiff(y, 12)),
          _ => throw StateError(transform),
        };
        expect(x.length, int.parse(tText));
        // The l4 rule gives the lag statsmodels was called with.
        expect(KpssLag.l4.lagFor(x.length), int.parse(lagText));
        final k = kpssLevelTest(x, lag: KpssLag.l4);
        expect(k.lag, int.parse(lagText));
        expect(_rel(k.statistic, double.parse(etaText)), lessThan(1e-8));
      });
    }
  });

  group(
    'KPSS η against R ur.kpss(type = "mu") (external tool, numbers only)',
    () {
      // out2.txt: η at use.lag 0..8, printed with 6 decimals; out.txt: η at
      // R's "short" lag with 12 decimals and the lag R used.
      final out2 = File('test/fixtures/auto_arima/r_probe/out2.txt');
      final out1 = File('test/fixtures/auto_arima/r_probe/out.txt');
      if (!out2.existsSync() || !out1.existsSync()) {
        // The R probe output is development material and is not shipped in
        // the published package (.pubignore).
        test(
          'R probe output',
          () {},
          skip: 'test/fixtures/auto_arima/r_probe/ is not in this checkout',
        );
        return;
      }
      final y = {
        for (final (n, _) in [
          ('air_passengers_log', 12),
          ('wwwusage', 1),
          ('lake_huron', 1),
          ('nile', 1),
          ('lynx', 1),
          ('us_acc_deaths', 12),
        ])
          n: _series(n),
      };
      final probes = <String, List<double>>{
        'air_log': y['air_passengers_log']!,
        'air_log_sdiff12': _sdiff(y['air_passengers_log']!, 12),
        'air_log_diff1': _diff(y['air_passengers_log']!),
        'wwwusage': y['wwwusage']!,
        'wwwusage_diff1': _diff(y['wwwusage']!),
        'lake_huron': y['lake_huron']!,
        'nile': y['nile']!,
        'lynx': y['lynx']!,
        'us_acc_deaths': y['us_acc_deaths']!,
        'us_acc_sdiff12': _sdiff(y['us_acc_deaths']!, 12),
      };
      final grid = out2.readAsLinesSync().where(
        (l) => l.contains('teststat by use.lag'),
      );
      for (final line in grid) {
        final name = line.split(' ').first;
        test('$name, use.lag 0..8', () {
          final x = probes[name]!;
          final values = line
              .split(':')
              .last
              .trim()
              .split(RegExp(r'\s+'))
              .map(double.parse)
              .toList();
          expect(values, hasLength(9));
          for (var l = 0; l <= 8; l++) {
            final ours = kpssLevelTest(x, lag: KpssLag.fixed(l)).statistic;
            // 6 printed decimals: half a unit in the last place, or 1e-6 rel.
            expect(
              (ours - values[l]).abs(),
              lessThanOrEqualTo(math.max(1e-6 * values[l], 5.000001e-7)),
              reason: 'lag $l: ours $ours, R ${values[l]}',
            );
          }
        });
      }
      final short = out1.readAsLinesSync().where((l) => l.startsWith('KPSS '));
      for (final line in short) {
        final name = line.split(' ')[1].replaceAll(':', '');
        test('$name, R "short" lag (12 decimals)', () {
          final x = probes[name]!;
          final eta = double.parse(
            RegExp(r'teststat=([0-9.]+)').firstMatch(line)!.group(1)!,
          );
          final lag = int.parse(
            RegExp(r' lag=(\d+)').firstMatch(line)!.group(1)!,
          );
          // urca's own lags = "short" option (unrelated to our KpssLag.short):
          // on these lengths it reported the l4 lags (observed in the output).
          // What is checked: the same lag gives the same η.
          expect(KpssLag.l4.lagFor(x.length), lag);
          final ours = kpssLevelTest(x, lag: KpssLag.fixed(lag)).statistic;
          expect(_rel(ours, eta), lessThan(1e-9), reason: 'ours $ours R $eta');
        });
      }
    },
  );

  group('KPSS details', () {
    test('Table 1 critical values and p-value brackets', () {
      expect(kpssLevelCriticalValues, {
        0.10: 0.347,
        0.05: 0.463,
        0.025: 0.574,
        0.01: 0.739,
      });
      expect(KpssPValueBracket.of(0.80), KpssPValueBracket.below001);
      expect(KpssPValueBracket.of(0.739), KpssPValueBracket.below001);
      expect(KpssPValueBracket.of(0.6), KpssPValueBracket.from001To0025);
      expect(KpssPValueBracket.of(0.5), KpssPValueBracket.from0025To005);
      expect(KpssPValueBracket.of(0.4), KpssPValueBracket.from005To010);
      expect(KpssPValueBracket.of(0.1), KpssPValueBracket.above010);
    });

    test('decision is η > crit(α)', () {
      final x = _series('wwwusage'); // η(l4 = 4) = 0.4542
      expect(kpssLevelTest(x, lag: KpssLag.l4).rejectsStationarity, isFalse);
      expect(
        kpssLevelTest(x, lag: KpssLag.l4, alpha: 0.10).rejectsStationarity,
        isTrue,
      );
      // Default short rule: l = 2 at T = 100, η = 0.722 > 0.463.
      final k = kpssLevelTest(x);
      expect(k.lag, 2);
      expect(k.lagRule, KpssLag.short);
      expect(k.rejectsStationarity, isTrue);
    });

    test('lag rules (exact integer boundaries)', () {
      expect(KpssLag.l4.lagFor(100), 4);
      expect(KpssLag.l4.lagFor(99), 3);
      expect(KpssLag.l4.lagFor(12), 2);
      expect(KpssLag.l4.lagFor(1000), 7);
      expect(KpssLag.l12.lagFor(100), 12);
      expect(KpssLag.l12.lagFor(132), 12);
      // short = floor(3·sqrt(T)/13): 0 below 19, boundaries 75.1 and 169.
      expect(KpssLag.short.lagFor(18), 0);
      expect(KpssLag.short.lagFor(19), 1);
      expect(KpssLag.short.lagFor(75), 1);
      expect(KpssLag.short.lagFor(76), 2);
      expect(KpssLag.short.lagFor(168), 2);
      expect(KpssLag.short.lagFor(169), 3);
      expect(KpssLag.short.lagFor(100), 2);
      expect(KpssLag.short.lagFor(132), 2);
      expect(const KpssLag.fixed(5).lagFor(1000), 5);
      expect(const KpssLag.fixed(5), const KpssLag.fixed(5));
    });

    test('refuses a constant series, bad α, bad lag, NaN', () {
      expect(() => kpssLevelTest(List.filled(20, 3.0)), throwsArgumentError);
      expect(
        () => kpssLevelTest(_series('nile'), alpha: 0.2),
        throwsArgumentError,
      );
      expect(
        () => kpssLevelTest(_series('nile'), lag: const KpssLag.fixed(100)),
        throwsArgumentError,
      );
      expect(() => kpssLevelTest([1, double.nan, 2]), throwsArgumentError);
    });
  });

  group('classical decomposition and F_S against statsmodels', () {
    for (final r in _rows('decomposition_reference.txt')) {
      final [name, periodText, fsText, indicesText] = r;
      final s = int.parse(periodText);
      test('$name, s = $s', () {
        final dec = classicalDecomposition(_series(name), period: s);
        expect(dec.isDegenerate, isFalse);
        expect(
          _rel(dec.seasonalStrength, double.parse(fsText)),
          lessThan(1e-10),
        );
        final ref = indicesText.split(',').map(double.parse).toList();
        expect(dec.seasonalIndices.length, s);
        final scale = ref.map((v) => v.abs()).reduce(math.max);
        for (var j = 0; j < s; j++) {
          expect(
            (dec.seasonalIndices[j] - ref[j]).abs(),
            lessThan(1e-10 * scale),
          );
        }
      });
    }

    test('hand calculation on a tiny series (s = 4)', () {
      // y = 0..11 plus pattern (1, -1, 2, -2): the 2x4 MA of the pattern is
      // 0, so T_t = t on t = 2..9 and y − T is the pattern exactly.
      final pattern = [1.0, -1.0, 2.0, -2.0];
      final y = [for (var t = 0; t < 12; t++) t + pattern[t % 4]];
      final dec = classicalDecomposition(y, period: 4);
      for (var t = 2; t < 10; t++) {
        expect(dec.trend[t], closeTo(t.toDouble(), 1e-12));
      }
      expect(dec.trend[1], isNull);
      expect(dec.trend[10], isNull);
      for (var j = 0; j < 4; j++) {
        expect(dec.seasonalIndices[j], closeTo(pattern[j], 1e-12));
      }
      // Remainder is 0, so F_S = 1.
      expect(dec.seasonalStrength, closeTo(1, 1e-12));
    });

    test('a straight line is degenerate (no variation around the trend)', () {
      final dec = classicalDecomposition([
        for (var t = 0; t < 48; t++) 3 + 0.1 * t,
      ], period: 12);
      expect(dec.isDegenerate, isTrue);
      expect(dec.seasonalStrength, 0);
    });

    test('refuses fewer than two seasons and period < 2', () {
      expect(
        () => classicalDecomposition(List.filled(23, 1.0), period: 12),
        throwsArgumentError,
      );
      expect(
        () => classicalDecomposition(List.filled(30, 1.0), period: 1),
        throwsArgumentError,
      );
    });
  });

  group(
    'root check on 10^4 random polynomials (numpy.roots; mpmath arbiter)',
    () {
      // Columns: coefficients | min |root| by mpmath at 60 digits | by numpy.
      final rows = _rows('roots_reference.txt');
      final polys = [
        for (final r in rows)
          (
            r[0].split(',').map(double.parse).toList(),
            double.parse(r[1]),
            double.parse(r[2]),
          ),
      ];

      test('fixture size', () => expect(polys.length, 10000));

      for (final radius in [1.0, 1.001, 1.01, math.pow(1.01, 12).toDouble()]) {
        test('step-down verdict at radius $radius', () {
          var compared = 0;
          var inBand = 0;
          var numpyWrong = 0;
          final mismatches = <String>[];
          for (final (a, exact, np) in polys) {
            if ((exact - radius).abs() <= 1e-9) {
              inBand++;
              continue;
            }
            compared++;
            if ((np > radius) != (exact > radius)) numpyWrong++;
            if (allRootsOutside(a, radius: radius) != (exact > radius)) {
              mismatches.add('$a min|z|=$exact');
            }
          }
          // ignore: avoid_print
          print(
            'radius $radius: $compared compared, $inBand in the ±1e-9 band, '
            '${mismatches.length} mismatches; numpy.roots itself on the wrong '
            'side in $numpyWrong',
          );
          expect(mismatches, isEmpty);
          expect(compared + inBand, 10000);
        });
      }

      test('minimum root modulus', () {
        var worst = 0.0;
        var worstHigh = 0.0;
        var worstNumpy = 0.0;
        for (final (a, exact, np) in polys) {
          final e = _rel(minRootModulus(a)!, exact);
          worstNumpy = math.max(worstNumpy, _rel(np, exact));
          if (a.length <= 2) {
            worst = math.max(worst, e);
          } else {
            worstHigh = math.max(worstHigh, e);
          }
        }
        // ignore: avoid_print
        print(
          'max rel. error vs 60 digits: degree<=2 $worst, degree 3..5 '
          '$worstHigh (numpy.roots: $worstNumpy)',
        );
        // Clustered (near-multiple) roots are ill-conditioned: a coefficient
        // perturbation of 1 ulp moves them by ~sqrt(ulp) or more, so no
        // double-precision method does better than ~1e-8 there.
        expect(worst, lessThan(1e-7));
        expect(worstHigh, lessThan(1e-6));
      });

      test('hand cases', () {
        expect(allRootsOutside(const []), isTrue);
        expect(minRootModulus(const []), isNull);
        expect(minRootModulus(const [0.0, 0.0]), isNull);
        expect(allRootsOutside(const [0.5]), isTrue); // root 2
        expect(allRootsOutside(const [2.0]), isFalse); // root 0.5
        expect(allRootsOutside(const [0.5], radius: 2.5), isFalse);
        // AR(2) stationarity triangle.
        expect(allRootsOutside(const [0.5, 0.3]), isTrue);
        expect(allRootsOutside(const [0.8, 0.3]), isFalse);
        // Complex pair: 1 − 0.25 z² → |z| = 2... with a_1 = 0.
        expect(minRootModulus(const [0.0, 0.25]), closeTo(2, 1e-15));
        expect(minRootModulus(const [0.0, -0.25]), closeTo(2, 1e-15));
        // Unit root exactly on the boundary is not "outside".
        expect(allRootsOutside(const [1.0]), isFalse);
      });
    },
  );

  group('information criteria (FPP3 §9.6)', () {
    test('formulas', () {
      const ll = -100.0;
      const k = 3;
      const t = 50;
      expect(
        InformationCriterion.aic.evaluate(
          logLikelihood: ll,
          parameters: k,
          observations: t,
        ),
        206,
      );
      expect(
        InformationCriterion.aicc.evaluate(
          logLikelihood: ll,
          parameters: k,
          observations: t,
        ),
        closeTo(206 + 24 / 46, 1e-12),
      );
      expect(
        InformationCriterion.bic.evaluate(
          logLikelihood: ll,
          parameters: k,
          observations: t,
        ),
        closeTo(200 + 3 * math.log(50), 1e-12),
      );
    });

    test('AICc undefined for T − K − 1 ≤ 0', () {
      expect(
        () => InformationCriterion.aicc.evaluate(
          logLikelihood: 1,
          parameters: 4,
          observations: 5,
        ),
        throwsArgumentError,
      );
    });
  });
}
