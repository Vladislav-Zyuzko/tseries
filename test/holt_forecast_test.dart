import 'dart:typed_data';

import 'package:test/test.dart';
import 'package:tseries/tseries.dart';

/// Raw IEEE-754 bit pattern of [v].
///
/// Used where the claim under test is *bit-for-bit* equality rather than
/// numerical closeness: `==` on doubles already conflates `+0.0` with `-0.0`,
/// and `closeTo` would silently accept a real drift. Comparing the 64-bit
/// pattern is the strongest statement available.
int _bits(double v) => (ByteData(8)..setFloat64(0, v)).getUint64(0);

void main() {
  group('holtForecast — alpha = 1, phi = 0 is exactly persistence', () {
    // The free correctness test from the methodology spec: with alpha = 1 the
    // level collapses onto the latest observation, and with phi = 0 the damped
    // horizon sum is 0, so the trend cannot contribute at any horizon. The
    // forecast must therefore be y_t repeated — bit for bit, not merely close.
    // If this drifts, the primitive is broken; it is not a discovery.
    const series = <double>[5.4, 6.1, 7.8, 7.2, 6.6, 6.9, 8.3, 9.1, 8.4, 7.7];

    test('reproduces y_t bit-for-bit at every horizon (beta = 0.1)', () {
      final f = holtForecast(series, alpha: 1, beta: 0.1, phi: 0, horizon: 6);
      expect(f.point, hasLength(6));
      for (var h = 0; h < 6; h++) {
        expect(
          _bits(f.point[h]),
          equals(_bits(series.last)),
          reason: 'point[$h] must be bit-identical to y_t',
        );
      }
    });

    test('holds for any beta — beta cannot leak in when phi = 0', () {
      // The spec says "alpha = 1, beta = anything, phi = 0". Sweep the whole
      // admitted beta range, endpoints included.
      for (final beta in <double>[0, 0.05, 0.1, 0.2, 0.5, 1]) {
        final f = holtForecast(
          series,
          alpha: 1,
          beta: beta,
          phi: 0,
          horizon: 6,
        );
        for (var h = 0; h < 6; h++) {
          expect(
            _bits(f.point[h]),
            equals(_bits(series.last)),
            reason: 'beta = $beta leaked into point[$h]',
          );
        }
      }
    });

    test('holds regardless of series length and shape', () {
      // A long ramp: persistence must still ignore the (large) trend entirely.
      final ramp = List<double>.generate(200, (i) => 3.0 + 0.25 * i);
      final f = holtForecast(ramp, alpha: 1, beta: 0.2, phi: 0, horizon: 12);
      for (var h = 0; h < 12; h++) {
        expect(_bits(f.point[h]), equals(_bits(ramp.last)));
      }
    });
  });

  group('holtForecast — hand-computed golden reference', () {
    // An independent by-hand trace of the full recursion, seed included. The
    // symmetry tests above (persistence / SES / classical Holt) all probe
    // *degenerate* corners of the parameter space, so none of them would catch
    // an off-by-one in the seed or a swapped level/trend update in the general
    // case. This one does: every intermediate below was computed by hand and
    // the arithmetic is exact in binary (all values are dyadic rationals).
    //
    // series = [1, 2, 3, 4], alpha = beta = phi = 0.5.
    // Seed OLS over min(12, 4) = 4 points of y = 1 + x → slope 1, intercept 1
    //   ⇒ l0 = 1 (fitted value at the first point), b0 = 1.
    // t=0: l = .5·1 + .5·(1 + .5·1)         = 1.25
    //      b = .5·(1.25 − 1) + .5·.5·1      = 0.375
    // t=1: l = .5·2 + .5·(1.25 + .5·.375)   = 1.71875
    //      b = .5·(1.71875 − 1.25) + .25·.375 = 0.328125
    // t=2: l = .5·3 + .5·(1.71875 + .5·.328125) = 2.44140625
    //      b = .5·(2.44140625 − 1.71875) + .25·.328125 = 0.443359375
    // t=3: l = .5·4 + .5·(2.44140625 + .5·.443359375) = 3.33154296875
    //      b = .5·(3.33154296875 − 2.44140625) + .25·.443359375
    //        = 0.555908203125
    // S_1 = .5·(1 − .5)/.5 = 0.5   → point[0] = l + b·0.5  = 3.6094970703125
    // S_2 = .5·(1 − .25)/.5 = 0.75 → point[1] = l + b·0.75 = 3.74847412109375
    test('reproduces a by-hand trace of the recursion exactly', () {
      final f = holtForecast(
        [1.0, 2.0, 3.0, 4.0],
        alpha: 0.5,
        beta: 0.5,
        phi: 0.5,
        horizon: 2,
      );

      expect(_bits(f.point[0]), equals(_bits(3.6094970703125)));
      expect(_bits(f.point[1]), equals(_bits(3.74847412109375)));
    });
  });

  group('holtForecast — phi = 1 is classical Holt', () {
    test('continues a linear ramp exactly', () {
      // y = 10 + 1·t. With alpha = beta = 1 the level locks onto y_t and the
      // trend onto the first difference (= 1), so the forecast must continue
      // the straight line: point[h-1] = y_last + h.
      final ramp = List<double>.generate(30, (i) => 10.0 + i);
      final f = holtForecast(ramp, alpha: 1, beta: 1, phi: 1, horizon: 6);

      for (var h = 1; h <= 6; h++) {
        expect(
          f.point[h - 1],
          closeTo(ramp.last + h, 1e-9),
          reason: 'phi = 1 must not damp the trend at h = $h',
        );
      }
    });

    test('undamped increments are constant (no saturation)', () {
      final ramp = List<double>.generate(30, (i) => 10.0 + 0.4 * i);
      final f = holtForecast(ramp, alpha: 0.5, beta: 0.3, phi: 1, horizon: 10);

      final first = f.point[1] - f.point[0];
      for (var h = 1; h < 9; h++) {
        expect(
          f.point[h + 1] - f.point[h],
          closeTo(first, 1e-9),
          reason: 'classical Holt must extrapolate a straight line',
        );
      }
    });

    test(
      'phi = 1 produces no NaN/Inf (the 1/(1 - phi) branch is not taken)',
      () {
        final ramp = List<double>.generate(30, (i) => 10.0 + 0.4 * i);
        final f = holtForecast(
          ramp,
          alpha: 0.5,
          beta: 0.3,
          phi: 1,
          horizon: 24,
        );
        for (final v in f.point) {
          expect(v.isFinite, isTrue);
        }
      },
    );
  });

  group('holtForecast — phi = 0 is SES (flat forecast at the level)', () {
    test('forecast is flat across the whole horizon', () {
      const series = <double>[5.4, 6.1, 7.8, 7.2, 6.6, 6.9, 8.3, 9.1, 8.4, 7.7];
      final f = holtForecast(
        series,
        alpha: 0.3,
        beta: 0.1,
        phi: 0,
        horizon: 12,
      );

      for (var h = 1; h < 12; h++) {
        expect(
          _bits(f.point[h]),
          equals(_bits(f.point[0])),
          reason: 'phi = 0 must give a flat forecast at l_t',
        );
      }
    });

    test('the flat value is the SES level, not the last observation', () {
      // With alpha = 0.3 the level lags the series, so a flat forecast at l_t
      // must NOT coincide with y_t — this guards against phi = 0 accidentally
      // degenerating into persistence for every alpha.
      const series = <double>[5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 9.0];
      final f = holtForecast(series, alpha: 0.3, beta: 0.1, phi: 0, horizon: 3);
      expect(f.point[0], lessThan(series.last));
      expect(f.point[0], greaterThan(5.0));
    });

    test('a flat series forecasts flat at its own level', () {
      final flat = List<double>.filled(40, 6.5);
      final f = holtForecast(flat, alpha: 0.3, beta: 0.1, phi: 0, horizon: 6);
      for (final v in f.point) {
        expect(v, closeTo(6.5, 1e-12));
      }
    });
  });

  group('holtForecast — damping saturates for 0 < phi < 1', () {
    // On a ramp with a live trend, forecast increments must decay geometrically
    // by exactly phi, and the forecast must converge to l_t + b_t·phi/(1 − phi).
    final ramp = List<double>.generate(60, (i) => 4.0 + 0.3 * i);
    const phi = 0.9;

    test('increments decay geometrically by exactly phi', () {
      final f = holtForecast(
        ramp,
        alpha: 0.3,
        beta: 0.1,
        phi: phi,
        horizon: 40,
      );

      for (var h = 1; h < 20; h++) {
        final prev = f.point[h] - f.point[h - 1];
        final next = f.point[h + 1] - f.point[h];
        expect(
          next / prev,
          closeTo(phi, 1e-9),
          reason: 'increment ratio at h = $h must equal phi',
        );
      }
    });

    test('converges to l_t + b_t·phi/(1 - phi) as h grows', () {
      final f = holtForecast(
        ramp,
        alpha: 0.3,
        beta: 0.1,
        phi: phi,
        horizon: 400,
      );

      // Recover the state algebraically from the forecast itself (black box):
      //   point[0] = l + b·phi
      //   point[1] = l + b·(phi + phi²)  ⇒  point[1] − point[0] = b·phi²
      final b = (f.point[1] - f.point[0]) / (phi * phi);
      final l = f.point[0] - b * phi;
      final limit = l + b * phi / (1 - phi);

      expect(b, greaterThan(0), reason: 'the ramp must leave a positive trend');
      // phi^400 ≈ 1e-19: the tail is numerically at the limit.
      expect(f.point[399], closeTo(limit, 1e-9));
      // ...and it approaches from below, never overshooting it.
      expect(f.point[0], lessThan(limit));
      expect(f.point[399], lessThanOrEqualTo(limit + 1e-12));
    });

    test('damped forecast stays below the undamped one on a rising ramp', () {
      final damped = holtForecast(
        ramp,
        alpha: 0.3,
        beta: 0.1,
        phi: 0.9,
        horizon: 12,
      );
      final undamped = holtForecast(
        ramp,
        alpha: 0.3,
        beta: 0.1,
        phi: 1,
        horizon: 12,
      );

      // The whole point of damping: it cures trend overshoot.
      for (var h = 1; h < 12; h++) {
        expect(
          damped.point[h],
          lessThan(undamped.point[h]),
          reason: 'damping must shrink the trend contribution at h = $h',
        );
      }
    });

    test('no NaN/Inf near the phi = 1 boundary', () {
      for (final p in <double>[0.95, 0.99, 0.999, 0.9999999]) {
        final f = holtForecast(
          ramp,
          alpha: 0.3,
          beta: 0.1,
          phi: p,
          horizon: 24,
        );
        for (final v in f.point) {
          expect(v.isFinite, isTrue, reason: 'phi = $p produced $v');
        }
      }
    });
  });

  group('holtForecast — determinism', () {
    // Trivially true for closed-form Dart arithmetic, but pinned by a test on
    // purpose: determinism is the property that justifies this primitive's
    // existence as a Tier-1 baseline, so it gets a regression guard.
    test('the same input twice yields a bit-identical forecast', () {
      final series = List<double>.generate(
        144,
        (i) => 6.0 + 0.02 * i + (i % 7) * 0.13,
      );

      final a = holtForecast(
        series,
        alpha: 0.3,
        beta: 0.1,
        phi: 0.9,
        horizon: 6,
      );
      final b = holtForecast(
        series,
        alpha: 0.3,
        beta: 0.1,
        phi: 0.9,
        horizon: 6,
      );

      expect(a.point, hasLength(6));
      for (var h = 0; h < 6; h++) {
        expect(
          _bits(a.point[h]),
          equals(_bits(b.point[h])),
          reason: 'point[$h] differed between two identical calls',
        );
      }
    });

    test('is stable across the parameter grid', () {
      final series = List<double>.generate(
        144,
        (i) => 6.0 + 0.02 * i + (i % 7) * 0.13,
      );

      for (final alpha in <double>[0.2, 0.35, 0.5, 0.7]) {
        for (final beta in <double>[0.05, 0.1, 0.2]) {
          for (final phi in <double>[0, 0.7, 0.8, 0.9, 0.95, 1.0]) {
            final a = holtForecast(
              series,
              alpha: alpha,
              beta: beta,
              phi: phi,
              horizon: 6,
            );
            final b = holtForecast(
              series,
              alpha: alpha,
              beta: beta,
              phi: phi,
              horizon: 6,
            );
            for (var h = 0; h < 6; h++) {
              expect(
                _bits(a.point[h]),
                equals(_bits(b.point[h])),
                reason: 'non-determinism at ($alpha, $beta, $phi) h = $h',
              );
              expect(
                a.point[h].isFinite,
                isTrue,
                reason: 'non-finite at ($alpha, $beta, $phi) h = $h',
              );
            }
          }
        }
      }
    });
  });

  group('holtForecast — output shape', () {
    test('point length equals the horizon', () {
      final series = List<double>.generate(20, (i) => 5.0 + i * 0.1);
      expect(
        holtForecast(series, alpha: 0.3, beta: 0.1, phi: 0.9, horizon: 1).point,
        hasLength(1),
      );
      expect(
        holtForecast(series, alpha: 0.3, beta: 0.1, phi: 0.9, horizon: 6).point,
        hasLength(6),
      );
      expect(
        holtForecast(
          series,
          alpha: 0.3,
          beta: 0.1,
          phi: 0.9,
          horizon: 50,
        ).point,
        hasLength(50),
      );
    });

    test('the shortest admissible series (2 points) works', () {
      final f = holtForecast(
        [5.0, 6.0],
        alpha: 0.3,
        beta: 0.1,
        phi: 0.9,
        horizon: 3,
      );
      expect(f.point, hasLength(3));
      for (final v in f.point) {
        expect(v.isFinite, isTrue);
      }
    });
  });

  group('holtForecast — input guards (ArgumentError)', () {
    final series = List<double>.generate(20, (i) => 5.0 + i * 0.1);

    test('horizon of 0 is rejected', () {
      expect(
        () => holtForecast(series, alpha: 0.3, beta: 0.1, phi: 0.9, horizon: 0),
        throwsArgumentError,
      );
    });

    test('negative horizon is rejected', () {
      expect(
        () =>
            holtForecast(series, alpha: 0.3, beta: 0.1, phi: 0.9, horizon: -1),
        throwsArgumentError,
      );
    });

    test('empty series is rejected', () {
      expect(
        () => holtForecast([], alpha: 0.3, beta: 0.1, phi: 0.9, horizon: 1),
        throwsArgumentError,
      );
    });

    test('single-point series is rejected (no slope to seed)', () {
      expect(
        () => holtForecast([5.0], alpha: 0.3, beta: 0.1, phi: 0.9, horizon: 1),
        throwsArgumentError,
      );
    });

    test('NaN in the series is rejected', () {
      expect(
        () => holtForecast(
          [1.0, 2.0, double.nan, 4.0],
          alpha: 0.3,
          beta: 0.1,
          phi: 0.9,
          horizon: 1,
        ),
        throwsArgumentError,
      );
    });

    test('Infinity in the series is rejected', () {
      expect(
        () => holtForecast(
          [1.0, 2.0, 3.0, double.infinity],
          alpha: 0.3,
          beta: 0.1,
          phi: 0.9,
          horizon: 1,
        ),
        throwsArgumentError,
      );
    });

    test('alpha outside (0, 1] is rejected', () {
      for (final alpha in <double>[0, -0.1, 1.1, 2]) {
        expect(
          () => holtForecast(
            series,
            alpha: alpha,
            beta: 0.1,
            phi: 0.9,
            horizon: 1,
          ),
          throwsArgumentError,
          reason: 'alpha = $alpha must be rejected',
        );
      }
    });

    test('beta outside [0, 1] is rejected', () {
      for (final beta in <double>[-0.1, 1.1, 2]) {
        expect(
          () => holtForecast(
            series,
            alpha: 0.3,
            beta: beta,
            phi: 0.9,
            horizon: 1,
          ),
          throwsArgumentError,
          reason: 'beta = $beta must be rejected',
        );
      }
    });

    test('phi outside [0, 1] is rejected', () {
      for (final phi in <double>[-0.1, -1, 1.1, 2]) {
        expect(
          () =>
              holtForecast(series, alpha: 0.3, beta: 0.1, phi: phi, horizon: 1),
          throwsArgumentError,
          reason: 'phi = $phi must be rejected',
        );
      }
    });

    test(
      'non-finite parameters are rejected (NaN defeats a bare range check)',
      () {
        // Every comparison against NaN is false, so `alpha <= 0 || alpha > 1`
        // would wave it straight through into the recursion.
        expect(
          () => holtForecast(
            series,
            alpha: double.nan,
            beta: 0.1,
            phi: 0.9,
            horizon: 1,
          ),
          throwsArgumentError,
        );
        expect(
          () => holtForecast(
            series,
            alpha: 0.3,
            beta: double.nan,
            phi: 0.9,
            horizon: 1,
          ),
          throwsArgumentError,
        );
        expect(
          () => holtForecast(
            series,
            alpha: 0.3,
            beta: 0.1,
            phi: double.nan,
            horizon: 1,
          ),
          throwsArgumentError,
        );
        expect(
          () => holtForecast(
            series,
            alpha: double.infinity,
            beta: 0.1,
            phi: 0.9,
            horizon: 1,
          ),
          throwsArgumentError,
        );
      },
    );

    test('the admitted boundary values are accepted', () {
      expect(
        () => holtForecast(series, alpha: 1, beta: 0, phi: 0, horizon: 1),
        returnsNormally,
      );
      expect(
        () => holtForecast(series, alpha: 1, beta: 1, phi: 1, horizon: 1),
        returnsNormally,
      );
    });
  });
}
