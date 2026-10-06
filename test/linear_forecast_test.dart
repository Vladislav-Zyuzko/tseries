import 'package:test/test.dart';
import 'package:tseries/tseries.dart';

void main() {
  group('linearForecast — exact line (residuals = 0)', () {
    // y = 2x + 1 sampled at x = 0..3 → [1, 3, 5, 7].
    // OLS must recover slope 2, intercept 1 exactly; residuals are all zero so
    // s = 0 and every standard error is 0. The forecast is the line continued.
    test('recovers slope, intercept and continues the line', () {
      final f = linearForecast([1, 3, 5, 7], window: 4, horizon: 3);

      expect(f.slope, closeTo(2.0, 1e-12));
      expect(f.intercept, closeTo(1.0, 1e-12));

      // x_last = 3, forecast at x = 4, 5, 6 → 9, 11, 13.
      expect(f.point, hasLength(3));
      expect(f.point[0], closeTo(9.0, 1e-9));
      expect(f.point[1], closeTo(11.0, 1e-9));
      expect(f.point[2], closeTo(13.0, 1e-9));

      // Exact fit → s = 0 → all SE are exactly 0.
      expect(f.se, hasLength(3));
      for (final e in f.se) {
        expect(e, closeTo(0.0, 1e-12));
      }
    });

    test('only the last `window` points enter the regression', () {
      // Garbage prefix, then the clean y = 2x + 1 tail. window = 4 must ignore
      // the prefix and fit only [1, 3, 5, 7].
      final f = linearForecast([100, 100, 1, 3, 5, 7], window: 4, horizon: 1);
      expect(f.slope, closeTo(2.0, 1e-12));
      expect(f.intercept, closeTo(1.0, 1e-12));
      expect(f.point.single, closeTo(9.0, 1e-9));
    });
  });

  group('linearForecast — standard error formula & growth', () {
    // Hand-computed reference on y = [1, 2, 2, 5], window = 4:
    //   x̄ = 1.5, ȳ = 2.5, Sxx = 5, Sxy = 6 → slope = 1.2, intercept = 0.7.
    //   RSS = 0.09 + 0.01 + 1.21 + 0.49 = 1.80 → s = √(1.80 / 2) = √0.9.
    //   se[1]: x = 4, term = 1 + 1/4 + 2.5²/5 = 2.5 → se = √0.9·√2.5 = √2.25 = 1.5
    //   se[2]: x = 5, term = 1 + 1/4 + 3.5²/5 = 3.7 → se = √0.9·√3.7 = √3.33
    test('matches the closed-form OLS fit', () {
      final f = linearForecast([1, 2, 2, 5], window: 4, horizon: 2);
      expect(f.slope, closeTo(1.2, 1e-12));
      expect(f.intercept, closeTo(0.7, 1e-12));
      expect(f.point[0], closeTo(5.5, 1e-9)); // 0.7 + 1.2·4
      expect(f.point[1], closeTo(6.7, 1e-9)); // 0.7 + 1.2·5
    });

    test('SE matches the prediction-interval formula and grows with h', () {
      final f = linearForecast([1, 2, 2, 5], window: 4, horizon: 2);
      expect(f.se[0], closeTo(1.5, 1e-9));
      expect(f.se[1], closeTo(1.8248287590, 1e-9)); // √3.33
      // SE strictly increases with the forecast horizon (points sit further from
      // the centre of the fitting window).
      expect(f.se[1], greaterThan(f.se[0]));
    });
  });

  group('linearForecast — output shape', () {
    test('point and se lengths equal the horizon', () {
      final f = linearForecast([1, 3, 5, 7, 9], window: 5, horizon: 6);
      expect(f.point, hasLength(6));
      expect(f.se, hasLength(6));
    });
  });

  group('linearForecast — input guards (ArgumentError)', () {
    test('window < 3 (SE undefined) is rejected', () {
      expect(
        () => linearForecast([1, 2, 3], window: 2, horizon: 1),
        throwsArgumentError,
      );
    });

    test('series shorter than window is rejected', () {
      expect(
        () => linearForecast([1, 2], window: 4, horizon: 1),
        throwsArgumentError,
      );
    });

    test('series shorter than the minimum 3 is rejected (default window)', () {
      expect(() => linearForecast([1, 2], horizon: 1), throwsArgumentError);
    });

    test('horizon of 0 is rejected', () {
      expect(
        () => linearForecast([1, 3, 5, 7], window: 4, horizon: 0),
        throwsArgumentError,
      );
    });

    test('negative horizon is rejected', () {
      expect(
        () => linearForecast([1, 3, 5, 7], window: 4, horizon: -1),
        throwsArgumentError,
      );
    });

    test('NaN in the fitting window is rejected', () {
      expect(
        () => linearForecast([1, 2, double.nan, 4], window: 4, horizon: 1),
        throwsArgumentError,
      );
    });

    test('Infinity in the fitting window is rejected', () {
      expect(
        () => linearForecast([1, 2, 3, double.infinity], window: 4, horizon: 1),
        throwsArgumentError,
      );
    });
  });
}
