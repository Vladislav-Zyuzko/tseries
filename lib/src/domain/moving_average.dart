/// Pure-Dart moving-average primitives.
///
/// This is domain-layer code: trivial, well-understood arithmetic that has no
/// business crossing the FFI boundary. Per the package's design rule, only
/// genuinely hard, validated numerics (ARIMA/SARIMA and friends) are delegated
/// to the native C core; simple rolling statistics stay here in idiomatic Dart.
library;

/// Computes the **simple moving average** (SMA) of [series] over a sliding
/// window of size [window].
///
/// The SMA at each position is the arithmetic mean of the `window` most recent
/// values. Only *complete* windows are emitted (the "valid" convention), so the
/// result has length `series.length - window + 1` when the series is at least
/// as long as the window, and is empty otherwise.
///
/// For a series `[x0, x1, ..., x(n-1)]` and window `w`, element `i` of the
/// result is `mean(x[i], x[i+1], ..., x[i+w-1])` for `i` in `0 .. n - w`.
///
/// ### Examples
///
/// ```dart
/// simpleMovingAverage([1, 2, 3, 4, 5], window: 3); // [2.0, 3.0, 4.0]
/// simpleMovingAverage([2, 4, 6, 8], window: 2);     // [3.0, 5.0, 7.0]
/// simpleMovingAverage([1, 2, 3], window: 1);        // [1.0, 2.0, 3.0]
/// simpleMovingAverage([1, 2], window: 5);           // []  (window too long)
/// ```
///
/// ### Errors
///
/// Throws [ArgumentError] if [window] is not a positive integer. A window
/// larger than the series is *not* an error — it simply yields no complete
/// window and returns an empty list, which composes cleanly in pipelines.
///
/// Runs in O(n) time via an incremental rolling sum (no per-window
/// re-summation). Integer inputs are promoted to `double` in the result.
List<double> simpleMovingAverage(List<num> series, {required int window}) {
  if (window < 1) {
    throw ArgumentError.value(window, 'window', 'must be a positive integer');
  }

  final n = series.length;
  if (window > n) return const <double>[];

  final result = List<double>.filled(n - window + 1, 0);

  // Seed the sum of the first full window.
  var sum = 0.0;
  for (var i = 0; i < window; i++) {
    sum += series[i];
  }
  result[0] = sum / window;

  // Slide: drop the value leaving the window, add the one entering it.
  for (var i = window; i < n; i++) {
    sum += series[i] - series[i - window];
    result[i - window + 1] = sum / window;
  }

  return result;
}
