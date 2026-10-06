import 'package:test/test.dart';
import 'package:tseries/tseries.dart';

void main() {
  group('simpleMovingAverage', () {
    test('window of 3 over a simple integer series', () {
      // Windows: mean(1,2,3)=2, mean(2,3,4)=3, mean(3,4,5)=4.
      expect(simpleMovingAverage([1, 2, 3, 4, 5], window: 3), [2.0, 3.0, 4.0]);
    });

    test('window of 2 yields non-integer means', () {
      expect(simpleMovingAverage([1, 2, 4, 8], window: 2), [1.5, 3.0, 6.0]);
    });

    test('window of 1 is the series itself, promoted to double', () {
      expect(simpleMovingAverage([3, 1, 4, 1, 5], window: 1), [
        3.0,
        1.0,
        4.0,
        1.0,
        5.0,
      ]);
    });

    test('window equal to the series length yields a single mean', () {
      expect(simpleMovingAverage([2, 4, 6, 8], window: 4), [5.0]);
    });

    test('handles double inputs', () {
      expect(simpleMovingAverage([1.0, 2.0], window: 2), [1.5]);
    });

    test('window longer than the series returns empty (not an error)', () {
      expect(simpleMovingAverage([1, 2], window: 5), isEmpty);
    });

    test('empty series returns empty', () {
      expect(simpleMovingAverage(<num>[], window: 3), isEmpty);
    });

    test('result length is n - window + 1', () {
      final out = simpleMovingAverage([10, 20, 30, 40, 50, 60], window: 4);
      expect(out.length, 6 - 4 + 1);
      expect(out, [25.0, 35.0, 45.0]);
    });

    test('rolling sum stays correct over a longer series', () {
      // Cross-check the O(n) rolling implementation against a naive mean.
      final series = List<num>.generate(50, (i) => (i * i) % 7);
      const window = 5;
      final expected = <double>[
        for (var i = 0; i + window <= series.length; i++)
          series.sublist(i, i + window).fold<num>(0, (a, b) => a + b) / window,
      ];
      expect(simpleMovingAverage(series, window: window), expected);
    });

    test('throws ArgumentError on zero window', () {
      expect(
        () => simpleMovingAverage([1, 2, 3], window: 0),
        throwsArgumentError,
      );
    });

    test('throws ArgumentError on negative window', () {
      expect(
        () => simpleMovingAverage([1, 2, 3], window: -2),
        throwsArgumentError,
      );
    });
  });
}
