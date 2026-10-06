/// Native-memory leak check for the SARIMAX path.
///
/// Repeats every distinct exit of `tseries_sarimax_fit` thousands of times and
/// asserts the process's resident set does not grow with the count. The
/// threshold is deliberately coarse (RSS is noisy): it catches leaks of the
/// sizes this path could have — a forgotten work/design buffer is n·r·8 bytes,
/// i.e. kilobytes per fit, megabytes per thousand — not single bytes. The
/// byte-exact measurement (private bytes over 20 000+ fits per path, from a C
/// harness against the same shim) is recorded in the CHANGELOG for 0.7.0.
library;

import 'dart:io';
import 'dart:math' as math;

import 'package:test/test.dart';
import 'package:tseries/tseries.dart';

void main() {
  test(
    'no native memory growth across all SARIMAX exits',
    () {
      final rnd = math.Random(1);
      const n = 96;
      final x1 = [
        for (var t = 0; t < n; t++) math.sin(t / 4) + rnd.nextDouble(),
      ];
      final x2 = [
        for (var t = 0; t < n; t++) math.cos(t / 7) + rnd.nextDouble(),
      ];
      var level = 0.0;
      final y = [
        for (var t = 0; t < n; t++)
          (level += rnd.nextDouble() - 0.5) + 0.8 * x1[t] - 0.4 * x2[t],
      ];
      final explosive = <double>[10];
      for (var i = 1; i < 80; i++) {
        explosive.add(explosive[i - 1] * 1.15 + rnd.nextDouble() - 0.5);
      }
      final ex = [for (var i = 0; i < 80; i++) math.sin(i / 3)];
      const order = ArimaOrder(p: 1, d: 1, q: 1);
      final fut = {
        'x1': [0.1, 0.2, 0.3],
        'x2': [0.0, 0.1, 0.2],
      };

      void oneRound() {
        // 1. success + forecast + covariance
        fitSarimax(
          y,
          {'x1': x1, 'x2': x2},
          order: order,
          horizon: 3,
          futureExog: fut,
        );
        // 2. drop policy (compaction buffers + dropped column)
        fitSarimax(
          y,
          {'x1': x1, 'z': List.filled(n, 0), 'x2': x2},
          order: order,
          exogPolicy: ExogPolicy.dropDegenerate,
        );
        // 3. strict defect: returns after screening allocations
        try {
          fitSarimax(y, {'x1': x1, 'c': List.filled(n, 1)}, order: order);
        } on ExogenousColumnError {
          // expected
        }
        // 4. ctsa early return after its CSS step (status 4, before MLE)
        fitSarimax(explosive, {
          'x': ex,
        }, order: const ArimaOrder(p: 1, d: 0, q: 0));
        // 5. CSS method
        fitSarimax(y, {'x1': x1}, order: order, method: SarimaxMethod.css);
        // 6. degrees-of-freedom refusal after screening
        try {
          fitSarimax(y.sublist(0, 6), {
            'x1': x1.sublist(0, 6),
            'x2': x2.sublist(0, 6),
          }, order: order);
        } on ArgumentError {
          // expected
        }
      }

      for (var i = 0; i < 100; i++) {
        oneRound();
      }
      final before = ProcessInfo.currentRss;
      const rounds = 1500;
      for (var i = 0; i < rounds; i++) {
        oneRound();
      }
      final grown = ProcessInfo.currentRss - before;
      // 9000 native calls. A leak of even 1 KB on any one path would add ~1.5 MB;
      // allow 4 MB for allocator/GC noise.
      expect(
        grown,
        lessThan(4 * 1024 * 1024),
        reason: 'RSS grew by $grown bytes over $rounds rounds',
      );
    },
    timeout: const Timeout(Duration(minutes: 5)),
  );
}
