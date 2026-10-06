/// Contract tests for the public native health check.
///
/// The point of [checkNativeCore] is to answer "is the native core deployed and
/// working *here*" independently of any data, so a caller can tell a broken
/// build apart from a hard series. These tests pin that contract on the host;
/// the same call is what a device/diagnostics screen runs to get the same answer
/// on an ABI we cannot reach from here.
library;

import 'package:test/test.dart';
import 'package:tseries/tseries.dart';

void main() {
  group('checkNativeCore', () {
    test('returns healthy on a platform where the native core is built', () {
      final health = checkNativeCore();
      expect(health.smokeValue, 42);
    });

    test('is idempotent and cheap enough to call repeatedly', () {
      // A diagnostics screen may call this on every rebuild; it must not depend
      // on being called once (no lazy one-shot init to get wrong).
      for (var i = 0; i < 100; i++) {
        expect(checkNativeCore().smokeValue, 42);
      }
    });

    test('describes itself usefully for a log/diagnostics line', () {
      expect(checkNativeCore().toString(), contains('smoke=42'));
    });

    test('a healthy native core means forecast failures are data failures', () {
      // The whole reason the check exists: with the core proven up, a
      // TseriesNumericException from a forecast can be trusted to mean "this
      // input has no good model" rather than "the build is broken".
      expect(checkNativeCore().smokeValue, 42);
      expect(
        () => arimaForecast(
          List<double>.filled(60, 5.5), // constant series: degenerate
          orders: const [ArimaOrder(p: 2, d: 1, q: 1)],
          horizon: 3,
        ),
        throwsA(isA<TseriesException>()),
      );
    });
  });

  group('TseriesNativeUnavailableException', () {
    test('is a TseriesException so callers can catch one family', () {
      const e = TseriesNativeUnavailableException('boom');
      expect(e, isA<TseriesException>());
    });

    test('surfaces its cause in toString for logs', () {
      final e = TseriesNativeUnavailableException(
        'could not load',
        cause: ArgumentError('no symbol'),
      );
      expect(e.toString(), contains('could not load'));
      expect(e.toString(), contains('no symbol'));
    });
  });
}
