@ffi.DefaultAsset('package:tseries/ctsa')
library;

import 'dart:ffi' as ffi;

import 'package:test/test.dart';

// Smoke-test binding (hand-written, test-only). Proves the native build hook
// compiled the vendored ctsa + shim, linked them into the `ctsa` code asset,
// and that dart:ffi can resolve and call an exported symbol on this host.
@ffi.Native<ffi.Int32 Function()>(symbol: 'tseries_smoke')
external int tseriesSmoke();

void main() {
  test('native ctsa pipeline links and runs (arima_init/free via shim)', () {
    expect(tseriesSmoke(), 42);
  });
}
