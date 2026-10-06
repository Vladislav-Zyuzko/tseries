/// Health check for the native core (domain layer).
///
/// The rest of the package hides the native boundary on purpose: consumers call
/// [arimaForecast] and cannot tell it is backed by C. That abstraction is right
/// for normal use and actively unhelpful when something is wrong on a device —
/// a missing `ctsa` code asset and a series that simply has no good model both
/// end as "no forecast", and no caller can tell which.
///
/// This library exposes the one thing the abstraction has to let through: a
/// direct, cheap answer to *"is the native core actually here and working on
/// this device?"*, separate from any question about data.
library;

import '../ffi/wrapper/ctsa_native.dart' as native;
import '../ffi/wrapper/tseries_exceptions.dart';

export '../ffi/wrapper/tseries_exceptions.dart'
    show TseriesNativeUnavailableException;

/// A successful native-core health check.
///
/// Only ever produced by [checkNativeCore]; its existence is itself the result
/// (the check throws rather than returning an unhealthy instance), so there is
/// no `isHealthy` flag to forget to test.
class NativeCoreHealth {
  /// @nodoc.
  const NativeCoreHealth({required this.smokeValue});

  /// The sentinel the native shim returned. Carried for diagnostics/logging;
  /// callers should not branch on it — a wrong value throws instead.
  final int smokeValue;

  @override
  String toString() =>
      'native ctsa core reachable and functional (smoke=$smokeValue)';
}

/// Verifies that the native ctsa core is present and working on this device.
///
/// Resolves and calls a native symbol that allocates and frees a real ctsa
/// model. It touches no user data and takes microseconds, so it is safe to call
/// from a diagnostics screen, a start-up check, or a test.
///
/// This is the check to run when a forecast silently does not appear: it splits
/// the two causes that otherwise look identical.
///
/// * Returns [NativeCoreHealth] — the native core is built for this ABI,
///   bundled, loaded, resolved, and ctsa can drive a model. Any forecast
///   failure is then about the **data or the model**, not the build.
/// * Throws [TseriesNativeUnavailableException] — the native core is not
///   reachable on this device: the code asset was not built for this ABI, was
///   not bundled, or its symbols did not resolve. Every native-backed API in
///   this package will fail until that is fixed.
///
/// It deliberately does not report on any *numeric* health: a passing check
/// says nothing about whether a given series will fit. Use the forecast APIs
/// for that, and read a [TseriesNumericException] from them as "this input has
/// no good model" — which a passing [checkNativeCore] proves is the real story.
NativeCoreHealth checkNativeCore() =>
    NativeCoreHealth(smokeValue: native.nativeSmoke());
