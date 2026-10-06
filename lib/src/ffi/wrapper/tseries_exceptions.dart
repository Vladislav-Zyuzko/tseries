/// Exceptions thrown by the tseries native wrapper.
///
/// These translate ctsa/shim C status codes (negative ints) into typed Dart
/// errors so callers never see a raw status code or a `Pointer`. Programmer
/// errors (bad arguments) surface as [ArgumentError]; genuine numerical
/// failures of the native model surface as [TseriesNumericException].
library;

/// Base class for failures originating in the native time-series core.
sealed class TseriesException implements Exception {
  const TseriesException(this.message);

  /// Human-readable description of what went wrong.
  final String message;

  @override
  String toString() => 'TseriesException: $message';
}

/// The native model produced a non-finite (NaN/Inf) coefficient or forecast, or
/// otherwise failed to converge to a usable result.
///
/// This is ctsa's most important known failure mode: the wrapper refuses to
/// hand back numbers it cannot vouch for. Treat this as "no forecast available
/// for this input", not as a crash.
final class TseriesNumericException extends TseriesException {
  const TseriesNumericException(super.message, {this.retval});

  /// The raw ctsa fit-status code (`retval`) of the attempt that failed, or
  /// null if the native fit never ran.
  ///
  /// Carried on the exception on purpose: a fit that the wrapper rejects still
  /// has a status worth observing, and without this the code would be lost
  /// exactly on the paths where it is most diagnostic. Interpret it via
  /// `ArimaFitStatus.fromNative`.
  final int? retval;

  @override
  String toString() => retval == null
      ? 'TseriesNumericException: $message'
      : 'TseriesNumericException: $message (ctsa retval: $retval)';
}

/// An internal invariant of the native boundary was violated (e.g. the shim
/// reported a NULL argument or failed to allocate). Indicates a bug in the
/// wrapper rather than in caller input; should not normally be reachable.
final class TseriesInternalException extends TseriesException {
  const TseriesInternalException(super.message);

  @override
  String toString() => 'TseriesInternalException: $message';
}

/// The requested model order is infeasible: ctsa's exact-likelihood state for
/// it would exceed the package's memory ceiling ([limitBytes], 64 MiB).
///
/// ctsa's state grows roughly with the **fourth power** of
/// `max(p + s·P, q + s·Q + 1)`, so it is the seasonal AR/MA terms at a long
/// period that blow up: a single seasonal MA term at `s = 288` would need ~7 GB
/// per likelihood evaluation. ctsa does not check its allocations, so without
/// this refusal the outcome would be an out-of-memory kill, not an error.
///
/// This is an [ArgumentError] — the order itself is the problem, independently
/// of the data — so the degradation ladders of `arimaForecast` /
/// `sarimaxForecast` treat it like any other order that cannot be fitted and
/// move on to the next candidate.
final class ModelTooLargeError extends ArgumentError {
  ModelTooLargeError({
    required this.order,
    required this.estimatedBytes,
    required this.limitBytes,
  }) : super(
         'the exact-likelihood state of $order would need about '
         '${_mib(estimatedBytes)} MiB (limit ${_mib(limitBytes)} MiB). Seasonal '
         'AR/MA terms at a long period are the usual cause: reduce P/Q, use a '
         'shorter period, or model the seasonality with regressors instead.',
       );

  /// Human-readable description of the offending order.
  final String order;

  /// The estimated peak state size in bytes (saturates at 2^63-1 for orders
  /// too large to even estimate).
  final int estimatedBytes;

  /// The ceiling it was compared against.
  final int limitBytes;

  static String _mib(int bytes) => (bytes / (1024 * 1024)).toStringAsFixed(1);
}

/// The series is too short for the requested order **and estimation method**:
/// the fit is refused before reaching the native core.
///
/// Raised where ctsa would otherwise run on an input it cannot handle — e.g.
/// an MLE/CSS fit whose differenced series is not longer than `p + s·P` (its
/// CSS step would have no residual left and read uninitialised memory, giving
/// run-to-run different results), or a seasonal Box-Jenkins fit on too few
/// seasons. [minimumLength] is the smallest series length the fit accepts and
/// [reason] names the rule.
///
/// This is an [ArgumentError] — the order is incompatible with this length, as
/// with [ModelTooLargeError] — so the degradation ladder of `arimaForecast`
/// skips such an order (with a note) and moves on.
final class SeriesTooShortError extends ArgumentError {
  SeriesTooShortError({
    required this.order,
    required this.length,
    required this.minimumLength,
    required this.reason,
  }) : super(
         'series of length $length is too short for $order: at least '
         '$minimumLength observations are required ($reason)',
       );

  /// Human-readable description of the order.
  final String order;

  /// Length of the series that was passed.
  final int length;

  /// Smallest series length accepted for this order and method.
  final int minimumLength;

  /// The rule that sets [minimumLength].
  final String reason;
}

/// Why a regressor (exogenous) column cannot be estimated. Columns are judged
/// **after** the same differencing (d, and D at the seasonal period) that the
/// model applies to them — which is why a column can be fine for one order and
/// defective for another.
enum ExogColumnDefect {
  /// Every value of the column is exactly zero.
  allZero,

  /// The column is not zero, but differences to (numerically) zero — e.g. a
  /// constant under `d ≥ 1`, a linear trend under `d ≥ 2`, a pattern repeating
  /// with the seasonal period under `D ≥ 1`. It carries no information about
  /// the differenced series.
  constantAfterDifferencing,

  /// A constant column while an intercept is also estimated (`d + D == 0` with
  /// `includeMean`): the two are the same regressor.
  duplicatesIntercept,

  /// The column is (numerically) a linear combination of the intercept and the
  /// columns before it. The **later** member of a collinear group is the one
  /// flagged, so reorder the columns to choose which one survives.
  collinear,
}

/// A regressor column is unusable and the caller asked for a strict check
/// (`ExogPolicy.strict`, the default).
///
/// Carries the first offending column ([column], [name], [defect]) and the
/// verdict for **every** column ([defects], null for a usable one), so all the
/// problems are visible in one round trip. Use `ExogPolicy.dropDegenerate` to
/// have such columns dropped (and reported) instead.
final class ExogenousColumnError extends ArgumentError {
  ExogenousColumnError({
    required this.column,
    required this.columnName,
    required this.defect,
    required List<ExogColumnDefect?> defects,
    required String orderDescription,
  }) : defects = List<ExogColumnDefect?>.unmodifiable(defects),
       super(
         'regressor "$columnName" (column $column) cannot be estimated under '
         '$orderDescription: ${_describe(defect)}',
       );

  /// Index of the first defective column, in the caller's column order.
  final int column;

  /// Name of that column. (Not `name`: [ArgumentError.name] is the name of
  /// the invalid *argument*, which here is `exog`.)
  final String columnName;

  /// What is wrong with it.
  final ExogColumnDefect defect;

  /// Per-column verdict, in the caller's column order; null = usable.
  final List<ExogColumnDefect?> defects;

  static String _describe(ExogColumnDefect d) => switch (d) {
    ExogColumnDefect.allZero => 'it is identically zero',
    ExogColumnDefect.constantAfterDifferencing =>
      'it differences to zero under this order (constant, trend or seasonal '
          'pattern absorbed by d/D)',
    ExogColumnDefect.duplicatesIntercept =>
      'it is constant and duplicates the estimated intercept',
    ExogColumnDefect.collinear =>
      'it is a linear combination of the intercept and/or earlier columns',
  };
}

/// The native core could not be reached at all: the `ctsa` code asset failed to
/// load, or its symbols failed to resolve, on this platform/ABI.
///
/// This is a **deployment/packaging** failure, categorically different from
/// [TseriesNumericException] (which means the native core ran fine but this
/// input has no good model). Distinguishing the two is the whole point: one is
/// "the build is broken on this device", the other is "this data is hard".
///
/// [cause] carries the underlying error (typically the `dart:ffi` lookup
/// failure) for logs.
final class TseriesNativeUnavailableException extends TseriesException {
  const TseriesNativeUnavailableException(super.message, {this.cause});

  /// The underlying error that prevented the native core from being reached.
  final Object? cause;

  @override
  String toString() => cause == null
      ? 'TseriesNativeUnavailableException: $message'
      : 'TseriesNativeUnavailableException: $message (cause: $cause)';
}
