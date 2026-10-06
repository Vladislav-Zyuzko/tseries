/// Safe wrapper (layer 2) over the raw ctsa/shim bindings.
///
/// This is the ONLY code that touches `Pointer`s. It owns every native
/// allocation (scoped with [using]/[Arena], freed on every path including
/// throws), validates arguments and native return values, and translates C
/// status codes into typed Dart exceptions. Nothing above this layer imports
/// `dart:ffi`, and no `Pointer` ever crosses the function boundaries below —
/// inputs and outputs are plain Dart types.
///
/// Native buffers are read and written ONLY element by element through the
/// `Pointer` (`ptr[i]`), never through a `Pointer.asTypedList` view. Dart VM
/// 3.10 JIT (`dart run`, `dart test`, Flutter debug builds) mis-materializes a
/// scalar-replaced external typed-data view on deoptimization: the object is
/// rebuilt with its data address stored as a tagged integer, i.e. pointing to
/// twice the real address. The next access through the view then reads (or,
/// for a write view, writes) wild memory. It was observed as a 0xC0000005 in
/// `sarimaxFit` when the first fit with a regressor lazily deoptimized the
/// optimized closure while its `diag` view was live. `Pointer` objects are
/// materialized correctly, so element access through them is safe. AOT
/// (release builds) never deoptimizes and was not affected. See CHANGELOG and
/// test/no_typed_data_views_test.dart, which keeps views out of this library.
library;

import 'dart:ffi' as ffi;
import 'dart:typed_data';

import 'package:ffi/ffi.dart';

import '../bindings/tseries_ctsa_bindings.g.dart' as b;
import 'tseries_exceptions.dart';

/// Raw result of a fixed-order (S)ARIMA fit. Field lists have the lengths
/// implied by the requested order; [forecast]/[standardErrors] are empty when
/// no forecast was requested.
///
/// [logLikelihood]/[aic] are null when ctsa did not actually compute them for
/// the requested estimator (see the shim header); [retval] is ctsa's own fit
/// status, reported verbatim and never acted on here.
typedef SarimaFitRaw = ({
  List<double> phi,
  List<double> theta,
  List<double> seasonalPhi,
  List<double> seasonalTheta,
  double mean,
  double sigma2,
  double? logLikelihood,
  double? aic,
  List<double> forecast,
  List<double> standardErrors,
  int retval,
});

/// Maps a native diagnostic to null when the shim marked it "not computed".
///
/// The shim writes NaN — never a plausible-looking substitute — for a
/// diagnostic ctsa did not genuinely produce, so NaN here means "unavailable",
/// not "zero".
double? _nullIfNaN(double v) => v.isNaN ? null : v;

/// The value `tseries_smoke()` returns when the native core is healthy. Defined
/// by the C shim (see native/tseries_ctsa.c); an arbitrary sentinel chosen so a
/// zeroed/garbage return cannot masquerade as success.
const int kSmokeExpected = 42;

/// Calls the native smoke check: exercises `arima_init`/`arima_free` inside the
/// vendored ctsa and returns its sentinel.
///
/// This is the lowest-cost proof that the whole native pipeline works on this
/// device: the `ctsa` code asset was built for this ABI, bundled, loaded, its
/// symbol resolved, and ctsa can allocate and free a model.
///
/// Throws [TseriesNativeUnavailableException] if the symbol cannot be resolved
/// (asset missing for this ABI / build hook did not run) or if ctsa reports it
/// could not initialise a model.
int nativeSmoke() {
  final int value;
  try {
    value = b.tseries_smoke();
  } on Object catch (error) {
    // A failed @Native lookup surfaces here (ArgumentError on most platforms).
    // Never let it escape as an opaque low-level error: this is precisely the
    // "native core is not deployed on this device" case callers must be able to
    // tell apart from a numeric failure.
    throw TseriesNativeUnavailableException(
      'the native ctsa core could not be loaded on this platform: the `ctsa` '
      'code asset is missing or its symbols failed to resolve',
      cause: error,
    );
  }
  if (value != kSmokeExpected) {
    throw TseriesNativeUnavailableException(
      'the native ctsa core loaded but its smoke check returned $value '
      '(expected $kSmokeExpected): ctsa could not initialise a model',
    );
  }
  return value;
}

/// Translate a native status code into an exception (or return normally on OK).
///
/// [retval] is ctsa's own fit status when the native fit ran; it is attached to
/// the numeric exception so callers can still observe the code of an attempt
/// that was rejected.
Never _throwForStatus(int status, String op, {int? retval}) {
  final message = switch (status) {
    b.TSERIES_ERR_NON_FINITE => throw TseriesNumericException(
      '$op: the native model returned a non-finite (NaN/Inf) result; '
      'the fit did not converge to a usable model for this input.',
      retval: retval,
    ),
    b.TSERIES_ERR_DEGENERATE => throw TseriesNumericException(
      '$op: the native model is degenerate (residual variance <= 0); it fits '
      'this input exactly and would report zero forecast uncertainty, so it '
      'is not a usable model. Typical cause: a constant/flatlined series.',
      retval: retval,
    ),
    b.TSERIES_ERR_NULL_ARG =>
      'a required native argument was null (internal wrapper bug)',
    b.TSERIES_ERR_INIT_FAILED =>
      'the native model failed to initialise/allocate',
    _ => 'unexpected native status code $status',
  };
  throw TseriesInternalException('$op: $message');
}

/// Fits a fixed-order (Seasonal) ARIMA model and optionally forecasts
/// [horizon] steps.
///
/// [series] must contain only finite values. Orders must be non-negative and
/// consistent with the seasonal period [s] (pass `s == 0` and
/// `seasonalP/D/Q == 0` for a non-seasonal ARIMA). [method] is 0 = MLE,
/// 1 = CSS, 2 = Box-Jenkins.
///
/// Throws [ArgumentError] for invalid arguments and [TseriesNumericException]
/// if the native fit produces non-finite output.
SarimaFitRaw sarimaFit(
  Float64List series, {
  required int p,
  required int d,
  required int q,
  int s = 0,
  int seasonalP = 0,
  int seasonalD = 0,
  int seasonalQ = 0,
  int method = 0,
  int horizon = 0,
}) {
  final n = series.length;
  if (n < 3) {
    throw ArgumentError.value(n, 'series.length', 'must be at least 3');
  }
  if (horizon < 0) {
    throw ArgumentError.value(horizon, 'horizon', 'must be non-negative');
  }
  for (final v in [p, d, q, s, seasonalP, seasonalD, seasonalQ]) {
    if (v < 0) throw ArgumentError.value(v, 'order', 'must be non-negative');
  }

  return using((arena) {
    final seriesPtr = arena<ffi.Double>(n);
    _fill(seriesPtr, series);

    final phiPtr = p > 0 ? arena<ffi.Double>(p) : ffi.nullptr;
    final thetaPtr = q > 0 ? arena<ffi.Double>(q) : ffi.nullptr;
    final bigPhiPtr = seasonalP > 0
        ? arena<ffi.Double>(seasonalP)
        : ffi.nullptr;
    final bigThetaPtr = seasonalQ > 0
        ? arena<ffi.Double>(seasonalQ)
        : ffi.nullptr;
    final diagPtr = arena<ffi.Double>(4);
    final forecastPtr = horizon > 0 ? arena<ffi.Double>(horizon) : ffi.nullptr;
    final stdErrPtr = horizon > 0 ? arena<ffi.Double>(horizon) : ffi.nullptr;
    final retvalPtr = arena<ffi.Int>(1);

    final status = b.tseries_sarima_fit(
      seriesPtr,
      n,
      p,
      d,
      q,
      s,
      seasonalP,
      seasonalD,
      seasonalQ,
      method,
      horizon,
      phiPtr,
      thetaPtr,
      bigPhiPtr,
      bigThetaPtr,
      diagPtr,
      forecastPtr,
      stdErrPtr,
      retvalPtr,
    );

    // Read the status even on failure paths: the shim writes it as soon as ctsa
    // returns, so a rejected fit still carries a meaningful code.
    final retval = retvalPtr.value;

    if (status != b.TSERIES_OK) {
      if (status == b.TSERIES_ERR_BAD_SIZE) {
        throw ArgumentError(
          'series too short ($n) for the requested order '
          '(p=$p d=$d q=$q P=$seasonalP D=$seasonalD Q=$seasonalQ s=$s)',
        );
      }
      if (status == b.TSERIES_ERR_TOO_LARGE) {
        throw _tooLarge(
          p: p,
          d: d,
          q: q,
          s: s,
          seasonalP: seasonalP,
          seasonalD: seasonalD,
          seasonalQ: seasonalQ,
          withForecast: horizon > 0,
        );
      }
      if (status == b.TSERIES_ERR_BAD_PARAM) {
        if (method == b.TSERIES_METHOD_BOX_JENKINS &&
            p + q + seasonalP + seasonalQ == 0 &&
            d + seasonalD > 0) {
          // The shim refuses this: with no ARMA term and no mean (the series
          // is differenced) Box-Jenkins has nothing to estimate, and ctsa's
          // estimator writes out of bounds on it (PROVENANCE.md, Local
          // modifications 12).
          throw ArgumentError(
            'Box-Jenkins estimation has nothing to estimate for the '
            'differenced white-noise order ${_describeOrder(p, d, q, s, seasonalP, seasonalD, seasonalQ)}: '
            'use the MLE or CSS method',
          );
        }
        throw ArgumentError('inconsistent model order for sarimaFit');
      }
      _throwForStatus(status, 'sarimaFit', retval: retval);
    }

    return (
      phi: _copy(phiPtr, p),
      theta: _copy(thetaPtr, q),
      seasonalPhi: _copy(bigPhiPtr, seasonalP),
      seasonalTheta: _copy(bigThetaPtr, seasonalQ),
      mean: diagPtr[0],
      sigma2: diagPtr[1],
      logLikelihood: _nullIfNaN(diagPtr[2]),
      aic: _nullIfNaN(diagPtr[3]),
      forecast: _copy(forecastPtr, horizon),
      standardErrors: _copy(stdErrPtr, horizon),
      retval: retval,
    );
  });
}

/// Estimated peak bytes of ctsa's exact-likelihood state for an order, as the
/// native guard computes it (see `TSERIES_MAX_STATE_BYTES` in the shim header).
int modelStateBytes({
  required int p,
  required int d,
  required int q,
  int s = 0,
  int seasonalP = 0,
  int seasonalD = 0,
  int seasonalQ = 0,
  bool withForecast = true,
}) => b.tseries_model_state_bytes(
  p,
  d,
  q,
  s,
  seasonalP,
  seasonalD,
  seasonalQ,
  withForecast ? 1 : 0,
);

/// The ceiling [modelStateBytes] is compared against.
const int maxModelStateBytes = b.TSERIES_MAX_STATE_BYTES;

ModelTooLargeError _tooLarge({
  required int p,
  required int d,
  required int q,
  required int s,
  required int seasonalP,
  required int seasonalD,
  required int seasonalQ,
  required bool withForecast,
}) => ModelTooLargeError(
  order: _describeOrder(p, d, q, s, seasonalP, seasonalD, seasonalQ),
  estimatedBytes: modelStateBytes(
    p: p,
    d: d,
    q: q,
    s: s,
    seasonalP: seasonalP,
    seasonalD: seasonalD,
    seasonalQ: seasonalQ,
    withForecast: withForecast,
  ),
  limitBytes: maxModelStateBytes,
);

String _describeOrder(int p, int d, int q, int s, int sp, int sd, int sq) =>
    (sp | sd | sq) != 0
    ? 'ARIMA($p,$d,$q)($sp,$sd,$sq)[$s]'
    : 'ARIMA($p,$d,$q)';

/// Raw result of a SARIMAX fit. Per-column lists ([beta], [columnStatus])
/// follow the caller's column order; [beta] is null for a column that was not
/// estimated. [vcov] is the row-major `k×k` covariance matrix over
/// `[phi, theta, PHI, THETA, mu?, beta_0..beta_{r-1}]` (rows/columns of a
/// dropped column are NaN), or null when ctsa's estimate is not trustworthy.
typedef SarimaxFitRaw = ({
  List<double> phi,
  List<double> theta,
  List<double> seasonalPhi,
  List<double> seasonalTheta,
  List<double?> beta,
  List<ExogColumnDefect?> columnDefects,
  bool meanEstimated,
  double mean,
  double sigma2,
  double? logLikelihood,
  double? aic,
  List<double>? vcov,
  List<double> forecast,
  List<double> standardErrors,
  int retval,
});

ExogColumnDefect? _defectFromNative(int code) => switch (code) {
  b.TSERIES_EXOG_USED => null,
  b.TSERIES_EXOG_ALL_ZERO => ExogColumnDefect.allZero,
  b.TSERIES_EXOG_CONSTANT_AFTER_DIFF =>
    ExogColumnDefect.constantAfterDifferencing,
  b.TSERIES_EXOG_COLLINEAR_INTERCEPT => ExogColumnDefect.duplicatesIntercept,
  b.TSERIES_EXOG_COLLINEAR => ExogColumnDefect.collinear,
  // Never trust the native side to stay within the documented set.
  _ => throw TseriesInternalException(
    'sarimaxFit: unknown regressor status code $code from the native shim',
  ),
};

/// Fits a regression with (seasonal) ARIMA errors and optionally forecasts
/// [horizon] steps.
///
/// [exog] holds the regressor columns (each of length `series.length`),
/// [futureExog] their future values (each of length [horizon]; may be empty
/// when [horizon] is 0), and [names] one label per column for error messages.
/// [method] is a `TSERIES_SARIMAX_METHOD_*` code. With [dropDegenerate] false a
/// defective column throws [ExogenousColumnError]; with it true the column is
/// dropped and reported in `columnDefects`.
///
/// Throws [ArgumentError] (incl. [ModelTooLargeError] / [ExogenousColumnError])
/// for unusable input, and [TseriesNumericException] when the native fit does
/// not produce a usable model.
SarimaxFitRaw sarimaxFit(
  Float64List series, {
  required List<Float64List> exog,
  required List<Float64List> futureExog,
  required List<String> names,
  required int p,
  required int d,
  required int q,
  int s = 0,
  int seasonalP = 0,
  int seasonalD = 0,
  int seasonalQ = 0,
  required int method,
  required bool includeMean,
  required bool dropDegenerate,
  int horizon = 0,
}) {
  final n = series.length;
  final r = exog.length;
  if (n < 3) {
    throw ArgumentError.value(n, 'series.length', 'must be at least 3');
  }
  if (horizon < 0) {
    throw ArgumentError.value(horizon, 'horizon', 'must be non-negative');
  }
  for (final v in [p, d, q, s, seasonalP, seasonalD, seasonalQ]) {
    if (v < 0) throw ArgumentError.value(v, 'order', 'must be non-negative');
  }
  if (names.length != r) {
    throw ArgumentError('names must have one entry per exog column');
  }
  for (var j = 0; j < r; j++) {
    if (exog[j].length != n) {
      throw ArgumentError(
        'exog column "${names[j]}" has ${exog[j].length} values, '
        'expected $n (the series length)',
      );
    }
    if (horizon > 0 &&
        (futureExog.length != r || futureExog[j].length != horizon)) {
      throw ArgumentError(
        'future values of exog column "${names[j]}" must have exactly '
        '$horizon values (the horizon)',
      );
    }
  }
  final orderText = _describeOrder(p, d, q, s, seasonalP, seasonalD, seasonalQ);
  final meanEstimated = includeMean && d + seasonalD == 0;
  final k = p + q + seasonalP + seasonalQ + (meanEstimated ? 1 : 0) + r;

  return using((arena) {
    final seriesPtr = arena<ffi.Double>(n);
    _fill(seriesPtr, series);

    // Column-major, as ctsa expects: column j occupies [j*n, (j+1)*n).
    final xregPtr = r > 0 ? arena<ffi.Double>(n * r) : ffi.nullptr;
    final futurePtr = r > 0 && horizon > 0
        ? arena<ffi.Double>(horizon * r)
        : ffi.nullptr;
    for (var j = 0; j < r; j++) {
      _fill(xregPtr, exog[j], at: j * n);
      if (horizon > 0) _fill(futurePtr, futureExog[j], at: j * horizon);
    }

    final phiPtr = p > 0 ? arena<ffi.Double>(p) : ffi.nullptr;
    final thetaPtr = q > 0 ? arena<ffi.Double>(q) : ffi.nullptr;
    final bigPhiPtr = seasonalP > 0
        ? arena<ffi.Double>(seasonalP)
        : ffi.nullptr;
    final bigThetaPtr = seasonalQ > 0
        ? arena<ffi.Double>(seasonalQ)
        : ffi.nullptr;
    final betaPtr = r > 0 ? arena<ffi.Double>(r) : ffi.nullptr;
    final colStatusPtr = r > 0 ? arena<ffi.Int>(r) : ffi.nullptr;
    final diagPtr = arena<ffi.Double>(4);
    final vcovPtr = k > 0 ? arena<ffi.Double>(k * k) : ffi.nullptr;
    final forecastPtr = horizon > 0 ? arena<ffi.Double>(horizon) : ffi.nullptr;
    final stdErrPtr = horizon > 0 ? arena<ffi.Double>(horizon) : ffi.nullptr;
    final retvalPtr = arena<ffi.Int>(1);
    final badColumnPtr = arena<ffi.Int>(1);

    final status = b.tseries_sarimax_fit(
      seriesPtr,
      n,
      xregPtr,
      r,
      futurePtr,
      horizon,
      p,
      d,
      q,
      s,
      seasonalP,
      seasonalD,
      seasonalQ,
      method,
      includeMean ? 1 : 0,
      dropDegenerate
          ? b.TSERIES_EXOG_POLICY_DROP
          : b.TSERIES_EXOG_POLICY_STRICT,
      phiPtr,
      thetaPtr,
      bigPhiPtr,
      bigThetaPtr,
      betaPtr,
      colStatusPtr,
      diagPtr,
      vcovPtr,
      forecastPtr,
      stdErrPtr,
      retvalPtr,
      badColumnPtr,
    );

    final retval = retvalPtr.value;
    final badColumn = badColumnPtr.value;

    if (status != b.TSERIES_OK) {
      switch (status) {
        case b.TSERIES_ERR_BAD_SIZE:
          throw ArgumentError(
            'series too short ($n) for $orderText with ${r > 0 ? '$r regressor(s)' : 'no regressors'}'
            '${meanEstimated ? ' and an intercept' : ''}: at least '
            '${b.TSERIES_MIN_RESIDUAL_DOF} residual degrees of freedom must '
            'remain after differencing and every estimated coefficient'
            '${seasonalD > 0 && d > 0 && r > 0 ? ', and with seasonal differencing d × (regressors used) must not exceed the differenced length ${n - d - s * seasonalD}' : ''}',
          );
        case b.TSERIES_ERR_BAD_PARAM:
          throw ArgumentError(
            'inconsistent model specification for sarimaxFit ($orderText, '
            'method code $method)',
          );
        case b.TSERIES_ERR_TOO_LARGE:
          throw _tooLarge(
            p: p,
            d: d,
            q: q,
            s: s,
            seasonalP: seasonalP,
            seasonalD: seasonalD,
            seasonalQ: seasonalQ,
            withForecast: horizon > 0,
          );
        case b.TSERIES_ERR_EXOG_DEFECT:
          if (badColumn < 0 || badColumn >= r) {
            throw TseriesInternalException(
              'sarimaxFit: regressor defect reported for out-of-range column '
              '$badColumn',
            );
          }
          final defects = [
            for (var j = 0; j < r; j++) _defectFromNative(colStatusPtr[j]),
          ];
          throw ExogenousColumnError(
            column: badColumn,
            columnName: names[badColumn],
            defect: defects[badColumn]!,
            defects: defects,
            orderDescription: orderText,
          );
        case b.TSERIES_ERR_NON_FINITE when retval == b.TSERIES_RETVAL_NOT_RUN:
          // Rejected before ctsa ran: an input value (or a differenced value)
          // is not finite.
          throw ArgumentError(
            badColumn >= 0 && badColumn < r
                ? 'exog column "${names[badColumn]}" contains a value that is non-finite or beyond '
                      '±${b.TSERIES_MAX_ABS_VALUE}, or exceeds that when differenced '
                      'under $orderText'
                : 'series contains a value that is non-finite or beyond '
                      '±${b.TSERIES_MAX_ABS_VALUE}, or exceeds that when '
                      'differenced under $orderText',
          );
        case b.TSERIES_ERR_NATIVE_COLLINEAR:
          throw TseriesNumericException(
            'sarimaxFit: ctsa rejected the regressors as collinear although '
            "tseries's own screening accepted them. ctsa's rank test is "
            'scale-dependent: regressors are rescaled to the series, but the '
            'intercept column is not, so a series whose differenced scale is '
            'many orders of magnitude away from 1 can trigger this. Rescale '
            'the series.',
            retval: retval,
          );
      }
      _throwForStatus(status, 'sarimaxFit', retval: retval);
    }

    List<double>? vcov;
    if (k > 0) {
      // The shim NaN-fills the whole matrix when it is not trustworthy, and
      // only the rows/columns of dropped regressors otherwise. A finite first
      // estimated diagonal entry therefore distinguishes the two.
      final firstEstimated = _firstEstimatedIndex(
        k: k,
        fixed: k - r,
        columnStatus: colStatusPtr,
        r: r,
      );
      if (firstEstimated != null &&
          !vcovPtr[firstEstimated * k + firstEstimated].isNaN) {
        vcov = _copy(vcovPtr, k * k);
      }
    }
    return (
      phi: _copy(phiPtr, p),
      theta: _copy(thetaPtr, q),
      seasonalPhi: _copy(bigPhiPtr, seasonalP),
      seasonalTheta: _copy(bigThetaPtr, seasonalQ),
      beta: [for (var j = 0; j < r; j++) _nullIfNaN(betaPtr[j])],
      columnDefects: [
        for (var j = 0; j < r; j++) _defectFromNative(colStatusPtr[j]),
      ],
      meanEstimated: meanEstimated,
      mean: diagPtr[0],
      sigma2: diagPtr[1],
      logLikelihood: _nullIfNaN(diagPtr[2]),
      aic: _nullIfNaN(diagPtr[3]),
      vcov: vcov,
      forecast: _copy(forecastPtr, horizon),
      standardErrors: _copy(stdErrPtr, horizon),
      retval: retval,
    );
  });
}

/// Index (into the k×k layout) of the first parameter that was actually
/// estimated, or null when nothing was (k == 0 or every column dropped and no
/// other parameter).
int? _firstEstimatedIndex({
  required int k,
  required int fixed,
  required ffi.Pointer<ffi.Int> columnStatus,
  required int r,
}) {
  if (fixed > 0) return 0;
  for (var j = 0; j < r; j++) {
    if (columnStatus[j] == b.TSERIES_EXOG_USED) return fixed + j;
  }
  return null;
}

/// Copies [len] doubles out of a native buffer into an owned Dart list. Returns
/// an empty list when [len] is zero (and the pointer is null).
///
/// Element access through the `Pointer`, not an `asTypedList` view: see the
/// library comment.
List<double> _copy(ffi.Pointer<ffi.Double> ptr, int len) {
  if (len <= 0) return const <double>[];
  return [for (var i = 0; i < len; i++) ptr[i]];
}

/// Writes [values] into the native buffer [dst] starting at element [at].
///
/// Element access through the `Pointer`, not an `asTypedList` view: see the
/// library comment.
void _fill(ffi.Pointer<ffi.Double> dst, List<double> values, {int at = 0}) {
  for (var i = 0; i < values.length; i++) {
    dst[at + i] = values[i];
  }
}
