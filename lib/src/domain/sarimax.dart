/// Regression with (seasonal) ARIMA errors — SARIMAX (domain layer).
///
/// The model is
///
/// ```text
/// y_t = μ + Σ_j β_j · x_{j,t} + u_t,    u_t ~ SARIMA(p,d,q)(P,D,Q)[s]
/// ```
///
/// with the regressors differenced together with `y` (the standard
/// "regression with ARIMA errors" formulation — the same model as statsmodels'
/// `SARIMAX(..., trend='n')` with the intercept as a column of ones, and as R's
/// `arima(xreg = ...)`). The intercept μ is estimated only when nothing is
/// differenced (`d + D == 0`) and [fitSarimax]'s `includeMean` is set.
///
/// Like the rest of the package this knows nothing about FFI and nothing about
/// any application domain: which regressors to use, and what they mean, is the
/// caller's business.
library;

import 'dart:math' as math;
import 'dart:typed_data';

import '../ffi/wrapper/ctsa_native.dart' as native;
import '../ffi/wrapper/tseries_exceptions.dart';
import 'arima.dart';
import 'internal/arima_support.dart';

export '../ffi/wrapper/tseries_exceptions.dart'
    show ExogColumnDefect, ExogenousColumnError, ModelTooLargeError;

/// Estimation method for [fitSarimax] / [sarimaxForecast].
///
/// Deliberately a separate enum from [ArimaMethod]: ctsa's SARIMAX estimators
/// are not the same set as its SARIMA ones (there is no Box-Jenkins here, and
/// there is a pure-MLE start that the SARIMA path does not offer), and ctsa
/// numbers them differently. The mapping onto ctsa is done once, in the native
/// shim.
enum SarimaxMethod {
  /// Conditional sum of squares for starting values, then exact maximum
  /// likelihood. The default, and the same estimator [ArimaMethod.mle] runs on
  /// the SARIMA path — with no regressors the two give identical fits.
  cssMle(0),

  /// Exact maximum likelihood started from zero ARMA coefficients (regression
  /// coefficients start from OLS). Slower to converge than [cssMle]; useful
  /// when the CSS start lands somewhere unhelpful (e.g. status 10/12).
  mle(1),

  /// Conditional sum of squares only. Fastest; no AIC and no coefficient
  /// covariance (see [SarimaxFitResult.aic], [SarimaxFitResult.covariance]).
  css(2);

  const SarimaxMethod(this.code);

  /// The tseries shim's method code (`TSERIES_SARIMAX_METHOD_*`).
  final int code;
}

/// What to do with a regressor column that cannot be estimated (see
/// [ExogColumnDefect]).
enum ExogPolicy {
  /// Throw an [ExogenousColumnError] naming the first defective column (and
  /// carrying the verdict for every column). The default: a defective
  /// regressor is almost always a data-preparation bug worth hearing about.
  strict,

  /// Drop defective columns, fit the rest, and report each dropped column in
  /// [SarimaxFitResult.regressors] with a null coefficient and its
  /// [SarimaxRegressor.defect]. If every column is dropped the result is a
  /// plain SARIMA fit.
  dropDegenerate,
}

/// One regressor of a SARIMAX fit, in the caller's column order.
class SarimaxRegressor {
  const SarimaxRegressor({
    required this.name,
    required this.coefficient,
    required this.standardError,
    required this.defect,
  });

  /// The column's name (its key in the `exog` map; `drift` for the drift
  /// column added by `includeDrift`).
  final String name;

  /// The estimated coefficient β, or null when the column was dropped
  /// ([defect] is non-null).
  final double? coefficient;

  /// Standard error of [coefficient] from the coefficient covariance, or null
  /// when the column was dropped or no trustworthy covariance exists (see
  /// [SarimaxFitResult.covariance]).
  final double? standardError;

  /// Why the column was dropped, or null when it was estimated.
  final ExogColumnDefect? defect;

  /// Whether the column took part in the fit.
  bool get isUsed => defect == null;

  @override
  String toString() => isUsed
      ? '$name: $coefficient (se ${standardError ?? 'n/a'})'
      : '$name: dropped (${defect!.name})';
}

/// Covariance matrix of the estimated coefficients of a SARIMAX fit.
///
/// Parameters are named `ar.L1..`, `ma.L1..`, `ar.S.L{s}..`, `ma.S.L{s}..`,
/// `intercept`, then the regressors by name — only the ones actually
/// estimated. Signs follow the reported coefficients (so the MA entries are in
/// ctsa's `1 − θB` convention, like [SarimaxFitResult.ma]).
///
/// It is the inverse of ctsa's finite-difference Hessian of the (σ²-profiled)
/// log-likelihood at the optimum — the inverse observed information, which is
/// what statsmodels reports with `cov_type='oim'`.
class CoefficientCovariance {
  CoefficientCovariance._(this.parameters, this._m)
    : _index = {for (var i = 0; i < parameters.length; i++) parameters[i]: i};

  /// Parameter names, in matrix order.
  final List<String> parameters;
  final List<double> _m;
  final Map<String, int> _index;

  int _at(String name) =>
      _index[name] ??
      (throw ArgumentError.value(name, 'name', 'not an estimated parameter'));

  /// Covariance of parameters [a] and [b].
  double covariance(String a, String b) =>
      _m[_at(a) * parameters.length + _at(b)];

  /// Variance of parameter [name].
  double variance(String name) => covariance(name, name);

  /// Standard error (square root of the variance) of parameter [name].
  double standardError(String name) => math.sqrt(variance(name));

  /// The matrix as rows, in [parameters] order (a fresh copy).
  List<List<double>> toRows() {
    final k = parameters.length;
    return [
      for (var i = 0; i < k; i++)
        List<double>.unmodifiable(_m.sublist(i * k, (i + 1) * k)),
    ];
  }
}

/// The result of a fixed-order SARIMAX fit: coefficients, regressors,
/// diagnostics and an optional forecast.
class SarimaxFitResult {
  const SarimaxFitResult({
    required this.order,
    required this.method,
    required this.ar,
    required this.ma,
    required this.seasonalAr,
    required this.seasonalMa,
    required this.intercept,
    required this.interceptStandardError,
    required this.regressors,
    required this.sigma2,
    required this.logLikelihood,
    required this.aic,
    required this.covariance,
    required this.forecast,
    required this.standardErrors,
    required this.retval,
  });

  /// The order of the ARIMA error process.
  final ArimaOrder order;

  /// The estimator that was requested.
  final SarimaxMethod method;

  /// Non-seasonal AR coefficients (φ), length [ArimaOrder.p].
  final List<double> ar;

  /// Non-seasonal MA coefficients (θ), length [ArimaOrder.q], in ctsa's
  /// `1 − θB` convention — the opposite sign to R and statsmodels, exactly as
  /// [SarimaFitResult.ma].
  final List<double> ma;

  /// Seasonal AR coefficients (Φ), length [ArimaOrder.seasonalP].
  final List<double> seasonalAr;

  /// Seasonal MA coefficients (Θ), length [ArimaOrder.seasonalQ], ctsa's sign
  /// convention as [ma].
  final List<double> seasonalMa;

  /// The estimated intercept μ, or null when none was estimated (`d + D > 0`,
  /// or `includeMean: false`).
  final double? intercept;

  /// Standard error of [intercept], or null (see [covariance]).
  final double? interceptStandardError;

  /// Every regressor the caller supplied, in the caller's order — estimated or
  /// dropped. This is the map of which columns were actually used.
  final List<SarimaxRegressor> regressors;

  /// Residual variance estimate (σ²).
  final double sigma2;

  /// Log-likelihood, or null when ctsa did not compute a trustworthy one: only
  /// for [status] == [ArimaFitStatus.probableSuccess]. Exact for
  /// [SarimaxMethod.cssMle]/[SarimaxMethod.mle]; a CSS log-likelihood for
  /// [SarimaxMethod.css].
  final double? logLikelihood;

  /// Akaike Information Criterion, or null when not trustworthy: non-null only
  /// for [SarimaxMethod.cssMle]/[SarimaxMethod.mle] with [status] ==
  /// [ArimaFitStatus.probableSuccess]. Counts every estimated coefficient
  /// (ARMA terms, intercept, used regressors) plus σ².
  final double? aic;

  /// Coefficient covariance, or null when it is not trustworthy: non-null only
  /// for [SarimaxMethod.cssMle]/[SarimaxMethod.mle] with [status] ==
  /// [ArimaFitStatus.probableSuccess] and a finite matrix with a positive
  /// diagonal. ctsa's CSS-path Hessian is mis-scaled, so [SarimaxMethod.css]
  /// never has one.
  final CoefficientCovariance? covariance;

  /// Optional forecast (empty if none was requested).
  final List<double> forecast;

  /// Standard errors of [forecast] (same length). These account for the ARIMA
  /// error process only; the future regressor values are treated as known.
  final List<double> standardErrors;

  /// The raw ctsa `retval` for this fit. Reported, never enforced.
  final int retval;

  /// [retval] interpreted.
  ArimaFitStatus get status => ArimaFitStatus.fromNative(retval);

  /// The regressors that took part in the fit.
  Iterable<SarimaxRegressor> get usedRegressors =>
      regressors.where((r) => r.isUsed);

  /// Coefficient by regressor name (null for a dropped column). Throws if
  /// [name] is not a regressor of this fit.
  double? coefficient(String name) => regressors
      .firstWhere(
        (r) => r.name == name,
        orElse: () =>
            throw ArgumentError.value(name, 'name', 'no such regressor'),
      )
      .coefficient;
}

/// The result of [sarimaxForecast]: point forecast, standard errors, the
/// confidence band, the fit that produced them and the ladder record.
class SarimaxForecast {
  const SarimaxForecast({
    required this.point,
    required this.standardErrors,
    required this.ciLower,
    required this.ciUpper,
    required this.confidenceLevel,
    required this.fit,
    required this.attempts,
  });

  /// Point forecast, one value per step ahead.
  final List<double> point;

  /// Standard error of each forecast step (future regressors taken as known).
  final List<double> standardErrors;

  /// `point[i] − z · standardErrors[i]` for the two-sided normal `z` of
  /// [confidenceLevel].
  final List<double> ciLower;

  /// `point[i] + z · standardErrors[i]`.
  final List<double> ciUpper;

  /// The confidence level of the band (e.g. 0.95 → z ≈ 1.96).
  final double confidenceLevel;

  /// The accepted fit: coefficients, regressors, diagnostics, covariance.
  final SarimaxFitResult fit;

  /// The order that actually produced the forecast (the accepted rung).
  ArimaOrder get order => fit.order;

  /// ctsa's status of the accepted fit (reported, not enforced).
  int get retval => fit.retval;

  /// [retval] interpreted.
  ArimaFitStatus get status => fit.status;

  /// Every rung of the ladder that was considered, in order; the accepted one
  /// is last. Same semantics as [ArimaForecast.attempts].
  final List<ArimaAttempt> attempts;
}

/// Name of the column [fitSarimax] adds for `includeDrift`.
const String sarimaxDriftColumn = 'drift';

final RegExp _reservedName = RegExp(r'^(ar|ma)\.(S\.)?L\d+$');

/// Validated, column-ordered regressor data.
typedef _Exog = ({
  List<String> names,
  List<Float64List> columns,
  List<Float64List> future,
});

_Exog _prepareExog({
  required Map<String, List<num>> exog,
  required Map<String, List<num>>? futureExog,
  required int n,
  required int horizon,
  required bool includeDrift,
}) {
  final names = <String>[];
  final columns = <Float64List>[];
  final future = <Float64List>[];
  for (final MapEntry(key: name, value: values) in exog.entries) {
    if (name.isEmpty ||
        name == 'intercept' ||
        _reservedName.hasMatch(name) ||
        (includeDrift && name == sarimaxDriftColumn)) {
      throw ArgumentError.value(
        name,
        'exog',
        'regressor name is empty or reserved (intercept, ar.L*, ma.L*, '
            'ar.S.L*, ma.S.L*${includeDrift ? ', $sarimaxDriftColumn' : ''})',
      );
    }
    if (values.length != n) {
      throw ArgumentError(
        'exog column "$name" has ${values.length} values, expected $n '
        '(the series length)',
      );
    }
    names.add(name);
    columns.add(toFiniteFloat64(values, what: 'exog["$name"]'));
  }

  if (exog.isNotEmpty && horizon > 0 && futureExog == null) {
    throw ArgumentError(
      'futureExog is required to forecast with regressors: the forecast needs '
      'the regressor values for each of the $horizon future steps',
    );
  }
  if (futureExog != null) {
    final extra = futureExog.keys.where((k) => !exog.containsKey(k)).toList();
    if (extra.isNotEmpty) {
      throw ArgumentError(
        'futureExog has columns not present in exog: ${extra.join(', ')}',
      );
    }
    for (final name in names) {
      final values = futureExog[name];
      if (values == null) {
        throw ArgumentError('futureExog is missing column "$name"');
      }
      if (values.length != horizon) {
        throw ArgumentError(
          'futureExog["$name"] has ${values.length} values, expected exactly '
          '$horizon (the horizon)',
        );
      }
      future.add(toFiniteFloat64(values, what: 'futureExog["$name"]'));
    }
  }

  if (includeDrift) {
    names.add(sarimaxDriftColumn);
    columns.add(
      Float64List.fromList([for (var t = 1; t <= n; t++) t.toDouble()]),
    );
    future.add(
      Float64List.fromList([
        for (var t = 1; t <= horizon; t++) (n + t).toDouble(),
      ]),
    );
  }
  return (names: names, columns: columns, future: future);
}

void _validateDrift(ArimaOrder order, bool includeDrift) {
  if (includeDrift && order.d + order.seasonalD > 1) {
    throw ArgumentError(
      'includeDrift requires d + D <= 1 (got $order): a linear time trend '
      'differences to a constant once and to zero twice, so with more '
      'differencing a drift term cannot be estimated',
    );
  }
}

List<String> _parameterNames(
  ArimaOrder o,
  bool meanEstimated,
  List<String> exogNames,
) => [
  for (var i = 1; i <= o.p; i++) 'ar.L$i',
  for (var i = 1; i <= o.q; i++) 'ma.L$i',
  for (var i = 1; i <= o.seasonalP; i++) 'ar.S.L${i * o.seasonalPeriod}',
  for (var i = 1; i <= o.seasonalQ; i++) 'ma.S.L${i * o.seasonalPeriod}',
  if (meanEstimated) 'intercept',
  ...exogNames,
];

/// Fits a regression with (seasonal) ARIMA errors to [series].
///
/// * [exog] — the regressors, one column per entry, each the same length as
///   [series]; the map's iteration order is the column order (which matters
///   for collinearity: the later member of a collinear group is the one
///   flagged). Pass `const {}` for a plain SARIMA fit.
/// * [futureExog] — the regressors' values for each of the [horizon] future
///   steps, with the same keys. Required when [horizon] > 0 and [exog] is not
///   empty: the forecast is conditional on them.
/// * [order] — the error process `(p,d,q)(P,D,Q)[s]`; seasonal orders are fully
///   supported, within the memory ceiling (see [ModelTooLargeError]).
/// * [includeMean] — estimate an intercept when `d + D == 0` (ignored
///   otherwise, as there is nothing to estimate it from).
/// * [includeDrift] — add a regressor named [sarimaxDriftColumn] holding the
///   time index `1..n` (future: `n+1..n+h`). With `d + D == 1` its coefficient
///   is the drift (the mean of the differenced series); with `d + D == 0` it is
///   a deterministic linear trend. Not allowed with `d + D > 1`.
/// * [exogPolicy] — what to do with columns that cannot be estimated under
///   this order (see [ExogPolicy], [ExogColumnDefect]).
///
/// Regressors are screened before the native fit: every value (of the series
/// too, also after differencing) must be finite and within ±1e150,
/// and each column — after the same differencing as the series — must be
/// neither zero nor (within 1e-8) a linear combination of the intercept and
/// earlier columns. At least 2 residual degrees of freedom must remain after
/// differencing and every estimated coefficient.
///
/// Throws [ArgumentError] for unusable input — including
/// [ExogenousColumnError] (strict policy) and [ModelTooLargeError] — and
/// [TseriesNumericException] when the native fit produces no usable model.
SarimaxFitResult fitSarimax(
  List<num> series,
  Map<String, List<num>> exog, {
  required ArimaOrder order,
  Map<String, List<num>>? futureExog,
  int horizon = 0,
  SarimaxMethod method = SarimaxMethod.cssMle,
  bool includeMean = true,
  bool includeDrift = false,
  ExogPolicy exogPolicy = ExogPolicy.strict,
}) {
  validateArimaOrder(order);
  _validateDrift(order, includeDrift);
  if (horizon < 0) {
    throw ArgumentError.value(horizon, 'horizon', 'must be non-negative');
  }
  final y = toFiniteFloat64(series);
  final x = _prepareExog(
    exog: exog,
    futureExog: futureExog,
    n: y.length,
    horizon: horizon,
    includeDrift: includeDrift,
  );
  return _fit(
    y,
    x,
    order: order,
    horizon: horizon,
    method: method,
    includeMean: includeMean,
    exogPolicy: exogPolicy,
  );
}

SarimaxFitResult _fit(
  Float64List y,
  _Exog x, {
  required ArimaOrder order,
  required int horizon,
  required SarimaxMethod method,
  required bool includeMean,
  required ExogPolicy exogPolicy,
}) {
  final raw = native.sarimaxFit(
    y,
    exog: x.columns,
    futureExog: horizon > 0 ? x.future : const [],
    names: x.names,
    p: order.p,
    d: order.d,
    q: order.q,
    s: order.seasonalPeriod,
    seasonalP: order.seasonalP,
    seasonalD: order.seasonalD,
    seasonalQ: order.seasonalQ,
    method: method.code,
    includeMean: includeMean,
    dropDegenerate: exogPolicy == ExogPolicy.dropDegenerate,
    horizon: horizon,
  );

  // Build the covariance over the parameters actually estimated. The native
  // layout covers every supplied column, with NaN rows for dropped ones.
  final allNames = _parameterNames(order, raw.meanEstimated, x.names);
  final fixed = allNames.length - x.names.length;
  final keep = <int>[
    for (var i = 0; i < fixed; i++) i,
    for (var j = 0; j < x.names.length; j++)
      if (raw.columnDefects[j] == null) fixed + j,
  ];
  CoefficientCovariance? cov;
  final v = raw.vcov;
  if (v != null && keep.isNotEmpty) {
    final k = allNames.length;
    final m = <double>[
      for (final a in keep)
        for (final b in keep) v[a * k + b],
    ];
    // Defensive: never hand out a covariance with a non-finite or
    // non-positive variance, whatever the native side claimed.
    var ok = allFinite(m);
    for (var i = 0; ok && i < keep.length; i++) {
      if (!(m[i * keep.length + i] > 0)) ok = false;
    }
    if (ok) {
      cov = CoefficientCovariance._(
        List<String>.unmodifiable([for (final i in keep) allNames[i]]),
        List<double>.unmodifiable(m),
      );
    }
  }

  return SarimaxFitResult(
    order: order,
    method: method,
    ar: raw.phi,
    ma: raw.theta,
    seasonalAr: raw.seasonalPhi,
    seasonalMa: raw.seasonalTheta,
    intercept: raw.meanEstimated ? raw.mean : null,
    interceptStandardError: raw.meanEstimated && cov != null
        ? cov.standardError('intercept')
        : null,
    regressors: List<SarimaxRegressor>.unmodifiable([
      for (var j = 0; j < x.names.length; j++)
        SarimaxRegressor(
          name: x.names[j],
          coefficient: raw.beta[j],
          standardError: raw.columnDefects[j] == null && cov != null
              ? cov.standardError(x.names[j])
              : null,
          defect: raw.columnDefects[j],
        ),
    ]),
    sigma2: raw.sigma2,
    logLikelihood: raw.logLikelihood,
    aic: raw.aic,
    covariance: cov,
    forecast: raw.forecast,
    standardErrors: raw.standardErrors,
    retval: raw.retval,
  );
}

/// Forecasts [horizon] steps of a regression with (seasonal) ARIMA errors,
/// with a confidence band and an order-degradation ladder — the SARIMAX
/// counterpart of [arimaForecast], with the same preconditions (a regular
/// series, [horizon] in steps) and the same ladder semantics.
///
/// [orders] are tried in the given order; the first that meets the numeric
/// length floor (≈10 observations per estimated coefficient, regressors
/// included, on top of what differencing consumes) and fits to a usable model
/// is accepted. Every rung is recorded in [SarimaxForecast.attempts]. A rung
/// whose order is infeasible ([ModelTooLargeError]) is skipped like a rung
/// below the floor.
///
/// Regressor problems are NOT a reason to fall back: under
/// [ExogPolicy.strict] an [ExogenousColumnError] from any rung is rethrown
/// immediately, because the ladder exists to trade model complexity for
/// robustness, not to route around bad input. Use
/// [ExogPolicy.dropDegenerate] to have each rung drop the columns that are
/// degenerate under its own differencing instead.
///
/// The band is `point ± z · standardErrors` with the two-sided normal `z` for
/// [confidenceLevel]; it is conditional on [futureExog] being the true future
/// regressor values.
///
/// Throws [ArgumentError] for bad arguments (including a series below the
/// floor of every rung) and [TseriesNumericException] when every attempted
/// rung failed.
SarimaxForecast sarimaxForecast(
  List<num> series,
  Map<String, List<num>> exog, {
  required List<ArimaOrder> orders,
  required int horizon,
  Map<String, List<num>>? futureExog,
  double confidenceLevel = 0.95,
  SarimaxMethod method = SarimaxMethod.cssMle,
  bool includeMean = true,
  bool includeDrift = false,
  ExogPolicy exogPolicy = ExogPolicy.strict,
}) {
  if (orders.isEmpty) {
    throw ArgumentError.value(orders, 'orders', 'must not be empty');
  }
  if (horizon < 1) {
    throw ArgumentError.value(horizon, 'horizon', 'must be at least 1');
  }
  validateConfidenceLevel(confidenceLevel);
  for (final order in orders) {
    validateArimaOrder(order);
    _validateDrift(order, includeDrift);
  }
  final y = toFiniteFloat64(series);
  final x = _prepareExog(
    exog: exog,
    futureExog: futureExog,
    n: y.length,
    horizon: horizon,
    includeDrift: includeDrift,
  );
  final n = y.length;
  final regressors = x.names.length;

  int floorOf(ArimaOrder o) => arimaMinLength(o, extraParameters: regressors);
  final smallestFloor = orders.map(floorOf).reduce(math.min);
  if (n < smallestFloor) {
    throw ArgumentError.value(
      n,
      'series.length',
      'below the numeric floor: the smallest candidate order with $regressors '
          'regressor(s) requires at least $smallestFloor points',
    );
  }

  final z = normalTwoSidedZ(confidenceLevel);
  final attempts = <ArimaAttempt>[];
  for (final order in orders) {
    final floor = floorOf(order);
    if (n < floor) {
      attempts.add(
        ArimaAttempt(
          order: order,
          accepted: false,
          note: 'series length $n below numeric floor $floor',
        ),
      );
      continue;
    }
    try {
      final fit = _fit(
        y,
        x,
        order: order,
        horizon: horizon,
        method: method,
        includeMean: includeMean,
        exogPolicy: exogPolicy,
      );
      if (!allFinite(fit.forecast) || !allFinite(fit.standardErrors)) {
        attempts.add(
          ArimaAttempt(
            order: order,
            accepted: false,
            retval: fit.retval,
            note: 'native returned non-finite forecast',
          ),
        );
        continue;
      }
      final point = fit.forecast;
      final se = fit.standardErrors;
      attempts.add(
        ArimaAttempt(order: order, accepted: true, retval: fit.retval),
      );
      return SarimaxForecast(
        point: point,
        standardErrors: se,
        ciLower: List<double>.unmodifiable([
          for (var i = 0; i < horizon; i++) point[i] - z * se[i],
        ]),
        ciUpper: List<double>.unmodifiable([
          for (var i = 0; i < horizon; i++) point[i] + z * se[i],
        ]),
        confidenceLevel: confidenceLevel,
        fit: fit,
        attempts: List<ArimaAttempt>.unmodifiable(attempts),
      );
    } on TseriesNumericException catch (e) {
      attempts.add(
        ArimaAttempt(
          order: order,
          accepted: false,
          retval: e.retval,
          note: 'fit did not converge to a usable model',
        ),
      );
    } on ExogenousColumnError {
      rethrow;
    } on ArgumentError catch (e) {
      attempts.add(
        ArimaAttempt(
          order: order,
          accepted: false,
          note: 'rejected by native fit (${e.message})',
        ),
      );
    }
  }

  throw TseriesNumericException(
    'no valid SARIMAX forecast for this series; all ${orders.length} candidate '
    'order(s) failed: ${attempts.join('; ')}',
  );
}
