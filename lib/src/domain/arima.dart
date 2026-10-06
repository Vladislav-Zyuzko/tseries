/// Idiomatic ARIMA / SARIMA forecasting API (domain layer).
///
/// This is the clean, framework-agnostic surface consumers use. It knows
/// nothing about FFI: it takes and returns plain Dart types and delegates the
/// heavy numerics to the native ctsa core through the safe wrapper. Trivial
/// statistics (moving averages, EWMA, ...) live elsewhere in pure Dart — only
/// genuinely hard econometrics come through here.
library;

import 'dart:math' as math;

import '../ffi/wrapper/ctsa_native.dart' as native;
import '../ffi/wrapper/tseries_exceptions.dart';
import 'internal/arima_support.dart';

export '../ffi/wrapper/tseries_exceptions.dart'
    show
        TseriesException,
        TseriesNumericException,
        TseriesInternalException,
        ModelTooLargeError,
        SeriesTooShortError;

/// ctsa's own fit-status code (`retval`), reported verbatim.
///
/// **tseries does not reject a fit on the strength of this code.** It reports
/// what the native library produced so callers can measure the real
/// distribution of outcomes before deciding what, if anything, to treat as a
/// failure. A fit can therefore be returned successfully with a status other
/// than [probableSuccess].
///
/// Only [probableSuccess] means the requested estimator ran to completion. The
/// codes are transcribed from the vendored ctsa sources (the "Error Codes"
/// comment in `ctsa.c` and `checkroots_cerr` in `emle.c`).
enum ArimaFitStatus {
  /// ctsa 0 — input error, or the fit never ran.
  notRun(0),

  /// ctsa 1 — upstream's own wording is "**Probable** Success". The only code
  /// for which the coefficients are a completed fit by the requested method.
  probableSuccess(1),

  /// ctsa 4 — the optimiser exhausted its iteration budget. The coefficients
  /// are wherever the search happened to stop, not a converged optimum.
  maxIterations(4),

  /// ctsa 7 — the exogenous regressors are collinear. Never the status of a
  /// returned fit: ctsa estimates nothing at all in this case. tseries screens
  /// regressors itself before ctsa sees them (see `fitSarimax`), so this code
  /// can only surface as the `retval` of a [TseriesNumericException] from the
  /// SARIMAX path, when ctsa's own (scale-dependent) rank test disagrees.
  collinearExogenous(7),

  /// ctsa 10 — the AR roots of the CSS pre-estimate were rejected as
  /// non-stationary. ctsa returns **before** the MLE step, so the coefficients
  /// are CSS estimates even when MLE was requested.
  nonStationaryAr(10),

  /// ctsa 12 — as [nonStationaryAr], but for the seasonal AR part, and with
  /// the same consequence: CSS-only coefficients.
  nonStationarySeasonalAr(12),

  /// ctsa 15 — the optimiser encountered Inf/NaN.
  nonFinite(15),

  /// A code the vendored ctsa does not document. The [code] of this value is a
  /// sentinel, **not** the real code — read the raw int (e.g.
  /// [ArimaAttempt.retval]) to see what ctsa actually returned.
  unknown(-1);

  const ArimaFitStatus(this.code);

  /// The numeric ctsa code, except for [unknown] (see its docs).
  final int code;

  /// Maps a raw ctsa `retval` onto this enum, never throwing: any code the
  /// vendored library did not document maps to [unknown] rather than crashing
  /// the caller. Never trust the native side to stay within a known set.
  static ArimaFitStatus fromNative(int code) => switch (code) {
    0 => notRun,
    1 => probableSuccess,
    4 => maxIterations,
    7 => collinearExogenous,
    10 => nonStationaryAr,
    12 => nonStationarySeasonalAr,
    15 => nonFinite,
    _ => unknown,
  };
}

/// One rung of the degradation ladder tried by [arimaForecast].
///
/// Every candidate order that the ladder considered produces exactly one of
/// these, in ladder order, whether it was accepted, rejected, or skipped. This
/// is what makes the ladder's behaviour measurable: group by [order] and
/// [status] for a per-order histogram of ctsa codes, and filter on [accepted]
/// for the mix of orders actually used.
class ArimaAttempt {
  const ArimaAttempt({
    required this.order,
    required this.accepted,
    this.retval,
    this.note,
  });

  /// The candidate order this rung tried.
  final ArimaOrder order;

  /// Whether this rung produced the returned forecast. Exactly one attempt on a
  /// successful [arimaForecast] has this set.
  final bool accepted;

  /// The raw ctsa `retval` for this rung, or null when the native fit never ran
  /// (the rung was skipped before reaching C — e.g. below the length floor).
  ///
  /// A non-null [retval] with `accepted == false` means ctsa did run and
  /// returned this code, but the wrapper then rejected the fit (non-finite
  /// output). That distinction is deliberate and worth measuring.
  final int? retval;

  /// [retval] interpreted; null when the native fit never ran.
  ArimaFitStatus? get status =>
      retval == null ? null : ArimaFitStatus.fromNative(retval!);

  /// Why this rung was skipped or rejected; null when [accepted].
  final String? note;

  @override
  String toString() {
    final outcome = accepted ? 'accepted' : (note ?? 'rejected');
    final code = retval == null ? 'not run' : '${status!.name}($retval)';
    return '$order: $code — $outcome';
  }
}

/// Estimation method for a fixed-order fit.
enum ArimaMethod {
  /// Maximum likelihood (default; most accurate).
  mle(0),

  /// Conditional sum of squares (faster, approximate).
  css(1),

  /// Box-Jenkins.
  boxJenkins(2);

  const ArimaMethod(this.code);

  /// Native method selector understood by the ctsa shim.
  final int code;
}

/// A (seasonal) ARIMA order `(p,d,q)(P,D,Q)[s]`.
class ArimaOrder {
  const ArimaOrder({
    required this.p,
    required this.d,
    required this.q,
    this.seasonalP = 0,
    this.seasonalD = 0,
    this.seasonalQ = 0,
    this.seasonalPeriod = 0,
  });

  final int p;
  final int d;
  final int q;
  final int seasonalP;
  final int seasonalD;
  final int seasonalQ;

  /// Seasonal period `s` (e.g. 12 for monthly data); 0 means non-seasonal.
  final int seasonalPeriod;

  /// Whether this order carries a seasonal component.
  bool get isSeasonal => seasonalPeriod > 0;

  /// Value equality, so an order can be used directly as a map key — e.g. to
  /// aggregate [ArimaAttempt]s into a per-order histogram of ctsa statuses.
  @override
  bool operator ==(Object other) =>
      identical(this, other) ||
      other is ArimaOrder &&
          p == other.p &&
          d == other.d &&
          q == other.q &&
          seasonalP == other.seasonalP &&
          seasonalD == other.seasonalD &&
          seasonalQ == other.seasonalQ &&
          seasonalPeriod == other.seasonalPeriod;

  @override
  int get hashCode =>
      Object.hash(p, d, q, seasonalP, seasonalD, seasonalQ, seasonalPeriod);

  @override
  String toString() => isSeasonal
      ? 'ARIMA($p,$d,$q)($seasonalP,$seasonalD,$seasonalQ)[$seasonalPeriod]'
      : 'ARIMA($p,$d,$q)';
}

/// The result of a robust ARIMA forecast produced by [arimaForecast].
///
/// Carries the point forecast, its per-step standard errors, the confidence
/// band derived from [confidenceLevel], the model order that actually produced
/// the forecast (see the degradation ladder in [arimaForecast]) and the model
/// diagnostics/coefficients. All lists have length equal to the requested
/// horizon.
///
/// This is a purely numeric result: it knows nothing about the domain the
/// series comes from. Any domain interpretation (thresholds, alerts,
/// probabilities of crossing a level) is the caller's job, computed from
/// [point] and [standardErrors].
class ArimaForecast {
  const ArimaForecast({
    required this.point,
    required this.standardErrors,
    required this.ciLower,
    required this.ciUpper,
    required this.confidenceLevel,
    required this.order,
    required this.ar,
    required this.ma,
    required this.seasonalAr,
    required this.seasonalMa,
    required this.mean,
    required this.sigma2,
    required this.logLikelihood,
    required this.aic,
    required this.retval,
    required this.attempts,
  });

  /// Point forecast, one value per step ahead (length == requested horizon).
  final List<double> point;

  /// Standard error of each forecast step (same length as [point]).
  final List<double> standardErrors;

  /// Lower confidence bound, `point[i] - z * standardErrors[i]`, where `z` is
  /// the two-sided normal critical value for [confidenceLevel].
  final List<double> ciLower;

  /// Upper confidence bound, `point[i] + z * standardErrors[i]`.
  final List<double> ciUpper;

  /// The confidence level of [ciLower]/[ciUpper] (e.g. 0.95 → z ≈ 1.96).
  final double confidenceLevel;

  /// The candidate order that actually produced this forecast. When a
  /// degradation ladder was supplied to [arimaForecast], this is the first
  /// order in the ladder that fitted successfully — not necessarily the first
  /// one requested.
  final ArimaOrder order;

  /// Non-seasonal AR coefficients (φ), length [ArimaOrder.p].
  final List<double> ar;

  /// Non-seasonal MA coefficients (θ), length [ArimaOrder.q].
  final List<double> ma;

  /// Seasonal AR coefficients (Φ), length [ArimaOrder.seasonalP].
  final List<double> seasonalAr;

  /// Seasonal MA coefficients (Θ), length [ArimaOrder.seasonalQ].
  final List<double> seasonalMa;

  /// Fitted mean / intercept term.
  final double mean;

  /// Residual variance estimate (σ²).
  final double sigma2;

  /// Maximised log-likelihood, or null when ctsa did not compute one for the
  /// requested method (see [SarimaFitResult.logLikelihood]).
  final double? logLikelihood;

  /// Akaike Information Criterion (lower is better), or null when it is not
  /// trustworthy (see [SarimaFitResult.aic] for exactly when).
  final double? aic;

  /// The raw ctsa `retval` of the accepted fit. See [status].
  final int retval;

  /// [retval] interpreted. Reported, not enforced: a forecast is returned even
  /// when this is not [ArimaFitStatus.probableSuccess].
  ArimaFitStatus get status => ArimaFitStatus.fromNative(retval);

  /// Every rung of the degradation ladder that was considered, in ladder order,
  /// including the accepted one (the last entry, with
  /// [ArimaAttempt.accepted] set) and any rungs that were skipped or failed
  /// before it.
  ///
  /// This is the raw material for measuring the ladder: a per-order histogram
  /// of ctsa statuses, and the mix of orders actually accepted.
  final List<ArimaAttempt> attempts;
}

/// The result of a fixed-order fit: estimated coefficients, diagnostics, and an
/// optional forecast.
class SarimaFitResult {
  const SarimaFitResult({
    required this.order,
    required this.ar,
    required this.ma,
    required this.seasonalAr,
    required this.seasonalMa,
    required this.mean,
    required this.sigma2,
    required this.logLikelihood,
    required this.aic,
    required this.forecast,
    required this.standardErrors,
    required this.retval,
  });

  final ArimaOrder order;

  /// Non-seasonal AR coefficients (φ), length [ArimaOrder.p].
  final List<double> ar;

  /// Non-seasonal MA coefficients (θ), length [ArimaOrder.q].
  final List<double> ma;

  /// Seasonal AR coefficients (Φ), length [ArimaOrder.seasonalP].
  final List<double> seasonalAr;

  /// Seasonal MA coefficients (Θ), length [ArimaOrder.seasonalQ].
  final List<double> seasonalMa;

  /// Fitted mean / intercept term.
  final double mean;

  /// Residual variance estimate (σ²).
  final double sigma2;

  /// Maximised log-likelihood, or **null when ctsa did not compute one** for
  /// the estimator that was requested.
  ///
  /// Non-null only for [ArimaMethod.mle] and [ArimaMethod.css] with
  /// [status] == [ArimaFitStatus.probableSuccess]. ctsa never computes a
  /// log-likelihood for [ArimaMethod.boxJenkins].
  final double? logLikelihood;

  /// Akaike Information Criterion (lower is better), or **null when it is not
  /// trustworthy**.
  ///
  /// Non-null only for [ArimaMethod.mle] with [status] ==
  /// [ArimaFitStatus.probableSuccess]. It is null in every other case, for two
  /// distinct reasons:
  ///
  /// * ctsa only ever computes an AIC on its MLE branch — for CSS and
  ///   Box-Jenkins the field is simply never assigned.
  /// * On [ArimaFitStatus.nonStationaryAr] / [ArimaFitStatus.nonStationarySeasonalAr]
  ///   ctsa returns before the MLE step, so the AIC it does compute is derived
  ///   from a **CSS** log-likelihood. It is not an MLE AIC and **is not
  ///   comparable with one** — which is precisely what an AIC is for.
  ///
  /// The null is the point: comparing models by AIC is the obvious next step,
  /// and a silently incomparable number would corrupt that comparison without
  /// ever looking wrong. Null forces the question at compile time.
  final double? aic;

  /// Optional forecast (empty if none was requested).
  final List<double> forecast;

  /// Standard errors of [forecast] (same length).
  final List<double> standardErrors;

  /// The raw ctsa `retval` for this fit. See [status].
  ///
  /// Reported, never enforced: this fit is returned regardless of the code.
  final int retval;

  /// [retval] interpreted.
  ArimaFitStatus get status => ArimaFitStatus.fromNative(retval);
}

/// Fits a fixed-order (seasonal) ARIMA model to [series].
///
/// Provide the full [order]. When [horizon] > 0 the result also carries a
/// forecast. [method] selects the estimation method (defaults to MLE).
///
/// Throws [ArgumentError] for invalid input — a [SeriesTooShortError] when
/// [series] is shorter than [order] and [method] allow (see its docs) — and
/// [TseriesNumericException] if the native fit produces non-finite output.
SarimaFitResult fitSarima(
  List<num> series, {
  required ArimaOrder order,
  int horizon = 0,
  ArimaMethod method = ArimaMethod.mle,
}) {
  final values = toFiniteFloat64(series);
  final minimum = sarimaMinimumLength(order, method);
  if (values.length < minimum.length) {
    throw SeriesTooShortError(
      order: '$order (${method.name})',
      length: values.length,
      minimumLength: minimum.length,
      reason: minimum.reason,
    );
  }
  final raw = native.sarimaFit(
    values,
    p: order.p,
    d: order.d,
    q: order.q,
    s: order.seasonalPeriod,
    seasonalP: order.seasonalP,
    seasonalD: order.seasonalD,
    seasonalQ: order.seasonalQ,
    method: method.code,
    horizon: horizon,
  );
  return SarimaFitResult(
    order: order,
    ar: raw.phi,
    ma: raw.theta,
    seasonalAr: raw.seasonalPhi,
    seasonalMa: raw.seasonalTheta,
    mean: raw.mean,
    sigma2: raw.sigma2,
    logLikelihood: raw.logLikelihood,
    aic: raw.aic,
    forecast: raw.forecast,
    standardErrors: raw.standardErrors,
    retval: raw.retval,
  );
}

/// Forecasts [horizon] steps ahead of a fixed-order ARIMA model, with a
/// confidence band and an automatic order-degradation ladder.
///
/// This is the forecast-oriented, robust entry point built on top of the
/// lower-level [fitSarima]. Prefer it when you want a ready-to-use forecast
/// (point + uncertainty band) that tolerates a hard series by falling back to
/// simpler models.
///
/// ## Preconditions (the package does NOT relax these)
///
/// * [series] must already be a **regular** series — evenly spaced, with no
///   gaps or missing samples. Resampling and imputation are the caller's job;
///   this package is time-agnostic and never invents data points.
/// * [horizon] is measured in **steps**, not in time. The caller owns the
///   sampling interval Δt. Far horizons are statistically unreliable regardless
///   of model; the package does not cap [horizon], but the caller should.
///
/// ## Degradation ladder
///
/// [orders] is a list of candidate orders tried **in the given order**. The
/// first one that (a) meets the numeric length floor and (b) fits to a finite
/// model is used, and its order is reported back in [ArimaForecast.order]. This
/// lets a caller express, e.g. `(2,1,1) → (1,1,1) → (0,1,1)`: prefer the richer
/// model, but fall back to a simpler one on a short or difficult series. The
/// package stays generic — it does not know or hardcode any particular ladder.
/// Pass a single-element list for a plain fixed-order forecast.
///
/// ## Confidence band
///
/// [ArimaForecast.ciLower]/[ArimaForecast.ciUpper] are `point ± z · stdErr`,
/// where `z` is the two-sided normal critical value for [confidenceLevel]
/// (0.95 → z ≈ 1.96). The band assumes Gaussian forecast errors, as usual for
/// ARIMA.
///
/// ## Failure modes
///
/// * [ArgumentError] — bad arguments: empty [orders], a malformed order,
///   non-finite values in [series], [horizon] < 1, [confidenceLevel] not in
///   `(0, 1)`, or a [series] shorter than the numeric floor of **every**
///   candidate order (nothing could be fitted).
/// * [TseriesNumericException] — every candidate order that was long enough to
///   attempt failed to converge to a finite model ("no valid forecast for this
///   input").
ArimaForecast arimaForecast(
  List<num> series, {
  required List<ArimaOrder> orders,
  required int horizon,
  double confidenceLevel = 0.95,
  ArimaMethod method = ArimaMethod.mle,
}) {
  if (orders.isEmpty) {
    throw ArgumentError.value(orders, 'orders', 'must not be empty');
  }
  if (horizon < 1) {
    throw ArgumentError.value(horizon, 'horizon', 'must be at least 1');
  }
  if (!(confidenceLevel > 0 && confidenceLevel < 1)) {
    throw ArgumentError.value(
      confidenceLevel,
      'confidenceLevel',
      'must be in the open interval (0, 1)',
    );
  }
  for (final order in orders) {
    validateArimaOrder(order);
  }

  // Validate finiteness once, up front, so a genuine non-finite input surfaces
  // as an ArgumentError immediately instead of being masked as "all candidates
  // failed" by the ladder below.
  toFiniteFloat64(series);

  final n = series.length;

  // If the series is shorter than the floor of even the smallest candidate,
  // nothing can be fitted — a clear ArgumentError, not a ladder exhaustion.
  // The floor of a candidate is its numeric floor or its structural minimum
  // for this method, whichever is larger.
  final smallestFloor = orders
      .map(
        (o) =>
            math.max(arimaMinLength(o), sarimaMinimumLength(o, method).length),
      )
      .reduce(math.min);
  if (n < smallestFloor) {
    throw ArgumentError.value(
      n,
      'series.length',
      'below the numeric floor: the smallest candidate order requires at least '
          '$smallestFloor points',
    );
  }

  final z = normalTwoSidedZ(confidenceLevel);
  final attempts = <ArimaAttempt>[];

  for (final order in orders) {
    final floor = arimaMinLength(order);
    if (n < floor) {
      // Skipped before reaching C, so there is no ctsa status to report.
      attempts.add(
        ArimaAttempt(
          order: order,
          accepted: false,
          note: 'series length $n below numeric floor $floor',
        ),
      );
      continue;
    }
    final structural = sarimaMinimumLength(order, method);
    if (n < structural.length) {
      // Refused before reaching C, like the numeric floor above: ctsa would
      // run on an input it cannot handle (see SeriesTooShortError).
      attempts.add(
        ArimaAttempt(
          order: order,
          accepted: false,
          note:
              'series length $n below structural minimum '
              '${structural.length} for ${method.name}: ${structural.reason}',
        ),
      );
      continue;
    }
    try {
      final fit = fitSarima(
        series,
        order: order,
        horizon: horizon,
        method: method,
      );
      // Defensive: the native shim already rejects non-finite output, but never
      // trust the boundary blindly before handing numbers back to callers.
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
      final stdErr = fit.standardErrors;
      final ciLower = List<double>.generate(
        horizon,
        (i) => point[i] - z * stdErr[i],
        growable: false,
      );
      final ciUpper = List<double>.generate(
        horizon,
        (i) => point[i] + z * stdErr[i],
        growable: false,
      );
      attempts.add(
        ArimaAttempt(order: order, accepted: true, retval: fit.retval),
      );
      return ArimaForecast(
        point: point,
        standardErrors: stdErr,
        ciLower: ciLower,
        ciUpper: ciUpper,
        confidenceLevel: confidenceLevel,
        order: order,
        ar: fit.ar,
        ma: fit.ma,
        seasonalAr: fit.seasonalAr,
        seasonalMa: fit.seasonalMa,
        mean: fit.mean,
        sigma2: fit.sigma2,
        logLikelihood: fit.logLikelihood,
        aic: fit.aic,
        retval: fit.retval,
        attempts: List<ArimaAttempt>.unmodifiable(attempts),
      );
    } on TseriesNumericException catch (e) {
      // The wrapper attaches ctsa's status to the exception precisely so a
      // rejected rung is still measurable rather than being lost here.
      attempts.add(
        ArimaAttempt(
          order: order,
          accepted: false,
          retval: e.retval,
          note: 'fit did not converge to a finite model',
        ),
      );
    } on ArgumentError catch (e) {
      // Only order-specific rejections can reach here (finiteness and order
      // validity were checked above): the native size floor, or a
      // [ModelTooLargeError] for an order whose likelihood state is infeasible.
      // Treat them as a ladder skip rather than a hard error.
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
    'no valid ARIMA forecast for this series; all ${orders.length} candidate '
    'order(s) failed: ${attempts.join('; ')}',
  );
}
