/// Evaluation of one `autoArima` candidate: the pre-fit length rule, the fit
/// through the single fit path, and the rejection rules (spec §4–§6).
///
/// Internal to the package (not exported from `lib/tseries.dart`); tests use
/// it to evaluate reference orders exactly as the search does.
library;

import 'dart:math' as math;
import 'dart:typed_data';

import '../arima.dart';
import '../auto_arima.dart';
import '../information_criterion.dart';
import '../polynomial_roots.dart';
import '../sarimax.dart';

/// A point of the search space at fixed (d, D): the ARMA orders and whether
/// the model carries a constant (mean when d + D = 0, drift when d + D = 1).
typedef CandidateKey = ({int p, int q, int sp, int sq, bool c});

/// What evaluating one candidate produced. [fitCalled] tells whether the fit
/// path was invoked (it counts towards the `maxModels` budget).
class CandidateOutcome {
  CandidateOutcome({
    required this.order,
    required this.constant,
    required this.verdict,
    required this.parameterCount,
    required this.effectiveLength,
    required this.fitCalled,
    this.ctsaStatus,
    this.logLikelihood,
    this.criterion,
    this.minArRootModulus,
    this.minMaRootModulus,
    this.fit,
    this.detail,
  });

  final ArimaOrder order;
  final bool constant;
  final CandidateVerdict verdict;
  final int parameterCount;
  final int effectiveLength;
  final bool fitCalled;
  final int? ctsaStatus;
  final double? logLikelihood;
  final double? criterion;
  final double? minArRootModulus;
  final double? minMaRootModulus;
  final SarimaxFitResult? fit;
  final String? detail;
}

/// The order a candidate is fitted with. A candidate without seasonal terms
/// (P = D = Q = 0) is fitted as a plain ARIMA (period 0).
ArimaOrder candidateOrder(
  CandidateKey key, {
  required int d,
  required int seasonalD,
  required int period,
}) {
  final seasonal = key.sp + key.sq + seasonalD > 0;
  return ArimaOrder(
    p: key.p,
    d: d,
    q: key.q,
    seasonalP: key.sp,
    seasonalD: seasonalD,
    seasonalQ: key.sq,
    seasonalPeriod: seasonal ? period : 0,
  );
}

/// Effective length T = N − d − s·D (spec §4).
int effectiveLength(
  int n, {
  required int d,
  required int seasonalD,
  required int period,
}) => n - d - (seasonalD > 0 ? period * seasonalD : 0);

/// K = p + q + P + Q + [c] + 1 (σ² counted; no regressors in v1; spec §4).
int parameterCount(CandidateKey key) =>
    key.p + key.q + key.sp + key.sq + (key.c ? 1 : 0) + 1;

/// The pre-fit length rule (spec §5, rule 1), or null when the candidate may
/// be fitted. Returns the reason otherwise.
String? lengthRuleViolation(
  CandidateKey key, {
  required int t,
  required int period,
}) {
  final k = parameterCount(key);
  if (t - k - 1 < 1) {
    return 'T − K − 1 = ${t - k - 1} < 1 (T = $t, K = $k)';
  }
  final orders = key.p + key.q + key.sp + key.sq;
  final floor = math.max(10 * orders, 8);
  if (t < floor) {
    return 'T = $t < max(10·(p+q+P+Q), 8) = $floor';
  }
  final span = math.max(key.p + period * key.sp, key.q + period * key.sq);
  if (!(t > span + k)) {
    return 'T = $t ≤ max(p + s·P, q + s·Q) + K = ${span + k}';
  }
  return null;
}

/// Evaluates [key] on [y] at fixed ([d], [seasonalD]) — spec §5 rules 1–5 in
/// order, the first that fires being the verdict.
CandidateOutcome evaluateCandidate(
  Float64List y,
  CandidateKey key, {
  required int d,
  required int seasonalD,
  required int period,
  required InformationCriterion criterion,
  required int horizon,
  required double rootMargin,
}) {
  final order = candidateOrder(key, d: d, seasonalD: seasonalD, period: period);
  final t = effectiveLength(
    y.length,
    d: d,
    seasonalD: seasonalD,
    period: period,
  );
  final k = parameterCount(key);

  CandidateOutcome outcome(
    CandidateVerdict verdict, {
    required bool fitCalled,
    int? status,
    double? logLik,
    double? crit,
    double? minAr,
    double? minMa,
    SarimaxFitResult? fit,
    String? detail,
  }) => CandidateOutcome(
    order: order,
    constant: key.c,
    verdict: verdict,
    parameterCount: k,
    effectiveLength: t,
    fitCalled: fitCalled,
    ctsaStatus: status,
    logLikelihood: logLik,
    criterion: crit,
    minArRootModulus: minAr,
    minMaRootModulus: minMa,
    fit: fit,
    detail: detail,
  );

  // Rule 1: pre-fit length.
  final short = lengthRuleViolation(key, t: t, period: period);
  if (short != null) {
    return outcome(
      CandidateVerdict.rejectedTooShort,
      fitCalled: false,
      detail: short,
    );
  }

  // The single fit path (spec §6): fitSarimax without regressors; the
  // constant is the mean when d + D = 0 and the drift when d + D = 1.
  final SarimaxFitResult fit;
  try {
    fit = fitSarimax(
      y,
      const {},
      order: order,
      horizon: horizon,
      includeMean: key.c && d + seasonalD == 0,
      includeDrift: key.c && d + seasonalD == 1,
    );
  } on ModelTooLargeError catch (e) {
    // Rule 2.
    return outcome(
      CandidateVerdict.rejectedTooLarge,
      fitCalled: true,
      detail: '${e.runtimeType}: ${e.message}',
    );
  } on SeriesTooShortError catch (e) {
    // Not expected: rule 1 is meant to be stricter than every length guard
    // of the fit path. Recorded loudly, not swallowed (spec §5, ред. 3).
    return outcome(
      CandidateVerdict.rejectedTooShort,
      fitCalled: true,
      detail: 'LENGTH-RULE DEFECT: ${e.runtimeType}: ${e.message}',
    );
  } on TseriesNumericException catch (e) {
    final status = e.retval;
    if (status != null && status != ArimaFitStatus.probableSuccess.code) {
      return outcome(
        CandidateVerdict.rejectedStatus,
        fitCalled: true,
        status: status,
        detail: '${e.runtimeType}: ${e.message}',
      );
    }
    return outcome(
      status == null
          ? CandidateVerdict.rejectedFitError
          : CandidateVerdict.rejectedDegenerate,
      fitCalled: true,
      status: status,
      detail: '${e.runtimeType}: ${e.message}',
    );
  } on ArgumentError catch (e) {
    // The shim's own length/size refusals are plain ArgumentErrors. Rule 1 is
    // meant to make them unreachable; if one fires it is a defect of rule 1
    // (spec §5, ред. 3), so it is recorded with a marker, not hidden.
    return outcome(
      CandidateVerdict.rejectedFitError,
      fitCalled: true,
      detail:
          'LENGTH-RULE DEFECT? ${e.runtimeType} from the fit path: '
          '${e.message}',
    );
  }

  final coefficients = [
    ...fit.ar,
    ...fit.ma,
    ...fit.seasonalAr,
    ...fit.seasonalMa,
    if (fit.intercept != null) fit.intercept!,
    for (final r in fit.regressors)
      if (r.coefficient != null) r.coefficient!,
  ];
  final coefficientsFinite = coefficients.every((v) => v.isFinite);

  double? minAr;
  double? minMa;
  if (coefficientsFinite) {
    minAr = _combinedMinModulus(fit.ar, fit.seasonalAr, period);
    // ctsa's MA convention is 1 − θB, so a_j = θ_j (see polynomial_roots).
    minMa = _combinedMinModulus(fit.ma, fit.seasonalMa, period);
  }
  final logLik = fit.logLikelihood;
  final crit = logLik != null && logLik.isFinite
      ? criterion.evaluate(
          logLikelihood: logLik,
          parameters: k,
          observations: t,
        )
      : null;

  CandidateOutcome fitted(CandidateVerdict v, {String? detail}) => outcome(
    v,
    fitCalled: true,
    status: fit.retval,
    logLik: logLik,
    crit: crit,
    minAr: minAr,
    minMa: minMa,
    fit: fit,
    detail: detail,
  );

  // Rule 3: status.
  if (fit.retval != ArimaFitStatus.probableSuccess.code) {
    return fitted(CandidateVerdict.rejectedStatus);
  }
  if (logLik == null) return fitted(CandidateVerdict.rejectedNoLogLik);

  // Rule 4: degenerate fit.
  final degenerate = <String>[
    if (!(fit.sigma2 > 0) || !fit.sigma2.isFinite) 'sigma2 = ${fit.sigma2}',
    if (!logLik.isFinite) 'logLik = $logLik',
    if (!coefficientsFinite) 'non-finite coefficient',
    if (!fit.forecast.every((v) => v.isFinite)) 'non-finite forecast',
    if (!fit.standardErrors.every((v) => v.isFinite)) 'non-finite forecast SE',
  ];
  if (degenerate.isNotEmpty) {
    return fitted(
      CandidateVerdict.rejectedDegenerate,
      detail: degenerate.join(', '),
    );
  }

  // Rule 5: roots near the unit circle, checked per factor with the
  // step-down recursion; the seasonal factor in w = z^s uses margin^s.
  final seasonalMargin = math.pow(rootMargin, period).toDouble();
  final failing = <String>[
    if (!allRootsOutside(fit.ar, radius: rootMargin)) 'AR',
    if (!allRootsOutside(fit.seasonalAr, radius: seasonalMargin)) 'seasonal AR',
    if (!allRootsOutside(fit.ma, radius: rootMargin)) 'MA',
    if (!allRootsOutside(fit.seasonalMa, radius: seasonalMargin)) 'seasonal MA',
  ];
  if (failing.isNotEmpty) {
    return fitted(
      CandidateVerdict.rejectedRootsNearUnit,
      detail: 'root with modulus < $rootMargin in: ${failing.join(', ')}',
    );
  }
  return fitted(CandidateVerdict.accepted);
}

/// Smallest root modulus in z of `a(z)·A(z^s)`: the roots of the product are
/// those of the factors, and a root w of A gives |z| = |w|^(1/s).
double? _combinedMinModulus(
  List<double> nonSeasonal,
  List<double> seasonal,
  int period,
) {
  final a = minRootModulus(nonSeasonal);
  final sm = minRootModulus(seasonal);
  final b = sm == null ? null : math.pow(sm, 1 / period).toDouble();
  if (a == null) return b;
  if (b == null) return a;
  return math.min(a, b);
}
