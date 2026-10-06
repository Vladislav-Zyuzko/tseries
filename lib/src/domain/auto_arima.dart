/// Automatic ARIMA order selection (pure Dart over the native fixed-order fit).
///
/// The method is that of Hyndman & Khandakar (2008), "Automatic time series
/// forecasting: the forecast package for R", *Journal of Statistical
/// Software* 27(3), §3.1–3.2, as specified for this package (spec
/// "autoArima", 2026-10-05, revisions 1–3), with two documented departures:
/// the seasonal difference is chosen by the strength of seasonality of a
/// classical decomposition (FPP3 §9.1/§4.3) instead of a Canova–Hansen test,
/// and every candidate is fitted by exact maximum likelihood through the one
/// fit path [fitSarimax] (no CSS approximation during the search).
///
/// 1. **Differencing** (HK08 §3.1): first D — by the strength of seasonality
///    F_S of the classical additive decomposition, with a threshold stepped
///    by the number of full cycles k = ⌊N/s⌋: no automatic D below 5 cycles,
///    F_S ≥ 0.70 for 5–9 cycles, F_S ≥ 0.35 from 10 cycles (spec §2.2,
///    §15.1) —, then d by successive
///    KPSS tests of level stationarity at the 5 % level ([kpssLevelTest],
///    KPSS 1992; lag rule [KpssLag.short] by default, re-applied to the T of
///    each step) on the (seasonally differenced) series (spec §2.1, §14.1).
///    A working series that is constant after differencing is refused before
///    the search (spec §2.3). d and D
///    are fixed before the order search: likelihoods of differently
///    differenced series are not comparable (HK08 §3.1, spec §2.0).
/// 2. **Search** over (p, q, P, Q) and the constant at fixed (d, D): the
///    stepwise procedure of HK08 §3.2 with the neighbourhood of FPP3 §9.8
///    (up to 17 neighbours: p and/or q ±1, P and/or Q ±1, the constant
///    toggled; spec §3.2, §16.3), or an exhaustive grid
///    (spec §3.3).
/// 3. **Criterion** AIC / AICc / BIC computed from the fit's log-likelihood
///    (FPP3 §9.6, Hurvich & Tsai 1989; spec §4) — never ctsa's own AIC.
/// 4. **Rejection** of candidates (spec §5): too short for the order, too
///    large for the native state, a fit status other than success, no
///    log-likelihood, a degenerate fit, or an AR/MA root within 1.01 of the
///    unit circle (HK08 §3.2), checked without root finding by the
///    step-down recursion ([allRootsOutside]).
library;

import 'dart:math' as math;

import 'arima.dart';
import 'decomposition.dart';
import 'information_criterion.dart';
import 'internal/arima_support.dart';
import 'internal/auto_arima_candidate.dart';
import 'kpss.dart';
import 'polynomial_roots.dart';
import 'sarimax.dart';

/// F_S (classical decomposition) at or above which one seasonal difference
/// is taken when the series has at least [autoArimaSeasonalManyCycles] full
/// cycles. FPP3 §9.1 uses 0.64 on an STL decomposition; the classical
/// decomposition gives systematically smaller F_S (spec §2.2, §15.1).
/// Provisional.
const double autoArimaSeasonalStrengthThreshold = 0.35;

/// F_S threshold for series with [autoArimaSeasonalMinCycles] to
/// [autoArimaSeasonalManyCycles] − 1 full cycles: only strong seasonality
/// gives D = 1 there (an extra seasonal difference costs more than a missed
/// one, and seasonal AR/MA terms stay in the search). Spec §15.1.
/// Provisional.
const double autoArimaSeasonalStrengthThresholdFewCycles = 0.70;

/// Below this many full cycles k = ⌊N/s⌋ no automatic seasonal difference is
/// taken (D = 0, source `skippedShort`, F_S not computed): on 2–4 cycles F_S
/// of white noise exceeds 0.35 in 73–95 % of series (spec §15.1). Pass
/// `seasonalD: 1` when the seasonality is known. Provisional.
const int autoArimaSeasonalMinCycles = 5;

/// From this many full cycles on, [autoArimaSeasonalStrengthThreshold]
/// applies; below it, [autoArimaSeasonalStrengthThresholdFewCycles].
/// Provisional.
const int autoArimaSeasonalManyCycles = 10;

/// A neighbour replaces the current model only if its criterion is lower by
/// more than this; also the tie width of the final tie-breaks (spec §3.2,
/// §3.5). Our choice, not taken from the literature; provisional.
const double autoArimaCriterionEpsilon = 1e-6;

/// Default bound on p + q + P + Q of the exhaustive grid (spec §3.3).
/// Our choice, not taken from the literature; provisional.
const int autoArimaDefaultMaxOrderSum = 5;

/// Below this many observations no KPSS test is run and differencing stops
/// (spec §2.1). Our choice, not taken from the literature; provisional.
const int autoArimaKpssMinLength = 12;

/// The order search of [autoArima].
enum ArimaSearch {
  /// The stepwise procedure of Hyndman & Khandakar (2008, §3.2).
  stepwise,

  /// Every (p, q, P, Q, constant) within the bounds and
  /// p + q + P + Q ≤ maxOrderSum, in lexicographic order (spec §3.3).
  exhaustive,
}

/// Where a differencing order of [autoArima] came from: why the procedure
/// that chose it stopped.
enum DifferencingSource {
  /// Given by the caller; no test was run.
  user,

  /// A test decided it: KPSS did not reject stationarity (d), or F_S was
  /// compared with the threshold (D).
  test,

  /// The series was too short to test (fewer than [autoArimaKpssMinLength]
  /// observations for KPSS, fewer than [autoArimaSeasonalMinCycles] full
  /// seasonal cycles for F_S — spec §15.1); differencing
  /// stopped there.
  skippedShort,

  /// The series (after the differences taken so far) is constant, or has no
  /// variation around its trend (`max|x − x̄| ≤ 1e-10·max|y|`, y the input):
  /// the statistic is undefined; differencing stopped there.
  constantSeries,

  /// The maximum was reached while the test still asked for more: for d, the
  /// KPSS test at maxD (computed and reported) still rejects stationarity;
  /// for D, maxSeasonalD = 0. (Spec §1, §2.1, revision 5.)
  limitedByMax,

  /// D only: the series is not seasonal (period 1), so D = 0 by definition.
  /// Our addition to the spec's list.
  notSeasonal,
}

/// How [autoArima] chose d and D (spec §1, §2).
class DifferencingReport {
  const DifferencingReport({
    required this.d,
    required this.dSource,
    required this.kpss,
    required this.seasonalD,
    required this.seasonalDSource,
    required this.seasonalStrength,
    required this.seasonalCycles,
    required this.seasonalThreshold,
    required this.period,
    required this.kpssLag,
  });

  /// The number of ordinary differences d.
  final int d;

  /// Where [d] came from.
  final DifferencingSource dSource;

  /// Every KPSS test run, in order (the first on the base series, each later
  /// one on the previous series differenced once more, the last possibly at
  /// maxD); each records its T, lag rule and actual lag. Empty when d was
  /// given or nothing was tested.
  final List<KpssResult> kpss;

  /// The KPSS lag rule (applied to the T of each step).
  final KpssLag kpssLag;

  /// The number of seasonal differences D (0 or 1).
  final int seasonalD;

  /// Where [seasonalD] came from.
  final DifferencingSource seasonalDSource;

  /// F_S of the classical decomposition of the series, or null when it was
  /// not computed (D given, period 1, fewer than
  /// [autoArimaSeasonalMinCycles] cycles, maxSeasonalD = 0, constant).
  final double? seasonalStrength;

  /// Number of full seasonal cycles k = ⌊N/s⌋ of the input, or null for a
  /// non-seasonal series (period 1).
  final int? seasonalCycles;

  /// The threshold F_S was compared with (0.70 for 5–9 cycles, 0.35 from 10
  /// cycles), or null when no comparison was made.
  final double? seasonalThreshold;

  /// The seasonal period s (1 for a non-seasonal series).
  final int period;

  /// A stable one-line-per-fact text form (part of the deterministic journal).
  String toJournal() {
    final b = StringBuffer()
      ..writeln(
        'D=$seasonalD source=${seasonalDSource.name} period=$period '
        'cycles=${seasonalCycles ?? '-'} F_S=${seasonalStrength ?? '-'} '
        'threshold=${seasonalThreshold ?? '-'}',
      )
      ..writeln('d=$d source=${dSource.name} lagRule=${kpssLag.description}');
    for (var i = 0; i < kpss.length; i++) {
      final k = kpss[i];
      b.writeln(
        'kpss[$i] T=${k.length} l=${k.lag} (${k.lagRule.description}) '
        'eta=${k.statistic} '
        's2=${k.longRunVariance} crit=${k.criticalValue} '
        'p=${k.pValue.label} '
        '${k.rejectsStationarity ? 'difference' : 'stop'}',
      );
    }
    return b.toString();
  }

  @override
  String toString() => toJournal();
}

/// The fate of one candidate of the order search (spec §1, §5). The first
/// rule of §5 that fires is the verdict.
enum CandidateVerdict {
  /// Fitted, passed every rule; competes on the criterion.
  accepted,

  /// ctsa's fit status was not "probable success" ([CandidateRecord.ctsaStatus]
  /// holds the code: 4 = iteration limit, 10/12 = non-stationary CSS start
  /// with no MLE step, 15 = non-finite, ...).
  rejectedStatus,

  /// The fit reported no log-likelihood.
  rejectedNoLogLik,

  /// σ² not positive/finite, or a non-finite log-likelihood, coefficient,
  /// forecast or forecast standard error — also a numeric failure reported
  /// by the fit with a successful status.
  rejectedDegenerate,

  /// An AR or MA factor has a root of modulus below the margin (1.01 by
  /// default, spec §5 revision 9; HK08 §3.2 use 1.001): a non-seasonal
  /// factor in z, a seasonal factor in w = z^s against margin^s.
  rejectedRootsNearUnit,

  /// The series is too short for the order (spec §5 rule 1); not fitted.
  rejectedTooShort,

  /// The order's exact-likelihood state is infeasible ([ModelTooLargeError]).
  rejectedTooLarge,

  /// The fit path raised an error not covered above (see
  /// [CandidateRecord.detail]).
  rejectedFitError,

  /// Already evaluated earlier in the search; not refitted
  /// ([CandidateRecord.duplicateOf] is the earlier record).
  cachedDuplicate,
}

/// The phase of the search a [CandidateRecord] belongs to.
enum CandidateStage {
  /// The starting models of the stepwise search (HK08 §3.2).
  initial,

  /// A neighbour of the current model in step [CandidateRecord.step].
  step,

  /// A point of the exhaustive grid.
  exhaustive,
}

/// One entry of the search journal ([AutoArimaResult.trace]), in the order of
/// evaluation.
class CandidateRecord {
  const CandidateRecord({
    required this.index,
    required this.stage,
    required this.step,
    required this.order,
    required this.constant,
    required this.ctsaStatus,
    required this.logLikelihood,
    required this.parameterCount,
    required this.effectiveLength,
    required this.criterion,
    required this.minArRootModulus,
    required this.minMaRootModulus,
    required this.verdict,
    required this.becameCurrent,
    required this.duplicateOf,
    required this.detail,
    required this.fit,
  });

  /// Position in the journal (0-based).
  final int index;

  /// Search phase.
  final CandidateStage stage;

  /// For [CandidateStage.step]: the step number (1-based: the k-th
  /// neighbourhood scan); 0 otherwise.
  final int step;

  /// The full order `(p, d, q)(P, D, Q)[s]`. A candidate without seasonal
  /// terms and D = 0 carries period 0.
  final ArimaOrder order;

  /// Whether the candidate has a constant: the mean when d + D = 0, the drift
  /// when d + D = 1.
  final bool constant;

  /// ctsa's fit status (`retval`), or null when the fit never ran or failed
  /// before reporting one.
  final int? ctsaStatus;

  /// The fit's log-likelihood, when it reported one.
  final double? logLikelihood;

  /// K = p + q + P + Q + [constant] + 1 (σ² counted).
  final int parameterCount;

  /// T = N − d − s·D, the number of observations in the likelihood.
  final int effectiveLength;

  /// The selected criterion, when a finite log-likelihood exists (computed
  /// for rejected candidates too, for transparency; only accepted ones
  /// compete).
  final double? criterion;

  /// Smallest root modulus of φ(z)Φ(z^s), or null when there are no AR terms
  /// or no finite coefficients. Exact for factor degree ≤ 2, bisection
  /// otherwise.
  final double? minArRootModulus;

  /// As [minArRootModulus], for θ(z)Θ(z^s).
  final double? minMaRootModulus;

  /// The verdict (first rule of spec §5 that fired).
  final CandidateVerdict verdict;

  /// Whether this candidate became the current model at this point (stepwise)
  /// or is the selected model (exhaustive).
  final bool becameCurrent;

  /// For [CandidateVerdict.cachedDuplicate]: the index of the earlier record
  /// of the same candidate.
  final int? duplicateOf;

  /// Why a candidate was rejected, in words (null when accepted or cached).
  final String? detail;

  /// The fit itself, when the fit path returned one (also for some rejected
  /// candidates). σ² is the profile MLE S/T.
  final SarimaxFitResult? fit;

  /// A stable text form of this record (part of the deterministic journal).
  String toJournalLine() {
    final stageText = switch (stage) {
      CandidateStage.initial => 'initial',
      CandidateStage.step => 'step$step',
      CandidateStage.exhaustive => 'exhaustive',
    };
    return '#$index $stageText $order c=${constant ? 1 : 0} '
        '${verdict.name}'
        '${duplicateOf != null ? ' of=#$duplicateOf' : ''} '
        'status=${ctsaStatus ?? '-'} ll=${logLikelihood ?? '-'} '
        'K=$parameterCount T=$effectiveLength crit=${criterion ?? '-'} '
        'arRoot=${minArRootModulus ?? '-'} maRoot=${minMaRootModulus ?? '-'}'
        '${becameCurrent ? ' current' : ''}'
        '${detail != null ? ' | $detail' : ''}';
  }

  @override
  String toString() => toJournalLine();
}

/// The selected model of [autoArima], its forecast, and the full journal.
class AutoArimaResult {
  const AutoArimaResult({
    required this.order,
    required this.constant,
    required this.fit,
    required this.criterion,
    required this.criterionValue,
    required this.search,
    required this.forecast,
    required this.lower,
    required this.upper,
    required this.level,
    required this.differencing,
    required this.trace,
    required this.fitsEvaluated,
    required this.truncated,
  });

  /// The selected order `(p, d, q)(P, D, Q)[s]`.
  final ArimaOrder order;

  /// Whether the selected model has a constant (mean if d + D = 0, drift if
  /// d + D = 1; see [SarimaxFitResult.intercept] and the regressor named
  /// [sarimaxDriftColumn] in [fit]).
  final bool constant;

  /// The fit of the selected model — the same fit that was compared during
  /// the search (no refit, spec §3.4). σ² is the profile MLE S/T, so
  /// intervals are slightly narrower than with a degrees-of-freedom
  /// corrected σ² (spec §4).
  final SarimaxFitResult fit;

  /// The criterion that was minimised.
  final InformationCriterion criterion;

  /// Its value for the selected model.
  final double criterionValue;

  /// The search that was run.
  final ArimaSearch search;

  /// Point forecast (empty when horizon was 0).
  final List<double> forecast;

  /// `forecast − z·SE` for the two-sided normal z of [level].
  final List<double> lower;

  /// `forecast + z·SE`.
  final List<double> upper;

  /// The level of [lower]/[upper].
  final double level;

  /// How d and D were chosen.
  final DifferencingReport differencing;

  /// Every candidate in evaluation order.
  final List<CandidateRecord> trace;

  /// Number of calls of the fit path (the unit of the maxModels budget).
  final int fitsEvaluated;

  /// True when the search stopped at maxModels: the result is the best model
  /// found until then.
  final bool truncated;

  /// The differencing report and journal as stable text — byte-identical
  /// across runs of the same build (gate E).
  String toJournal() {
    final b = StringBuffer()
      ..write(differencing.toJournal())
      ..writeln(
        'search=${search.name} criterion=${criterion.name} '
        'fits=$fitsEvaluated truncated=$truncated',
      );
    for (final r in trace) {
      b.writeln(r.toJournalLine());
    }
    b.writeln('selected $order c=${constant ? 1 : 0} crit=$criterionValue');
    return b.toString();
  }
}

/// Why [autoArima] returned no model (spec §1, revision 5).
enum AutoArimaNoModelReason {
  /// The working series (after the chosen d and D) is constant: every model
  /// would be degenerate or spurious (spec §2.3). Raised before the search;
  /// the journal is empty.
  constantSeries,

  /// Stepwise search: none of the starting models was accepted (spec §3.2).
  noStartModelAccepted,

  /// No candidate was accepted (spec §5 rule 6).
  noCandidateAccepted,
}

/// No model was selected by [autoArima]. The journal is carried in the
/// fields, not only in the message (spec §1, §5 rule 6): there is no silent
/// fallback model — choosing one is the caller's decision.
final class AutoArimaNoModelException implements Exception {
  const AutoArimaNoModelException({
    required this.reason,
    required this.message,
    required this.differencing,
    required this.trace,
    required this.fitsEvaluated,
    required this.truncated,
  });

  /// Why no model was selected.
  final AutoArimaNoModelReason reason;

  /// Summary of why nothing was accepted.
  final String message;

  /// How d and D were chosen.
  final DifferencingReport differencing;

  /// Every candidate evaluated, with its verdict.
  final List<CandidateRecord> trace;

  /// Number of calls of the fit path.
  final int fitsEvaluated;

  /// Whether the maxModels budget was exhausted.
  final bool truncated;

  @override
  String toString() => 'AutoArimaNoModelException(${reason.name}): $message';
}

/// Chooses a (seasonal) ARIMA order for [series] automatically and fits it.
///
/// See the library documentation for the method and its sources. Parameters:
///
/// * [period] — the seasonal period s; 1 for a non-seasonal series.
/// * [d], [seasonalD] — fix the differencing orders; null chooses them by
///   test (KPSS for d, F_S for D). [seasonalD] is only meaningful when
///   [period] > 1.
/// * [maxD], [maxSeasonalD] — bounds of the automatic choice (maxSeasonalD
///   is 0 or 1).
/// * [maxP], [maxQ], [maxSeasonalP], [maxSeasonalQ] — bounds of the search
///   (HK08 §3.2 defaults 5, 5, 2, 2).
/// * [maxOrderSum] — exhaustive search only: p + q + P + Q ≤ this
///   (default [autoArimaDefaultMaxOrderSum]).
/// * [criterion] — AICc by default.
/// * [allowConstant] — consider a constant (mean if d + D = 0, drift if
///   d + D = 1); never with d + D ≥ 2.
/// * [horizon], [level] — forecast of the selected model and its band.
/// * [kpssAlpha] — level of the KPSS tests; one of the Table 1 levels
///   (see [kpssLevelCriticalValues]). The decision compares η with the table
///   value also at 1 %, where R's `ndiffs` never differences (its p-value is
///   bounded below) — a deliberate difference.
/// * [kpssLag] — the KPSS lag rule, re-applied to the T of every step
///   (default [KpssLag.short], which reproduces R's `ndiffs`; [KpssLag.l4] is
///   the KPSS 1992 rule).
/// * [rootMargin] — candidates with an AR or MA root of modulus below this
///   are rejected (default 1.01, spec §5 revision 9; HK08 §3.2 use 1.001).
///   A seasonal factor Φ(w), w = z^s, is checked against rootMargin^s.
/// * [maxModels] — budget of fits; when reached the search stops and the
///   best model so far is returned with [AutoArimaResult.truncated] set.
///   There is no time limit (it would make the choice machine-dependent).
///
/// Synchronous and deterministic: the same input gives the same journal.
/// A search can take seconds for seasonal series; run it on a background
/// isolate (`Isolate.run`) where that matters.
///
/// Throws [ArgumentError] for invalid arguments (including a non-finite or
/// empty [series]) and [AutoArimaNoModelException] when the working series is
/// constant or no candidate is accepted (see [AutoArimaNoModelReason]).
AutoArimaResult autoArima(
  List<num> series, {
  int period = 1,
  int? d,
  int? seasonalD,
  int maxD = 2,
  int maxSeasonalD = 1,
  int maxP = 5,
  int maxQ = 5,
  int maxSeasonalP = 2,
  int maxSeasonalQ = 2,
  int? maxOrderSum,
  InformationCriterion criterion = InformationCriterion.aicc,
  ArimaSearch search = ArimaSearch.stepwise,
  bool allowConstant = true,
  int horizon = 0,
  double level = 0.95,
  double kpssAlpha = 0.05,
  KpssLag kpssLag = KpssLag.short,
  double rootMargin = 1.01,
  int maxModels = 100,
}) {
  // ---- Arguments ---------------------------------------------------------
  final y = toFiniteFloat64(series);
  if (y.isEmpty) {
    throw ArgumentError.value(series, 'series', 'must not be empty');
  }
  if (period < 1) {
    throw ArgumentError.value(period, 'period', 'must be at least 1');
  }
  if (d != null && d < 0) {
    throw ArgumentError.value(d, 'd', 'must be non-negative');
  }
  if (seasonalD != null) {
    if (seasonalD < 0) {
      throw ArgumentError.value(seasonalD, 'seasonalD', 'must be non-negative');
    }
    if (period == 1 && seasonalD != 0) {
      throw ArgumentError.value(seasonalD, 'seasonalD', 'requires period > 1');
    }
  }
  for (final (name, value) in [
    ('maxD', maxD),
    ('maxP', maxP),
    ('maxQ', maxQ),
    ('maxSeasonalP', maxSeasonalP),
    ('maxSeasonalQ', maxSeasonalQ),
    ('horizon', horizon),
  ]) {
    if (value < 0) {
      throw ArgumentError.value(value, name, 'must be non-negative');
    }
  }
  if (maxSeasonalD < 0 || maxSeasonalD > 1) {
    throw ArgumentError.value(maxSeasonalD, 'maxSeasonalD', 'must be 0 or 1');
  }
  if (maxOrderSum != null && maxOrderSum < 0) {
    throw ArgumentError.value(maxOrderSum, 'maxOrderSum', 'must be ≥ 0');
  }
  validateConfidenceLevel(level);
  if (!kpssLevelCriticalValues.containsKey(kpssAlpha)) {
    throw ArgumentError.value(
      kpssAlpha,
      'kpssAlpha',
      'must be one of ${kpssLevelCriticalValues.keys.join(', ')}',
    );
  }
  if (!(rootMargin >= 1) || !rootMargin.isFinite) {
    throw ArgumentError.value(rootMargin, 'rootMargin', 'must be ≥ 1');
  }
  if (maxModels < 1) {
    throw ArgumentError.value(maxModels, 'maxModels', 'must be at least 1');
  }

  // ---- Differencing (spec §2) -------------------------------------------
  final differencing = _chooseDifferencing(
    y,
    period: period,
    userD: d,
    userSeasonalD: seasonalD,
    maxD: maxD,
    maxSeasonalD: maxSeasonalD,
    kpssAlpha: kpssAlpha,
    kpssLag: kpssLag,
  );
  final dd = differencing.d;
  final sD = differencing.seasonalD;
  final seasonal = period > 1;
  final constantAllowed = allowConstant && dd + sD < 2;

  // ---- Degenerate working series (spec §2.3) ----------------------------
  var working = List<double>.of(y);
  for (var i = 0; i < sD; i++) {
    working = _seasonalDiff(working, period);
  }
  for (var i = 0; i < dd; i++) {
    working = _diff(working);
  }
  if (isConstantSeries(working, maxAbs(y))) {
    throw AutoArimaNoModelException(
      reason: AutoArimaNoModelReason.constantSeries,
      message:
          'the series is constant after d=$dd, D=$sD differences '
          '(max|x − mean| ≤ $constantSeriesTolerance·max|y|): no ARIMA model '
          'is meaningful',
      differencing: differencing,
      trace: const [],
      fitsEvaluated: 0,
      truncated: false,
    );
  }

  // ---- Search state ------------------------------------------------------
  final records = <_Record>[];
  final memo = <CandidateKey, int>{};
  var fits = 0;
  var truncated = false;

  /// Evaluates [key] (or records a cache hit). Returns null when the budget
  /// is exhausted, which ends the search.
  _Record? visit(CandidateKey key, CandidateStage stage, int step) {
    final cached = memo[key];
    if (cached != null) {
      final orig = records[cached];
      final rec = _Record(
        index: records.length,
        stage: stage,
        step: step,
        key: key,
        outcome: orig.outcome,
        duplicateOf: cached,
      );
      records.add(rec);
      return rec;
    }
    final t = effectiveLength(y.length, d: dd, seasonalD: sD, period: period);
    final needsFit = lengthRuleViolation(key, t: t, period: period) == null;
    if (needsFit && fits >= maxModels) {
      truncated = true;
      return null;
    }
    final outcome = evaluateCandidate(
      y,
      key,
      d: dd,
      seasonalD: sD,
      period: period,
      criterion: criterion,
      horizon: horizon,
      rootMargin: rootMargin,
    );
    if (outcome.fitCalled) fits++;
    final rec = _Record(
      index: records.length,
      stage: stage,
      step: step,
      key: key,
      outcome: outcome,
    );
    records.add(rec);
    memo[key] = rec.index;
    return rec;
  }

  _Record? selected;
  var noStart = false;
  switch (search) {
    case ArimaSearch.stepwise:
      selected = _stepwise(
        visit,
        seasonal: seasonal,
        constantAllowed: constantAllowed,
        maxP: maxP,
        maxQ: maxQ,
        maxSP: seasonal ? maxSeasonalP : 0,
        maxSQ: seasonal ? maxSeasonalQ : 0,
        records: records,
      );
      noStart = !records.any((r) => r.isAccepted);
    case ArimaSearch.exhaustive:
      final sum = maxOrderSum ?? autoArimaDefaultMaxOrderSum;
      outer:
      for (var p = 0; p <= maxP; p++) {
        for (var q = 0; q <= maxQ; q++) {
          for (var sp = 0; sp <= (seasonal ? maxSeasonalP : 0); sp++) {
            for (var sq = 0; sq <= (seasonal ? maxSeasonalQ : 0); sq++) {
              if (p + q + sp + sq > sum) continue;
              for (final c in constantAllowed ? [false, true] : [false]) {
                final r = visit(
                  (p: p, q: q, sp: sp, sq: sq, c: c),
                  CandidateStage.exhaustive,
                  0,
                );
                if (r == null) break outer;
              }
            }
          }
        }
      }
      for (final r in records) {
        if (r.duplicateOf != null || !r.isAccepted) continue;
        if (selected == null || _better(r, selected)) selected = r;
      }
      selected?.becameCurrent = true;
  }

  final trace = List<CandidateRecord>.unmodifiable([
    for (final r in records) r.freeze(),
  ]);
  if (selected == null) {
    final counts = <CandidateVerdict, int>{};
    for (final r in records) {
      counts[r.outcome.verdict] = (counts[r.outcome.verdict] ?? 0) + 1;
    }
    throw AutoArimaNoModelException(
      reason: noStart
          ? AutoArimaNoModelReason.noStartModelAccepted
          : AutoArimaNoModelReason.noCandidateAccepted,
      message:
          '${noStart ? 'no starting model' : 'no candidate'} accepted at '
          'd=$dd, D=$sD after $fits fit(s)'
          '${truncated ? ' (budget maxModels=$maxModels exhausted)' : ''}; '
          'verdicts: '
          '${counts.entries.map((e) => '${e.key.name}=${e.value}').join(', ')}',
      differencing: differencing,
      trace: trace,
      fitsEvaluated: fits,
      truncated: truncated,
    );
  }

  final fit = selected.outcome.fit!;
  final z = normalTwoSidedZ(level);
  return AutoArimaResult(
    order: selected.outcome.order,
    constant: selected.key.c,
    fit: fit,
    criterion: criterion,
    criterionValue: selected.outcome.criterion!,
    search: search,
    forecast: fit.forecast,
    lower: List<double>.unmodifiable([
      for (var i = 0; i < fit.forecast.length; i++)
        fit.forecast[i] - z * fit.standardErrors[i],
    ]),
    upper: List<double>.unmodifiable([
      for (var i = 0; i < fit.forecast.length; i++)
        fit.forecast[i] + z * fit.standardErrors[i],
    ]),
    level: level,
    differencing: differencing,
    trace: trace,
    fitsEvaluated: fits,
    truncated: truncated,
  );
}

/// The stepwise search of HK08 §3.2, neighbourhood FPP3 §9.8 (spec §3.2,
/// §16.3). Returns the final current
/// model, or null when no starting model was accepted or the budget ran out
/// before one was.
_Record? _stepwise(
  _Record? Function(CandidateKey, CandidateStage, int) visit, {
  required bool seasonal,
  required bool constantAllowed,
  required int maxP,
  required int maxQ,
  required int maxSP,
  required int maxSQ,
  required List<_Record> records,
}) {
  CandidateKey clip(int p, int q, int sp, int sq, bool c) => (
    p: math.min(p, maxP),
    q: math.min(q, maxQ),
    sp: seasonal ? math.min(sp, maxSP) : 0,
    sq: seasonal ? math.min(sq, maxSQ) : 0,
    c: c,
  );

  // Step 1: starting models (HK08 §3.2; the extra model without a constant
  // is FPP3 §9.8). Duplicates after clipping are dropped.
  final c = constantAllowed;
  final starts = <CandidateKey>{
    clip(2, 2, 1, 1, c),
    clip(0, 0, 0, 0, c),
    clip(1, 0, 1, 0, c),
    clip(0, 1, 0, 1, c),
    if (c) clip(0, 0, 0, 0, false),
  };
  _Record? current;
  for (final key in starts) {
    final r = visit(key, CandidateStage.initial, 0);
    if (r == null) break;
    if (!r.isAccepted) continue;
    if (current == null || _better(r, current)) current = r;
  }
  if (current == null) return null;
  current.becameCurrent = true;

  // Step 2+: first improvement among at most 17 neighbours (FPP3 §9.8:
  // p and/or q ±1, P and/or Q ±1, the constant toggled), in the fixed order
  // of spec revision 9 (observed in R's stepwise traces): the seasonal block
  // first — P−1, Q−1, P+1, Q+1, (−,−), (−,+), (+,−), (+,+) —, then the same
  // for (p, q), then the constant. maxOrderSum does not apply here (it bounds
  // only the exhaustive grid).
  var step = 0;
  while (true) {
    step++;
    final k = current!.key;
    final neighbours = <CandidateKey>[
      if (seasonal) ...[
        (p: k.p, q: k.q, sp: k.sp - 1, sq: k.sq, c: k.c),
        (p: k.p, q: k.q, sp: k.sp, sq: k.sq - 1, c: k.c),
        (p: k.p, q: k.q, sp: k.sp + 1, sq: k.sq, c: k.c),
        (p: k.p, q: k.q, sp: k.sp, sq: k.sq + 1, c: k.c),
        (p: k.p, q: k.q, sp: k.sp - 1, sq: k.sq - 1, c: k.c),
        (p: k.p, q: k.q, sp: k.sp - 1, sq: k.sq + 1, c: k.c),
        (p: k.p, q: k.q, sp: k.sp + 1, sq: k.sq - 1, c: k.c),
        (p: k.p, q: k.q, sp: k.sp + 1, sq: k.sq + 1, c: k.c),
      ],
      (p: k.p - 1, q: k.q, sp: k.sp, sq: k.sq, c: k.c),
      (p: k.p, q: k.q - 1, sp: k.sp, sq: k.sq, c: k.c),
      (p: k.p + 1, q: k.q, sp: k.sp, sq: k.sq, c: k.c),
      (p: k.p, q: k.q + 1, sp: k.sp, sq: k.sq, c: k.c),
      (p: k.p - 1, q: k.q - 1, sp: k.sp, sq: k.sq, c: k.c),
      (p: k.p - 1, q: k.q + 1, sp: k.sp, sq: k.sq, c: k.c),
      (p: k.p + 1, q: k.q - 1, sp: k.sp, sq: k.sq, c: k.c),
      (p: k.p + 1, q: k.q + 1, sp: k.sp, sq: k.sq, c: k.c),
      if (constantAllowed) (p: k.p, q: k.q, sp: k.sp, sq: k.sq, c: !k.c),
    ];
    var improved = false;
    for (final n in neighbours) {
      if (n.p < 0 || n.p > maxP || n.q < 0 || n.q > maxQ) continue;
      if (n.sp < 0 || n.sp > maxSP || n.sq < 0 || n.sq > maxSQ) continue;
      final r = visit(n, CandidateStage.step, step);
      if (r == null) return current;
      final target = r.duplicateOf == null ? r : records[r.duplicateOf!];
      if (!target.isAccepted) continue;
      if (target.outcome.criterion! <
          current!.outcome.criterion! - autoArimaCriterionEpsilon) {
        r.becameCurrent = true;
        current = target;
        improved = true;
        break;
      }
    }
    if (!improved) return current;
  }
}

/// Tie-breaks of spec §3.5: lower criterion; within ε, smaller K; then the
/// lexicographically smaller (p, q, P, Q); then no constant.
bool _better(_Record a, _Record b) {
  final ca = a.outcome.criterion!;
  final cb = b.outcome.criterion!;
  if ((ca - cb).abs() > autoArimaCriterionEpsilon) return ca < cb;
  final ka = a.outcome.parameterCount;
  final kb = b.outcome.parameterCount;
  if (ka != kb) return ka < kb;
  final x = a.key;
  final y = b.key;
  for (final (u, v) in [(x.p, y.p), (x.q, y.q), (x.sp, y.sp), (x.sq, y.sq)]) {
    if (u != v) return u < v;
  }
  if (x.c != y.c) return !x.c;
  return false;
}

class _Record {
  _Record({
    required this.index,
    required this.stage,
    required this.step,
    required this.key,
    required this.outcome,
    this.duplicateOf,
  });

  final int index;
  final CandidateStage stage;
  final int step;
  final CandidateKey key;
  final CandidateOutcome outcome;
  final int? duplicateOf;
  bool becameCurrent = false;

  bool get isAccepted =>
      duplicateOf == null && outcome.verdict == CandidateVerdict.accepted;

  CandidateRecord freeze() {
    final o = outcome;
    final dup = duplicateOf != null;
    return CandidateRecord(
      index: index,
      stage: stage,
      step: step,
      order: o.order,
      constant: o.constant,
      ctsaStatus: o.ctsaStatus,
      logLikelihood: o.logLikelihood,
      parameterCount: o.parameterCount,
      effectiveLength: o.effectiveLength,
      criterion: o.criterion,
      minArRootModulus: o.minArRootModulus,
      minMaRootModulus: o.minMaRootModulus,
      verdict: dup ? CandidateVerdict.cachedDuplicate : o.verdict,
      becameCurrent: becameCurrent,
      duplicateOf: duplicateOf,
      detail: dup ? null : o.detail,
      fit: o.fit,
    );
  }
}

// ---- Differencing ----------------------------------------------------------

DifferencingReport _chooseDifferencing(
  List<double> y, {
  required int period,
  required int? userD,
  required int? userSeasonalD,
  required int maxD,
  required int maxSeasonalD,
  required double kpssAlpha,
  required KpssLag kpssLag,
}) {
  // D first (HK08 §3.1; spec §2.2).
  int sD;
  DifferencingSource sSource;
  double? strength;
  double? threshold;
  final int? cycles = period > 1 ? y.length ~/ period : null;
  if (period == 1) {
    sD = 0;
    sSource = DifferencingSource.notSeasonal;
  } else if (userSeasonalD != null) {
    sD = userSeasonalD;
    sSource = DifferencingSource.user;
  } else if (maxSeasonalD == 0) {
    sD = 0;
    sSource = DifferencingSource.limitedByMax;
  } else if (cycles! < autoArimaSeasonalMinCycles) {
    // Spec §15.1: F_S is indistinguishable from noise on 2–4 cycles.
    sD = 0;
    sSource = DifferencingSource.skippedShort;
  } else {
    final dec = classicalDecomposition(y, period: period);
    if (dec.isDegenerate) {
      sD = 0;
      sSource = DifferencingSource.constantSeries;
    } else {
      strength = dec.seasonalStrength;
      threshold = cycles < autoArimaSeasonalManyCycles
          ? autoArimaSeasonalStrengthThresholdFewCycles
          : autoArimaSeasonalStrengthThreshold;
      sD = strength >= threshold ? 1 : 0;
      sSource = DifferencingSource.test;
    }
  }

  var x = y;
  for (var i = 0; i < sD; i++) {
    x = _seasonalDiff(x, period);
  }

  // Then d by successive KPSS tests (spec §2.1). The lag rule is re-applied
  // to the T of each step; at maxD the test is still run and reported, and
  // limitedByMax is set only when it still rejects stationarity.
  final scale = maxAbs(y);
  final tests = <KpssResult>[];
  int d;
  DifferencingSource dSource;
  if (userD != null) {
    d = userD;
    dSource = DifferencingSource.user;
  } else {
    d = 0;
    final crit = kpssLevelCriticalValues[kpssAlpha]!;
    while (true) {
      if (x.length < autoArimaKpssMinLength) {
        dSource = DifferencingSource.skippedShort;
        break;
      }
      // A fixed lag cannot exceed T − 1; clamp (the actual lag is reported).
      final lag = math.min(kpssLag.lagFor(x.length), x.length - 1);
      final stat = kpssStatistic(x, lag, scale: scale);
      if (stat == null) {
        dSource = DifferencingSource.constantSeries;
        break;
      }
      final result = KpssResult(
        length: x.length,
        lagRule: kpssLag,
        lag: lag,
        statistic: stat.eta,
        longRunVariance: stat.s2,
        alpha: kpssAlpha,
        criticalValue: crit,
      );
      tests.add(result);
      if (!result.rejectsStationarity) {
        dSource = DifferencingSource.test;
        break;
      }
      if (d == maxD) {
        dSource = DifferencingSource.limitedByMax;
        break;
      }
      x = _diff(x);
      d++;
    }
  }

  return DifferencingReport(
    d: d,
    dSource: dSource,
    kpss: List<KpssResult>.unmodifiable(tests),
    seasonalD: sD,
    seasonalDSource: sSource,
    seasonalStrength: strength,
    seasonalCycles: cycles,
    seasonalThreshold: threshold,
    period: period,
    kpssLag: kpssLag,
  );
}

List<double> _diff(List<double> x) => [
  for (var i = 1; i < x.length; i++) x[i] - x[i - 1],
];

List<double> _seasonalDiff(List<double> x, int s) => [
  for (var i = s; i < x.length; i++) x[i] - x[i - s],
];
