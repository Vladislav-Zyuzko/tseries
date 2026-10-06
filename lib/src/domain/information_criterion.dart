/// Information criteria computed from a log-likelihood (pure Dart).
///
/// Sources: Hyndman & Athanasopoulos, *Forecasting: Principles and Practice*
/// (3rd ed.), §9.6; the small-sample correction AICc is Hurvich & Tsai (1989).
library;

import 'dart:math' as math;

/// An information criterion for comparing models fitted to the **same** data
/// (same series, same differencing, hence the same number of observations in
/// the likelihood). Lower is better.
///
/// With log-likelihood `logL`, `K` estimated parameters (σ² included) and `T`
/// observations in the likelihood (FPP3 §9.6, spec §4):
///
/// ```text
/// AIC  = −2·logL + 2K
/// AICc = AIC + 2K(K+1)/(T − K − 1)
/// BIC  = AIC + (ln T − 2)·K
/// ```
enum InformationCriterion {
  /// Akaike's information criterion.
  aic,

  /// AIC with the small-sample correction of Hurvich & Tsai (1989).
  aicc,

  /// Schwarz's Bayesian information criterion.
  bic;

  /// The criterion for a model with log-likelihood [logLikelihood],
  /// [parameters] estimated parameters (K, counting σ²) and [observations]
  /// observations in the likelihood (T, after differencing).
  ///
  /// Throws [ArgumentError] for a non-finite [logLikelihood], `K < 1`,
  /// `T < 1`, or — for [aicc] — `T − K − 1 ≤ 0`, where the correction is
  /// undefined.
  double evaluate({
    required double logLikelihood,
    required int parameters,
    required int observations,
  }) {
    if (!logLikelihood.isFinite) {
      throw ArgumentError.value(
        logLikelihood,
        'logLikelihood',
        'must be finite',
      );
    }
    final k = parameters;
    final t = observations;
    if (k < 1) {
      throw ArgumentError.value(k, 'parameters', 'must be at least 1 (σ²)');
    }
    if (t < 1) {
      throw ArgumentError.value(t, 'observations', 'must be at least 1');
    }
    final aicValue = -2 * logLikelihood + 2 * k;
    switch (this) {
      case InformationCriterion.aic:
        return aicValue;
      case InformationCriterion.aicc:
        if (t - k - 1 <= 0) {
          throw ArgumentError(
            'AICc is undefined for T − K − 1 ≤ 0 (T = $t, K = $k)',
          );
        }
        return aicValue + 2 * k * (k + 1) / (t - k - 1);
      case InformationCriterion.bic:
        return aicValue + (math.log(t) - 2) * k;
    }
  }
}
