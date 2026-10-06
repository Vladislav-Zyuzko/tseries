// SPDX-License-Identifier: BSD-3-Clause
//
// tseries_ctsa.c — implementation of the thin C shim over vendored `ctsa`.
// See tseries_ctsa.h for the rationale (symbol export, no struct across FFI,
// paranoid validation).

#include "tseries_ctsa.h"

// <float.h> must precede "ctsa.h": the vendored erfunc.h provides fallback
// DBL_MAX_EXP/DBL_MIN_EXP behind #ifndef, so including the system header first
// lets that guard skip them. The other order makes MSVC warn C4005 when a
// later ctsa header pulls <float.h> in and redefines both macros.
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ctsa.h"  // vendored: third_party/ctsa/src/ctsa.h (added to includes)

// Returns 1 iff every element of `a` (length `len`) is finite (no NaN/Inf).
// A NULL buffer with len == 0 is vacuously finite.
static int all_finite(const double* a, int len) {
  if (len <= 0) return 1;
  if (a == NULL) return 0;
  for (int i = 0; i < len; ++i) {
    if (!isfinite(a[i])) return 0;
  }
  return 1;
}

// ---------------------------------------------------------------------------
// Memory estimate for ctsa's exact-likelihood machinery (see
// TSERIES_MAX_STATE_BYTES in the header for the derivation and the ceiling).
//
// Mirrors the allocations ctsa makes, in 64-bit arithmetic so that the very
// configurations that overflow ctsa's own `int` sizes are measured correctly:
//   * starma (per likelihood evaluation): rbar nrbar + xnext/thetab/xrow 3*np
//   * fas154*_seas (per evaluation): phi/theta/A 3*ir + P/V 2*np
//   * forkal (forecast only, id = d + s*D): A/store 2*ird + P irz + V/xrow 2*np
// Returns LLONG_MAX when the configuration is so large that even the estimate
// would overflow -- it is then certainly over any ceiling.
static long long arma_state_bytes(long long p, long long q, long long s,
                                  long long P, long long Q, long long id,
                                  int with_forecast) {
  const long long ip = p + s * P;
  const long long iq = q + s * Q;
  const long long ir = ip > iq + 1 ? ip : iq + 1;
  // ir = 16384 already means nrbar ~ 9e15 doubles; beyond it np*np overflows.
  if (ir > 16384 || id > 16384) return 0x7fffffffffffffffLL;
  const long long np = ir * (ir + 1) / 2;
  const long long nrbar = np * (np - 1) / 2;
  long long doubles = nrbar + 5 * np + 3 * ir;
  if (with_forecast) {
    const long long ird = ir + id;
    const long long irz = ird * (ird + 1) / 2;
    doubles += irz + 2 * ird;
  }
  return doubles * (long long)sizeof(double);
}

long long tseries_model_state_bytes(int p, int d, int q, int s, int P, int D,
                                    int Q, int with_forecast) {
  const long long pp = p > 0 ? p : 0, dd = d > 0 ? d : 0, qq = q > 0 ? q : 0;
  const long long ss = s > 0 ? s : 0, PP = P > 0 ? P : 0, DD = D > 0 ? D : 0;
  const long long QQ = Q > 0 ? Q : 0;
  return arma_state_bytes(pp, qq, ss, PP, QQ, dd + ss * DD, with_forecast != 0);
}

// In-place differencing exactly as ctsa applies it: D seasonal differences at
// lag s, then d regular differences. `buf` holds `len` values on entry and the
// differenced values in its first `len - s*D - d` slots on exit. Returns the
// differenced length (may be <= 0 if the series is too short; the caller has
// already checked it is not).
static int difference_in_place(double* buf, int len, int d, int s, int D) {
  int m = len;
  for (int k = 0; k < D; ++k) {
    // Forward pass is safe in place: slot i is written only after every read
    // of it (reads are at i and i + s, both >= i).
    for (int i = 0; i + s < m; ++i) buf[i] = buf[i + s] - buf[i];
    m -= s;
  }
  for (int k = 0; k < d; ++k) {
    for (int i = 0; i + 1 < m; ++i) buf[i] = buf[i + 1] - buf[i];
    m -= 1;
  }
  return m;
}

// Returns 1 iff every element is finite AND within +/-TSERIES_MAX_ABS_VALUE.
static int all_in_range(const double* a, int len) {
  if (len <= 0) return 1;
  if (a == NULL) return 0;
  for (int i = 0; i < len; ++i) {
    if (!(fabs(a[i]) <= TSERIES_MAX_ABS_VALUE)) return 0;  // also catches NaN
  }
  return 1;
}

static double max_abs(const double* a, int len) {
  double m = 0.0;
  for (int i = 0; i < len; ++i) {
    const double v = fabs(a[i]);
    if (v > m) m = v;
  }
  return m;
}

// Root-mean-square of `a`, computed against its max so no square overflows.
static double rms_scaled(const double* a, int len) {
  const double m = max_abs(a, len);
  if (m == 0.0 || len <= 0) return 0.0;
  double acc = 0.0;
  for (int i = 0; i < len; ++i) acc += (a[i] / m) * (a[i] / m);
  return m * sqrt(acc / (double)len);
}

static double dot(const double* a, const double* b, int len) {
  double acc = 0.0;
  for (int i = 0; i < len; ++i) acc += a[i] * b[i];
  return acc;
}

// Screens the r regressor columns (column-major, length n each) the way ctsa
// will see them: differenced with (d, s, D), and alongside an intercept when
// `with_intercept`. Writes a TSERIES_EXOG_* code per column into `status` and
// returns the number of usable columns; returns -1 (with *bad_column set) if a
// differenced column overflowed to a non-finite value. Columns are judged in
// order, each against the intercept and the columns ACCEPTED before it, so the
// later member of a collinear group is the one flagged (and, under the drop
// policy, the one dropped).
//
// Needs scratch: `work` (n doubles) and `basis` ((r + 1) * m doubles), where
// m = n - d - s*D.
static int screen_regressors(const double* xreg, int n, int r, int d, int s,
                             int D, int with_intercept, int* status,
                             double* work, double* basis, int* bad_column) {
  const int m = n - d - s * D;
  int nbasis = 0;
  if (with_intercept) {
    // The intercept column is all ones (d + D == 0, so it is not differenced).
    const double inv = 1.0 / sqrt((double)m);
    for (int i = 0; i < m; ++i) basis[i] = inv;
    nbasis = 1;
  }
  int used = 0;
  for (int j = 0; j < r; ++j) {
    const double* col = xreg + (size_t)j * (size_t)n;
    const double raw_max = max_abs(col, n);
    if (raw_max == 0.0) {
      status[j] = TSERIES_EXOG_ALL_ZERO;
      continue;
    }
    memcpy(work, col, sizeof(double) * (size_t)n);
    difference_in_place(work, n, d, s, D);
    if (!all_in_range(work, m)) {
      *bad_column = j;
      return -1;
    }
    const double diff_max = max_abs(work, m);
    if (diff_max <= TSERIES_EXOG_TOL_ZERO * raw_max) {
      status[j] = TSERIES_EXOG_CONSTANT_AFTER_DIFF;
      continue;
    }
    // Normalise by the max first so no square below can overflow.
    for (int i = 0; i < m; ++i) work[i] /= diff_max;
    const double norm0 = sqrt(dot(work, work, m));

    if (with_intercept) {
      // A constant column duplicates the intercept: residual after removing
      // the mean, measured on its own, so the message can say exactly that.
      double mean = 0.0;
      for (int i = 0; i < m; ++i) mean += work[i];
      mean /= (double)m;
      double ss = 0.0;
      for (int i = 0; i < m; ++i) ss += (work[i] - mean) * (work[i] - mean);
      if (sqrt(ss) <= TSERIES_EXOG_TOL_COLLINEAR * norm0) {
        status[j] = TSERIES_EXOG_COLLINEAR_INTERCEPT;
        continue;
      }
    }
    // Two-pass modified Gram-Schmidt against everything accepted so far. The
    // second pass ("twice is enough", Kahan/Parlett) removes the component the
    // first pass leaves behind through cancellation, so the residual is a
    // reliable measure of distance from the span even for nearly dependent
    // columns.
    for (int pass = 0; pass < 2; ++pass) {
      for (int b = 0; b < nbasis; ++b) {
        const double* qb = basis + (size_t)b * (size_t)m;
        const double c = dot(qb, work, m);
        for (int i = 0; i < m; ++i) work[i] -= c * qb[i];
      }
    }
    const double res = sqrt(dot(work, work, m));
    if (!(res > TSERIES_EXOG_TOL_COLLINEAR * norm0)) {
      status[j] = TSERIES_EXOG_COLLINEAR;
      continue;
    }
    double* qn = basis + (size_t)nbasis * (size_t)m;
    for (int i = 0; i < m; ++i) qn[i] = work[i] / res;
    ++nbasis;
    status[j] = TSERIES_EXOG_USED;
    ++used;
  }
  return used;
}

// ---------------------------------------------------------------------------
// Forecasting through AS 182 (`forkal`), in place of ctsa's sarima_predict()
// and sarimax_predict(). See PROVENANCE.md, Local modifications 12.
//
// This is those two functions' preparation re-done here, statement for
// statement and in the same floating-point order, with two deliberate
// differences:
//
//  1. ctsa copies the fit's residuals with `resid[i] = obj->res[i]` for
//     i < N, but `obj->res` holds only N - d - s*D of them. For every model
//     with d + s*D > (number of estimated coefficients)^2 -- which includes
//     every differenced white-noise order, where that number is 0 -- this reads
//     past the end of the fit object's heap block (AddressSanitizer:
//     heap-buffer-overflow, ctsa.c sarima_predict/sarimax_predict). The copied
//     values are never used: karma() overwrites resid[0 .. N-id) before forkal()
//     reads any of it, and nothing reads beyond. Here the buffer is simply
//     zero-filled, so the results are bit-identical and the read is gone.
//
//  2. forkal() refuses ip == iq == 0 (ifault 4) and returns WITHOUT writing
//     the forecasts or their MSEs; ctsa's wrappers ignore the code and hand
//     back whatever the caller's buffers held. A pure white-noise model is the
//     AR(1) model with phi = 0 -- the same state-space form AS 182 already
//     uses for every ir == 1 model (P = 1/(1 - phi^2) = 1, no starma call) --
//     so it is passed to forkal() as exactly that. Point forecasts and MSEs then
//     come out of the same algorithm, with the same sigma^2 (forkal's own MLE
//     from the one-step residuals), as for every other order.
//
// The parameters mirror the ctsa object fields: `ophi`/`otheta`/`oPHI`/
// `oTHETA` as reported (MA blocks with ctsa's own sign), `mean` is the fit's
// mean (used only when d == D == 0), and `exog`/`xreg`/`newxreg` the r
// regression coefficients and regressors exactly as ctsa saw them (r == 0 on
// the SARIMA path). Returns TSERIES_OK or TSERIES_ERR_INIT_FAILED.
static int shim_forecast(const double* inp, int N, int p, int d, int q, int s,
                         int P, int D, int Q, const double* ophi,
                         const double* otheta, const double* oPHI,
                         const double* oTHETA, double mean, int r,
                         const double* exog, const double* xreg,
                         const double* newxreg, int L, double* xpred,
                         double* amse) {
  const int ip = p + s * P;
  const int iq = q + s * Q;
  int ir = ip;
  if (ir < 1 + iq) ir = 1 + iq;
  const int id = d + D * s;

  double* coef1 = (double*)malloc(sizeof(double) * (size_t)(d + 1));
  double* coef2 = (double*)malloc(sizeof(double) * (size_t)(D * s + 1));
  double* delta = (double*)malloc(sizeof(double) * (size_t)(id + 1));
  double* W = (double*)malloc(sizeof(double) * (size_t)N);
  double* resid = (double*)calloc((size_t)N, sizeof(double));
  double* phi = (double*)malloc(sizeof(double) * (size_t)ir);
  double* theta = (double*)malloc(sizeof(double) * (size_t)ir);
  int status = TSERIES_OK;
  if (coef1 == NULL || coef2 == NULL || delta == NULL || W == NULL ||
      resid == NULL || phi == NULL || theta == NULL) {
    status = TSERIES_ERR_INIT_FAILED;
    goto done;
  }

  double wmean = 0.0;
  coef1[0] = coef2[0] = 1.0;
  if (d == 0 && D == 0) {
    delta[0] = 1.0;
    wmean = mean;
  }
  if (d > 0) deld(d, coef1);
  if (D > 0) delds(D, s, coef2);
  conv(coef1, d + 1, coef2, D * s + 1, delta);
  for (int i = 1; i <= id; ++i) delta[i] = -1.0 * delta[i];

  for (int i = 0; i < N; ++i) {
    W[i] = inp[i];
    if (d == 0 && D == 0) W[i] -= wmean;
    for (int j = 0; j < r; ++j) W[i] -= exog[j] * xreg[(size_t)j * N + i];
  }
  for (int i = 0; i < ir; ++i) phi[i] = theta[i] = 0.0;
  for (int i = 0; i < p; ++i) phi[i] = ophi[i];
  for (int i = 0; i < q; ++i) theta[i] = -1.0 * otheta[i];
  for (int j = 0; j < P; ++j) {
    phi[(j + 1) * s - 1] += oPHI[j];
    for (int i = 0; i < p; ++i) phi[(j + 1) * s + i] -= ophi[i] * oPHI[j];
  }
  for (int j = 0; j < Q; ++j) {
    theta[(j + 1) * s - 1] -= oTHETA[j];
    for (int i = 0; i < q; ++i) theta[(j + 1) * s + i] += otheta[i] * oTHETA[j];
  }

  // White noise: AR(1) with phi[0] == 0 (already zeroed above, and ir == 1).
  const int ip_forkal = (ip == 0 && iq == 0) ? 1 : ip;
  const int ifault = forkal(ip_forkal, iq, id, phi, theta, delta + 1, N, W,
                            resid, L, xpred, amse);
  if (ifault != 0) {
    // Unreachable with the arguments built above (ip_forkal >= 1, id >= 0,
    // L >= 1); never hand back buffers forkal did not write.
    status = TSERIES_ERR_INIT_FAILED;
    goto done;
  }
  for (int i = 0; i < L; ++i) {
    xpred[i] += wmean;
    for (int j = 0; j < r; ++j) xpred[i] += exog[j] * newxreg[(size_t)j * L + i];
  }

done:
  free(coef1);
  free(coef2);
  free(delta);
  free(W);
  free(resid);
  free(phi);
  free(theta);
  return status;
}

int tseries_smoke(void) {
  arima_object obj = arima_init(1, 1, 1, 50);
  if (obj == NULL) return -1;
  const int ok = (obj->p == 1 && obj->d == 1 && obj->q == 1);
  arima_free(obj);
  return ok ? 42 : -2;
}

int tseries_sarima_fit(const double* series, int n, int p, int d, int q, int s,
                       int P, int D, int Q, int method, int horizon,
                       double* out_phi, double* out_theta, double* out_bigphi,
                       double* out_bigtheta, double* out_diag,
                       double* out_forecast, double* out_stderr,
                       int* out_retval) {
  // --- Argument validation (trust nothing) -------------------------------
  if (series == NULL || out_diag == NULL || out_retval == NULL) {
    return TSERIES_ERR_NULL_ARG;
  }
  // Establish the "ctsa never ran" value before anything below can return.
  *out_retval = TSERIES_RETVAL_NOT_RUN;
  if (p > 0 && out_phi == NULL) return TSERIES_ERR_NULL_ARG;
  if (q > 0 && out_theta == NULL) return TSERIES_ERR_NULL_ARG;
  if (P > 0 && out_bigphi == NULL) return TSERIES_ERR_NULL_ARG;
  if (Q > 0 && out_bigtheta == NULL) return TSERIES_ERR_NULL_ARG;
  if (horizon > 0 && (out_forecast == NULL || out_stderr == NULL)) {
    return TSERIES_ERR_NULL_ARG;
  }
  if (n < 3 || horizon < 0) return TSERIES_ERR_BAD_SIZE;
  if (p < 0 || d < 0 || q < 0 || P < 0 || D < 0 || Q < 0 || s < 0) {
    return TSERIES_ERR_BAD_PARAM;
  }
  // Enough observations to difference away d + s*D and still fit. Computed in
  // 64-bit so an absurd s*D cannot wrap around and slip past the check.
  if ((long long)n <= (long long)d + (long long)s * D + p + q + P + Q + 1) {
    return TSERIES_ERR_BAD_SIZE;
  }
  // Refuse configurations whose exact-likelihood state ctsa would try (and
  // fail, unchecked) to allocate -- e.g. a seasonal MA term at s = 288.
  if (arma_state_bytes(p, q, s, P, Q, (long long)d + (long long)s * D,
                       horizon > 0) > TSERIES_MAX_STATE_BYTES) {
    return TSERIES_ERR_TOO_LARGE;
  }
  if (!all_finite(series, n)) return TSERIES_ERR_NON_FINITE;

  // MLE and CSS both start with ctsa's CSS fit (css_seas), which conditions on
  // the first p + s*P differenced observations. With no observation after
  // them the CSS sum is empty, its value is log(0/0) = NaN at the starting
  // point, and bfgs_min() then reads its gradient buffer uninitialised
  // (MemorySanitizer, secant.c): status and estimates differ from run to run,
  // including fits reported as successful. Refused like the same condition on
  // the SARIMAX CSS path. The MA analogue (Nd <= q + s*Q with Nd > p + s*P)
  // was checked and is NOT affected: the CSS sum is non-empty, no
  // MemorySanitizer report, identical results run to run. PROVENANCE.md,
  // Local modifications 15.
  const long long nused_all = (long long)n - d - (long long)s * D;
  if ((method == TSERIES_METHOD_MLE || method == TSERIES_METHOD_CSS) &&
      nused_all <= (long long)p + (long long)s * P) {
    return TSERIES_ERR_BAD_SIZE;
  }

  // Box-Jenkins (ctsa's nlalsms) guards -- PROVENANCE.md, Local modifications
  // 12. Each one refuses an input for which ctsa reads or writes outside its
  // own heap buffers (found with AddressSanitizer):
  //  * nothing to estimate (no ARMA term, and no mean because the series is
  //    differenced): nlalsms solves a 0x0 system and linsolve() writes past
  //    its buffers -- heap corruption, then a crash;
  //  * fewer differenced observations than the expanded AR or MA lag
  //    (p + s*P or q + s*Q): avaluem() reads past the series and residuals.
  // And one argument normalisation: for a non-seasonal model ctsa's seasonal
  // start-value step (USPE_seasonal) reads cv[0] of a zero-length buffer when
  // s == 0. With P = D = Q = 0 the period is otherwise unused, so s = 1 is
  // passed instead: the read is in bounds and no result changes.
  int s_ctsa = s;
  if (method == TSERIES_METHOD_BOX_JENKINS) {
    if (p + q + P + Q == 0 && d + D > 0) return TSERIES_ERR_BAD_PARAM;
    const long long nused_bj = (long long)n - d - (long long)s * D;
    if ((long long)p + (long long)s * P > nused_bj ||
        (long long)q + (long long)s * Q > nused_bj) {
      return TSERIES_ERR_BAD_SIZE;
    }
    // Seasonal start values (USPE_seasonal) read the autocovariances at lags
    // 0, s, ..., (P + Q)*s; autocovar() computes them only up to lag Nd - 2
    // when (P + Q + 1)*s > Nd, so for Nd <= (P + Q)*s + 1 the start values,
    // and from them the whole fit, depend on uninitialised memory
    // (MemorySanitizer; results differ run to run). Measured on a sweep over
    // s in {4, 7, 12}: no report above that bound, a report at it. Refused
    // with one season of margin: Nd >= (P + Q + 1)*s + 2. PROVENANCE.md, Local
    // modifications 15.
    if (P + Q > 0 &&
        nused_bj < (long long)(P + Q + 1) * (long long)s + 2) {
      return TSERIES_ERR_BAD_SIZE;
    }
    if (s == 0 && P + D + Q == 0) s_ctsa = 1;
  }

  sarima_object obj = sarima_init(p, d, q, s_ctsa, P, D, Q, n);
  if (obj == NULL) return TSERIES_ERR_INIT_FAILED;

  sarima_setMethod(obj, method);
  sarima_exec(obj, (double*)series);

  // Surface ctsa's own status VERBATIM, before any gate below can divert the
  // flow. The caller must be able to see the real code even for a fit we go on
  // to reject, otherwise a rejection is indistinguishable from a code we never
  // observed. This shim does not act on the value -- that is the caller's call.
  const int retval = obj->retval;
  *out_retval = retval;

  // Gather coefficients and diagnostics; validate finiteness before trusting.
  int status = TSERIES_OK;
  if (!all_finite(obj->phi, p) || !all_finite(obj->theta, q) ||
      !all_finite(obj->PHI, P) || !all_finite(obj->THETA, Q)) {
    status = TSERIES_ERR_NON_FINITE;
    goto cleanup;
  }
  {
    // Gate ONLY on fields we actually trust and actually hand back as the
    // product of the fit. `mean`/`var` qualify: sarima_init seeds them and
    // every estimation method assigns them.
    //
    // `loglik`/`aic` deliberately do NOT take part in this gate. They are
    // diagnostics, not the product, and they are not assigned on every path
    // (see below). Gating on them would let an unrelated, untrusted field veto
    // an otherwise good fit -- e.g. rejecting every Box-Jenkins fit, for which
    // ctsa never assigns either field.
    const double trusted[2] = {obj->mean, obj->var};
    if (!all_finite(trusted, 2)) {
      status = TSERIES_ERR_NON_FINITE;
      goto cleanup;
    }
    // A residual variance of zero (or below) is a degenerate model, not a
    // confident one: it comes from a series ctsa can fit exactly, e.g. a
    // constant/flatlined input. Every forecast standard error is then 0, so the
    // confidence band collapses to a point and the model claims perfect
    // certainty. Reject it here, on a field we trust and actually return,
    // rather than relying on loglik/aic blowing up to +/-Inf as a side effect
    // (which is how this used to be caught, by accident, through fields that
    // ctsa does not even assign for every method).
    if (!(obj->var > 0.0)) {
      status = TSERIES_ERR_DEGENERATE;
      goto cleanup;
    }

    // Report loglik/aic ONLY where ctsa demonstrably computed them for the
    // estimator that was requested; NaN otherwise, so "not computed" cannot be
    // mistaken for a real value.
    //
    // This is also why the reads below are guarded rather than the ctsa struct
    // being pre-initialised: sarima_init leaves loglik/aic uninitialised, so
    // for a method that never assigns them they hold malloc garbage that can
    // easily look like a plausible finite number. We are the only reader of
    // those fields, so simply NOT READING them unless ctsa assigned them is
    // both sufficient and strictly cheaper than patching vendored code.
    //   * MLE : sarima_exec assigns both.
    //   * CSS : assigns loglik only -- aic is never written on that branch.
    //   * Box-Jenkins: assigns neither -- both are read as garbage. Never read.
    //
    // The subtle case is retval 10/12: ctsa returns from as154_seas BEFORE the
    // MLE step, leaving a CSS log-likelihood in place, yet sarima_exec computes
    // `aic` from it unconditionally. That AIC is not an MLE AIC and is not
    // comparable across fits -- so it must not be reported as one.
    double loglik = NAN;
    double aic = NAN;
    if (retval == TSERIES_RETVAL_SUCCESS) {
      if (method == TSERIES_METHOD_MLE) {
        loglik = obj->loglik;
        aic = obj->aic;
      } else if (method == TSERIES_METHOD_CSS) {
        // CSS assigns a (CSS) log-likelihood -- which is what the caller asked
        // for -- but ctsa never assigns an aic on this branch.
        loglik = obj->loglik;
      }
      // Box-Jenkins: ctsa assigns neither; both stay NaN.
    }
    // Collapse any non-finite diagnostic to NaN, so that "NaN" is the single,
    // uniform way this shim says "not trustworthy" and the Dart layer has
    // exactly one value to map to null. An infinite AIC is not a usable number
    // and must not be handed back as if it were one.
    if (!isfinite(loglik)) loglik = NAN;
    if (!isfinite(aic)) aic = NAN;

    for (int i = 0; i < p; ++i) out_phi[i] = obj->phi[i];
    for (int i = 0; i < q; ++i) out_theta[i] = obj->theta[i];
    for (int i = 0; i < P; ++i) out_bigphi[i] = obj->PHI[i];
    for (int i = 0; i < Q; ++i) out_bigtheta[i] = obj->THETA[i];
    out_diag[0] = obj->mean;
    out_diag[1] = obj->var;
    out_diag[2] = loglik;
    out_diag[3] = aic;
  }

  if (horizon > 0) {
    double* amse = (double*)malloc(sizeof(double) * (size_t)horizon);
    if (amse == NULL) {
      status = TSERIES_ERR_INIT_FAILED;
      goto cleanup;
    }
    // Not sarima_predict(): see shim_forecast() for the two ctsa defects
    // (heap over-read; white-noise forecasts never written) it avoids.
    const int fstatus = shim_forecast(
        series, n, p, d, q, obj->s, P, D, Q, obj->phi, obj->theta, obj->PHI,
        obj->THETA, obj->mean, 0, NULL, NULL, NULL, horizon, out_forecast, amse);
    if (fstatus != TSERIES_OK) {
      free(amse);
      status = fstatus;
      goto cleanup;
    }
    if (!all_finite(out_forecast, horizon) || !all_finite(amse, horizon)) {
      free(amse);
      status = TSERIES_ERR_NON_FINITE;
      goto cleanup;
    }
    for (int i = 0; i < horizon; ++i) {
      // ctsa returns per-step MSE; expose the standard error.
      out_stderr[i] = sqrt(amse[i]);
    }
    free(amse);
  }

cleanup:
  sarima_free(obj);
  return status;
}

// Sign of parameter `a` in ctsa's internal vector relative to the reported
// coefficient: the MA blocks are stored negated (ctsa optimises -theta).
static int ctsa_param_sign(int a, int p, int q, int P, int Q) {
  if (a >= p && a < p + q) return -1;
  if (a >= p + q + P && a < p + q + P + Q) return -1;
  return 1;
}

// ---------------------------------------------------------------------------
// SARIMAX (regression with seasonal ARIMA errors). See the header for the
// contract; the comments here explain the HOW.
// ---------------------------------------------------------------------------
int tseries_sarimax_fit(const double* series, int n, const double* xreg, int r,
                        const double* future_xreg, int horizon, int p, int d,
                        int q, int s, int P, int D, int Q, int method,
                        int include_mean, int exog_policy, double* out_phi,
                        double* out_theta, double* out_bigphi,
                        double* out_bigtheta, double* out_beta,
                        int* out_column_status, double* out_diag,
                        double* out_vcov, double* out_forecast,
                        double* out_stderr, int* out_retval,
                        int* out_bad_column) {
  // --- Argument validation (trust nothing) -------------------------------
  if (series == NULL || out_diag == NULL || out_retval == NULL ||
      out_bad_column == NULL) {
    return TSERIES_ERR_NULL_ARG;
  }
  *out_retval = TSERIES_RETVAL_NOT_RUN;
  *out_bad_column = -1;
  if (n < 3 || horizon < 0 || r < 0) return TSERIES_ERR_BAD_SIZE;
  if (p > 0 && out_phi == NULL) return TSERIES_ERR_NULL_ARG;
  if (q > 0 && out_theta == NULL) return TSERIES_ERR_NULL_ARG;
  if (P > 0 && out_bigphi == NULL) return TSERIES_ERR_NULL_ARG;
  if (Q > 0 && out_bigtheta == NULL) return TSERIES_ERR_NULL_ARG;
  if (r > 0 &&
      (xreg == NULL || out_beta == NULL || out_column_status == NULL)) {
    return TSERIES_ERR_NULL_ARG;
  }
  if (r > 0 && horizon > 0 && future_xreg == NULL) return TSERIES_ERR_NULL_ARG;
  if (horizon > 0 && (out_forecast == NULL || out_stderr == NULL)) {
    return TSERIES_ERR_NULL_ARG;
  }
  if (p < 0 || d < 0 || q < 0 || P < 0 || D < 0 || Q < 0 || s < 0) {
    return TSERIES_ERR_BAD_PARAM;
  }
  // A seasonal term needs a real period. With no seasonal term ctsa forces
  // s = 0 itself (sarimax_init); normalise here too so every size computed
  // below agrees with what ctsa will actually do.
  if (P + D + Q > 0 && s < 1) return TSERIES_ERR_BAD_PARAM;
  if (P + D + Q == 0) s = 0;
  if (include_mean != 0 && include_mean != 1) return TSERIES_ERR_BAD_PARAM;
  if (exog_policy != TSERIES_EXOG_POLICY_STRICT &&
      exog_policy != TSERIES_EXOG_POLICY_DROP) {
    return TSERIES_ERR_BAD_PARAM;
  }
  // Explicit mapping onto ctsa's sarimax numbering. An unknown method must be
  // refused here: sarimax_exec() calls exit(-1) on one.
  int ctsa_method;
  switch (method) {
    case TSERIES_SARIMAX_METHOD_CSS_MLE:
      ctsa_method = 0;
      break;
    case TSERIES_SARIMAX_METHOD_MLE:
      ctsa_method = 1;
      break;
    case TSERIES_SARIMAX_METHOD_CSS:
      ctsa_method = 2;
      break;
    default:
      return TSERIES_ERR_BAD_PARAM;
  }

  const int nd = d + D;
  // ctsa estimates an intercept only when nothing is differenced.
  const int mean_term = (include_mean == 1 && nd == 0) ? 1 : 0;
  const long long nused_ll = (long long)n - d - (long long)s * D;
  const int arma = p + q + P + Q;
  // Coarse size floor before any allocation; the exact degrees-of-freedom
  // check (after screening may have dropped columns) follows below.
  if (nused_ll - arma - mean_term < TSERIES_MIN_RESIDUAL_DOF) {
    return TSERIES_ERR_BAD_SIZE;
  }
  if (arma_state_bytes(p, q, s, P, Q, (long long)d + (long long)s * D,
                       horizon > 0) > TSERIES_MAX_STATE_BYTES) {
    return TSERIES_ERR_TOO_LARGE;
  }
  const int nused = (int)nused_ll;
  // ctsa's CSS step for SARIMAX (cssx/fcssx, run by CSS_MLE and CSS)
  // conditions on the first p + s*P differenced observations and zeroes that
  // many slots starting N doubles into a buffer sized for (3 + regressors) * N
  // values. With no observation left after the conditioning lags it has
  // nothing to fit, and once p + s*P exceeds (2 + regressors) * N it writes past
  // the end of that heap block (AddressSanitizer: heap-buffer-overflow WRITE,
  // emle.c fcssx). Refuse it -- PROVENANCE.md, Local modifications 12.
  if ((ctsa_method == 0 || ctsa_method == 2) &&
      nused_ll <= (long long)p + (long long)s * P) {
    return TSERIES_ERR_BAD_SIZE;
  }

  // --- Every input value finite and within +/-TSERIES_MAX_ABS_VALUE ---------
  if (!all_in_range(series, n)) return TSERIES_ERR_NON_FINITE;
  for (int j = 0; j < r; ++j) {
    if (!all_in_range(xreg + (size_t)j * (size_t)n, n) ||
        (horizon > 0 &&
         !all_in_range(future_xreg + (size_t)j * (size_t)horizon, horizon))) {
      *out_bad_column = j;
      return TSERIES_ERR_NON_FINITE;
    }
  }

  // Every buffer below is owned by this function and released at `cleanup`
  // on every path (free(NULL) is a no-op).
  int status = TSERIES_OK;
  double* work = NULL;     // n: differencing scratch
  double* basis = NULL;    // (r+1)*nused: Gram-Schmidt basis
  double* x_used = NULL;   // n*r_used: compacted regressors (column-major)
  double* fx_used = NULL;  // horizon*r_used: compacted future regressors
  int* used_index = NULL;  // r_used: original column of each used regressor
  int* scale_exp = NULL;   // r_used: column k is passed to ctsa as x * 2^-e_k
  double* amse = NULL;     // horizon
  sarimax_object obj = NULL;
  int r_used = 0;
  int retval = TSERIES_RETVAL_NOT_RUN;

  work = (double*)malloc(sizeof(double) * (size_t)n);
  if (work == NULL) return TSERIES_ERR_INIT_FAILED;

  // Differencing can overflow perfectly finite input (1e308 minus -1e308);
  // ctsa would carry the resulting Inf straight into its optimiser.
  memcpy(work, series, sizeof(double) * (size_t)n);
  difference_in_place(work, n, d, s, D);
  if (!all_in_range(work, nused)) {
    status = TSERIES_ERR_NON_FINITE;
    goto cleanup;
  }
  // Scale of the differenced series: the target each regressor is brought to.
  const double rms_y = rms_scaled(work, nused);

  // --- Regressor screening ---------------------------------------------------
  if (r > 0) {
    basis = (double*)malloc(sizeof(double) * (size_t)(r + 1) * (size_t)nused);
    used_index = (int*)malloc(sizeof(int) * (size_t)r);
    scale_exp = (int*)malloc(sizeof(int) * (size_t)r);
    if (basis == NULL || used_index == NULL || scale_exp == NULL) {
      status = TSERIES_ERR_INIT_FAILED;
      goto cleanup;
    }
    int bad = -1;
    r_used = screen_regressors(xreg, n, r, d, s, D, mean_term,
                               out_column_status, work, basis, &bad);
    if (r_used < 0) {
      *out_bad_column = bad;
      status = TSERIES_ERR_NON_FINITE;
      goto cleanup;
    }
    if (r_used < r && exog_policy == TSERIES_EXOG_POLICY_STRICT) {
      for (int j = 0; j < r; ++j) {
        if (out_column_status[j] != TSERIES_EXOG_USED) {
          *out_bad_column = j;
          break;
        }
      }
      status = TSERIES_ERR_EXOG_DEFECT;
      goto cleanup;
    }
  }

  // --- Exact degrees-of-freedom guard ---------------------------------------
  // ctsa's regression step evaluates tinv(df = nused - regressors) and calls
  // exit(1) when df <= 0. Demand a margin over every estimated coefficient.
  if (nused - arma - mean_term - r_used < TSERIES_MIN_RESIDUAL_DOF) {
    status = TSERIES_ERR_BAD_SIZE;
    goto cleanup;
  }
  // ctsa's regression set-up (cssx() and as154x(), emle.c) seasonally
  // differences each regressor into XX + (n - s*D)*i, i < r_used, but XX holds
  // only nused*(r_used + 1) doubles. With D > 0 that runs past the heap block
  // once (n - s*D)*r_used > nused*(r_used + 1), i.e. d*r_used > nused
  // (AddressSanitizer: heap-buffer-overflow WRITE in diffs(), talg.c, from
  // cssx; e.g. 3 regressors under (0,2,0)(0,1,0)[4] on 11 points). Refuse it
  // -- PROVENANCE.md, Local modifications 16.
  if (D > 0 && r_used > 0 && (long long)d * r_used > nused) {
    status = TSERIES_ERR_BAD_SIZE;
    goto cleanup;
  }

  // --- Compact and scale the usable regressors (column-major) ---------------
  // ctsa optimises the regression coefficients in raw units alongside the ARMA
  // coefficients (BFGS from a unit Hessian, unit max step, unscaled
  // finite-difference Hessian). It is NOT scale-invariant: with a regressor in
  // units 1e6 times too large it was measured returning retval 1 at an optimum
  // whose log-likelihood was 52 worse, and at 1e-6 it stopped on the iteration
  // cap. Each column is therefore handed to ctsa multiplied by an exact power
  // of two, 2^-e, chosen so its differenced RMS matches the differenced
  // series' -- which puts the coefficient ctsa sees at O(1) whenever the
  // regressor explains anything. A power of two changes no bit of the data's
  // mantissas, the model is identical (beta_ctsa = beta * 2^e), and the
  // coefficients and their covariance are mapped back exactly below.
  if (r_used > 0) {
    x_used = (double*)malloc(sizeof(double) * (size_t)n * (size_t)r_used);
    if (horizon > 0) {
      fx_used =
          (double*)malloc(sizeof(double) * (size_t)horizon * (size_t)r_used);
    }
    if (x_used == NULL || (horizon > 0 && fx_used == NULL)) {
      status = TSERIES_ERR_INIT_FAILED;
      goto cleanup;
    }
    int k = 0;
    for (int j = 0; j < r; ++j) {
      if (out_column_status[j] != TSERIES_EXOG_USED) continue;
      used_index[k] = j;
      const double* col = xreg + (size_t)j * (size_t)n;
      memcpy(work, col, sizeof(double) * (size_t)n);
      difference_in_place(work, n, d, s, D);
      const double rms_x = rms_scaled(work, nused);
      int e = 0;
      if (rms_x > 0.0 && rms_y > 0.0) {
        e = (int)lround(log2(rms_x / rms_y));
      }
      scale_exp[k] = e;
      double* dst = x_used + (size_t)k * (size_t)n;
      for (int i = 0; i < n; ++i) dst[i] = ldexp(col[i], -e);
      if (horizon > 0) {
        const double* fcol = future_xreg + (size_t)j * (size_t)horizon;
        double* fdst = fx_used + (size_t)k * (size_t)horizon;
        for (int i = 0; i < horizon; ++i) fdst[i] = ldexp(fcol[i], -e);
      }
      ++k;
    }
  }

  // --- Fit ------------------------------------------------------------------
  // NOTE the argument order: sarimax_init takes s AFTER Q, unlike sarima_init.
  obj = sarimax_init(p, d, q, P, D, Q, s, r_used, mean_term, n);
  if (obj == NULL) {
    status = TSERIES_ERR_INIT_FAILED;
    goto cleanup;
  }
  sarimax_setMethod(obj, ctsa_method);
  sarimax_exec(obj, (double*)series, x_used);

  // Surface ctsa's status verbatim before any gate, as on the SARIMA path.
  retval = obj->retval;
  *out_retval = retval;

  // Code 7 means ctsa's regression step refused the design and NOTHING was
  // estimated (every coefficient is still the zero sarimax_init wrote). After
  // our own screening and column scaling this is only reachable through
  // ctsa's scale-dependent rank tolerance against the (unscaled) intercept
  // column, so it is a distinct failure, never a fit.
  if (retval == TSERIES_RETVAL_COLLINEAR_EXOG) {
    status = TSERIES_ERR_NATIVE_COLLINEAR;
    goto cleanup;
  }

  if (!all_finite(obj->phi, p) || !all_finite(obj->theta, q) ||
      !all_finite(obj->PHI, P) || !all_finite(obj->THETA, Q) ||
      !all_finite(obj->exog, r_used)) {
    status = TSERIES_ERR_NON_FINITE;
    goto cleanup;
  }
  {
    const double trusted[2] = {obj->mean, obj->var};
    if (!all_finite(trusted, 2)) {
      status = TSERIES_ERR_NON_FINITE;
      goto cleanup;
    }
    if (!(obj->var > 0.0)) {
      status = TSERIES_ERR_DEGENERATE;
      goto cleanup;
    }

    const int exact = (method == TSERIES_SARIMAX_METHOD_CSS_MLE ||
                       method == TSERIES_SARIMAX_METHOD_MLE);
    double loglik = NAN;
    double aic = NAN;
    if (retval == TSERIES_RETVAL_SUCCESS) {
      // sarimax_exec assigns loglik for every method and aic for the two exact
      // ones. On any other retval (4/10/12/15) CSS_MLE returns before its MLE
      // step and both would describe a CSS fit, so they are withheld.
      loglik = obj->loglik;
      if (exact) aic = obj->aic;
    }
    if (!isfinite(loglik)) loglik = NAN;
    if (!isfinite(aic)) aic = NAN;

    for (int i = 0; i < p; ++i) out_phi[i] = obj->phi[i];
    for (int i = 0; i < q; ++i) out_theta[i] = obj->theta[i];
    for (int i = 0; i < P; ++i) out_bigphi[i] = obj->PHI[i];
    for (int i = 0; i < Q; ++i) out_bigtheta[i] = obj->THETA[i];
    for (int j = 0; j < r; ++j) out_beta[j] = NAN;
    // Undo the column scaling: beta = beta_ctsa * 2^-e.
    for (int k = 0; k < r_used; ++k) {
      out_beta[used_index[k]] = ldexp(obj->exog[k], -scale_exp[k]);
    }
    out_diag[0] = obj->mean;
    out_diag[1] = obj->var;
    out_diag[2] = loglik;
    out_diag[3] = aic;

    if (out_vcov != NULL) {
      // Caller layout: [phi | theta | PHI | THETA | mu? | beta_0..beta_{r-1}].
      // ctsa layout:   [phi | -theta | PHI | -THETA | mu? | used betas].
      const int fixed = arma + mean_term;
      const int kc = fixed + r;
      const int ki = fixed + r_used;
      for (int i = 0; i < kc * kc; ++i) out_vcov[i] = NAN;
      int ok = exact && retval == TSERIES_RETVAL_SUCCESS &&
               obj->lvcov == ki * ki && all_finite(obj->vcov, ki * ki);
      for (int i = 0; ok && i < ki; ++i) {
        if (!(obj->vcov[i * ki + i] > 0.0)) ok = 0;
      }
      if (ok) {
        for (int a = 0; a < ki; ++a) {
          const int ca = a < fixed ? a : fixed + used_index[a - fixed];
          const int sa = ctsa_param_sign(a, p, q, P, Q);
          for (int b = 0; b < ki; ++b) {
            const int cb = b < fixed ? b : fixed + used_index[b - fixed];
            const int sb = ctsa_param_sign(b, p, q, P, Q);
            // Undo the column scaling on every beta index (exact, 2^-e).
            const int ea = a < fixed ? 0 : scale_exp[a - fixed];
            const int eb = b < fixed ? 0 : scale_exp[b - fixed];
            out_vcov[ca * kc + cb] =
                ldexp((double)(sa * sb) * obj->vcov[a * ki + b], -(ea + eb));
          }
        }
      }
    }
  }

  if (horizon > 0) {
    amse = (double*)malloc(sizeof(double) * (size_t)horizon);
    if (amse == NULL) {
      status = TSERIES_ERR_INIT_FAILED;
      goto cleanup;
    }
    // Not sarimax_predict(): see shim_forecast().
    status = shim_forecast(series, n, p, d, q, obj->s, P, D, Q, obj->phi,
                           obj->theta, obj->PHI, obj->THETA, obj->mean, obj->r,
                           obj->exog, x_used, fx_used, horizon, out_forecast,
                           amse);
    if (status != TSERIES_OK) goto cleanup;
    if (!all_finite(out_forecast, horizon) || !all_finite(amse, horizon)) {
      status = TSERIES_ERR_NON_FINITE;
      goto cleanup;
    }
    for (int i = 0; i < horizon; ++i) out_stderr[i] = sqrt(amse[i]);
  }

cleanup:
  if (obj != NULL) sarimax_free(obj);
  free(amse);
  free(fx_used);
  free(x_used);
  free(used_index);
  free(scale_exp);
  free(basis);
  free(work);
  return status;
}
