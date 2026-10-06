// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// Clean-room replacements for the two LGPL-3.0 sections of ctsa's
// `src/dist.c` (`r8_max` and `betainv`), which were cut from the vendored
// file. Written from published sources only; the upstream bodies were not
// consulted. Process and sources: third_party/ctsa/PROVENANCE.md,
// "Clean-room replacements".
//
// betainv(alpha, a, b) returns x in [0, 1] with I_x(a, b) = alpha, where I is
// the regularized incomplete beta function as computed by ctsa's own `ibeta`
// (dist.c; public-domain/CC0 section, continued fraction of Abramowitz &
// Stegun 26.5.8). Inverting the very function the rest of ctsa uses keeps
// `tinv`, `finv` and `tcdf`/`fcdf` mutually consistent.
//
// Method:
//  * Starting value: Abramowitz & Stegun, Handbook of Mathematical Functions
//    (1964), 26.5.22 (normal-quantile based) when a > 1 and b > 1; otherwise
//    the leading terms of the two tail expansions of I_x(a, b),
//      I_x ~ x^a / (a B(a, b))            as x -> 0,
//      1 - I_x ~ (1 - x)^b / (b B(a, b))  as x -> 1,
//    solved for x, keeping whichever reproduces alpha better.
//  * Refinement: Halley's method on f(x) = I_x(a, b) - alpha, with
//    f' = x^(a-1) (1-x)^(b-1) / B(a, b) and f''/f' = (a-1)/x - (b-1)/(1-x),
//    safeguarded by a bracket [lo, hi] that always contains the root: any step
//    that leaves the bracket is replaced by bisection (geometric near 0 and 1,
//    so extreme tails converge in a bounded number of steps), then a bounded
//    ulp-level polish against ctsa's computed ibeta (see polish()).

#include <float.h>
#include <math.h>

#include "dist.h"

#define TSERIES_BETAINV_MAX_ITER 200

double r8_max(double x, double y) { return x < y ? y : x; }

// Upper-tail standard normal quantile: y with Q(y) = 1 - Phi(y) = p.
static double normal_upper_quantile(double p) {
  return 1.4142135623730950488 * erfcinv(2.0 * p);
}

// Initial approximation; may lie outside (0, 1) — the caller checks.
static double betainv_start(double alpha, double a, double b, double lbeta) {
  double x_lo, x_hi, r_lo, r_hi;
  if (a > 1.0 && b > 1.0) {
    // A&S 26.5.22.
    const double y = normal_upper_quantile(alpha);
    const double lambda = (y * y - 3.0) / 6.0;
    const double ia = 1.0 / (2.0 * a - 1.0);
    const double ib = 1.0 / (2.0 * b - 1.0);
    const double h = 2.0 / (ia + ib);
    const double w = y * sqrt(h + lambda) / h -
                     (ib - ia) * (lambda + 5.0 / 6.0 - 2.0 / (3.0 * h));
    const double x = a / (a + b * exp(2.0 * w));
    if (x > 0.0 && x < 1.0) return x;
  }
  // Tail expansions, each solved for x.
  x_lo = exp((log(alpha) + log(a) + lbeta) / a);
  x_hi = -expm1((log1p(-alpha) + log(b) + lbeta) / b);
  r_lo = (x_lo > 0.0 && x_lo < 1.0) ? fabs(ibeta(x_lo, a, b) - alpha) : HUGE_VAL;
  r_hi = (x_hi > 0.0 && x_hi < 1.0) ? fabs(ibeta(x_hi, a, b) - alpha) : HUGE_VAL;
  if (r_lo <= r_hi && r_lo < HUGE_VAL) return x_lo;
  if (r_hi < HUGE_VAL) return x_hi;
  return 0.5;
}

// Bisection point of (lo, hi): geometric when the bracket spans orders of
// magnitude towards 0 (or towards 1, in 1 - x), arithmetic otherwise.
static double bisect(double lo, double hi) {
  if (lo == 0.0 && hi < 0.25) return hi * 0.0625;
  if (lo > 0.0 && hi < 0.5 && hi > 4.0 * lo) return sqrt(lo * hi);
  if (hi == 1.0 && lo > 0.75) return 1.0 - (1.0 - lo) * 0.0625;
  if (hi < 1.0 && lo > 0.5 && (1.0 - lo) > 4.0 * (1.0 - hi)) {
    return 1.0 - sqrt((1.0 - lo) * (1.0 - hi));
  }
  return 0.5 * (lo + hi);
}

// Final polish on the *computed* residual r(y) = ibeta(y) - alpha: ctsa's
// ibeta is a staircase at ulp level, so the iteration above can stop on a
// plateau next to a better double. Walk one ulp at a time towards the root
// (I is increasing in x) while |r| does not get worse, and return the best
// point seen. Bounded: at most TSERIES_BETAINV_POLISH_STEPS evaluations.
#define TSERIES_BETAINV_POLISH_STEPS 32

static double polish(double x, double r, double alpha, double a, double b) {
  double best_x = x, best_r = fabs(r);
  int k;
  for (k = 0; k < TSERIES_BETAINV_POLISH_STEPS && r != 0.0; ++k) {
    const double y = nextafter(x, r > 0.0 ? 0.0 : 1.0);
    double ry;
    if (!(y > 0.0 && y < 1.0)) break;
    ry = ibeta(y, a, b) - alpha;
    if (fabs(ry) > fabs(r) && (ry > 0.0) == (r > 0.0)) break;  // worse, same side
    if (fabs(ry) < best_r) {
      best_r = fabs(ry);
      best_x = y;
    }
    if ((ry > 0.0) != (r > 0.0) && ry != 0.0) break;  // crossed the root
    x = y;
    r = ry;
  }
  return best_x;
}

double betainv(double alpha, double a, double b) {
  double lbeta, x, lo = 0.0, hi = 1.0;
  int iter;

  if (isnan(alpha) || isnan(a) || isnan(b)) return NAN;
  if (!(a > 0.0) || !(b > 0.0) || !isfinite(a) || !isfinite(b)) return NAN;
  if (alpha <= 0.0) return 0.0;
  if (alpha >= 1.0) return 1.0;

  lbeta = beta_log(a, b);
  x = betainv_start(alpha, a, b, lbeta);
  if (!(x > 0.0 && x < 1.0)) x = 0.5;

  for (iter = 0; iter < TSERIES_BETAINV_MAX_ITER; ++iter) {
    const double f = ibeta(x, a, b) - alpha;
    double x_new = NAN;
    double log_pdf, pdf;

    if (f == 0.0) return x;
    if (f < 0.0) {
      lo = x;
    } else {
      hi = x;
    }

    log_pdf = (a - 1.0) * log(x) + (b - 1.0) * log1p(-x) - lbeta;
    pdf = exp(log_pdf);
    if (pdf > 0.0 && isfinite(pdf)) {
      double step = f / pdf;  // Newton
      const double curvature = (a - 1.0) / x - (b - 1.0) / (1.0 - x);
      const double denom = 1.0 - 0.5 * step * curvature;
      if (denom > 0.5 && denom < 2.0) step /= denom;  // Halley, when tame
      x_new = x - step;
    }
    if (!(x_new > lo && x_new < hi)) x_new = bisect(lo, hi);

    if (fabs(x_new - x) <= 2.0 * DBL_EPSILON * x_new ||
        hi - lo <= 2.0 * DBL_EPSILON * lo) {
      x = x_new;
      break;
    }
    x = x_new;
  }
  return polish(x, ibeta(x, a, b) - alpha, alpha, a, b);
}
