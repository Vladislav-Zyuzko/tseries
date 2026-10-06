// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// Clean-room replacements for ctsa's LGPL-3.0 `erfinv` / `erfcinv` (upstream
// `src/erfunc.c`, removed from this tree). Written from published sources
// only; the upstream body was not consulted. Process and sources:
// third_party/ctsa/PROVENANCE.md, "Clean-room replacements".
//
// `erf` / `erfc` themselves are no longer defined by ctsa: the vendored code
// now calls the C99 <math.h> functions from the platform's libm (the upstream
// definitions shadowed them under the same names).
//
// Method:
//  1. Starting value from P. J. Acklam's rational approximation of the
//     standard normal quantile (relative error < 1.15e-9), "An algorithm for
//     computing the inverse normal cumulative distribution function" (2003),
//     using erfinv(x) = -Phi^-1((1 - x) / 2) / sqrt(2) and
//     erfcinv(c) = -Phi^-1(c / 2) / sqrt(2).
//  2. Two Halley iterations against libm's erf / erfc, which take the
//     starting value to full double precision. With f(y) = erf(y) - x we have
//     f'(y) = 2/sqrt(pi) * exp(-y^2) and f''(y) = -2y f'(y), so Halley's step
//     y <- y - 2 f f' / (2 f'^2 - f f'') reduces to y <- y - r / (1 + y r),
//     r = f / f'. The same holds for g(y) = erfc(y) - c with g' = -f'.
//  Which residual is iterated is chosen so that it is computed without
//  cancellation: erf(y) - x where |x| <= 1/2 (erf is relatively accurate near
//  0), erfc(y) - c in the tail, where c = 1 - |x| is exact for |x| >= 1/2
//  (Sterbenz) and erfc keeps full relative accuracy down to underflow.

#include <float.h>
#include <math.h>

#include "erfunc.h"

#define TSERIES_TWO_OVER_SQRT_PI 1.1283791670955125739
#define TSERIES_SQRT1_2 0.70710678118654752440

// Acklam's coefficients (central region and tail).
static const double kA[6] = {-3.969683028665376e+01, 2.209460984245205e+02,
                             -2.759285104469687e+02, 1.383577518672690e+02,
                             -3.066479806614716e+01, 2.506628277459239e+00};
static const double kB[5] = {-5.447609879822406e+01, 1.615858368580409e+02,
                             -1.556989798598866e+02, 6.680131188771972e+01,
                             -1.328068155288572e+01};
static const double kC[6] = {-7.784894002430293e-03, -3.223964580411365e-01,
                             -2.400758277161838e+00, -2.549732539343734e+00,
                             4.374664141464968e+00, 2.938163982698783e+00};
static const double kD[4] = {7.784695709041462e-03, 3.224671290700398e-01,
                             2.445134137142996e+00, 3.754408661907416e+00};

// Approximate standard normal quantile Phi^-1(p) for 0 < p <= 1/2 (result
// <= 0). Starting value only — refined below.
static double normal_quantile_lower(double p) {
  if (p < 0.02425) {
    const double q = sqrt(-2.0 * log(p));
    return (((((kC[0] * q + kC[1]) * q + kC[2]) * q + kC[3]) * q + kC[4]) * q +
            kC[5]) /
           ((((kD[0] * q + kD[1]) * q + kD[2]) * q + kD[3]) * q + 1.0);
  } else {
    const double q = p - 0.5;
    const double r = q * q;
    return (((((kA[0] * r + kA[1]) * r + kA[2]) * r + kA[3]) * r + kA[4]) * r +
            kA[5]) *
           q /
           (((((kB[0] * r + kB[1]) * r + kB[2]) * r + kB[3]) * r + kB[4]) * r +
            1.0);
  }
}

// One Halley step for erf(y) = target (use_erfc == 0) or
// erfc(y) = target (use_erfc != 0). Returns the updated y.
static double halley_step(double y, double target, int use_erfc) {
  const double slope = TSERIES_TWO_OVER_SQRT_PI * exp(-y * y);  // erf'(y)
  double r;
  if (!(slope > 0.0) || !isfinite(slope)) return y;  // deep tail: keep y
  if (use_erfc) {
    r = -(erfc(y) - target) / slope;
  } else {
    r = (erf(y) - target) / slope;
  }
  return y - r / (1.0 + y * r);
}

// erfcinv for 0 < c <= 1 (result >= 0), iterating on erfc.
static double erfcinv_upper(double c) {
  double y = -normal_quantile_lower(0.5 * c) * TSERIES_SQRT1_2;
  y = halley_step(y, c, 1);
  y = halley_step(y, c, 1);
  return y;
}

double erfinv(double x) {
  double ax, y;
  if (isnan(x)) return x;
  if (x == 0.0) return x;  // keeps the sign of zero
  ax = fabs(x);
  if (ax > 1.0) return NAN;
  if (ax == 1.0) return x > 0.0 ? INFINITY : -INFINITY;
  if (ax <= 0.5) {
    y = -normal_quantile_lower(0.5 * (1.0 - ax)) * TSERIES_SQRT1_2;
    y = halley_step(y, ax, 0);
    y = halley_step(y, ax, 0);
  } else {
    y = erfcinv_upper(1.0 - ax);  // 1 - ax is exact here
  }
  return x > 0.0 ? y : -y;
}

double erfcinv(double c) {
  if (isnan(c)) return c;
  if (c < 0.0 || c > 2.0) return NAN;
  if (c == 0.0) return INFINITY;
  if (c == 2.0) return -INFINITY;
  if (c > 1.0) return -erfcinv(2.0 - c);  // 2 - c is exact for 1 <= c <= 2
  if (c >= 0.5) return erfinv(1.0 - c);  // 1 - c is exact, in [0, 1/2]
  return erfcinv_upper(c);
}
