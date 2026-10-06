// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// Replacements for ctsa's `archeck()` and `invertroot()` (upstream
// `src/talg.c`), whose upstream bodies were C ports of R's `arCheck()` and
// `maInvert()` (R `stats` package, GPL-2+) and have been removed from the
// vendored tree. See third_party/ctsa/PROVENANCE.md, "Clean-room
// replacements" (round 3). This file was written from a textbook
// specification (Box, Jenkins, Reinsel & Ljung, "Time Series Analysis", 5th
// ed., Wiley 2015, sec. 3.2-3.4; Brockwell & Davis, "Time Series: Theory and
// Methods", 2nd ed., Springer 1991, sec. 3.1 and 4.4) without reading the
// upstream bodies or R's sources. Sign conventions, return values and edge
// cases were established by black-box tests of the upstream binary.
//
// Polynomials (coefficients at lags 1..n, as the ctsa callers store them):
//   AR  phi(z)   = 1 - ar[0] z - ar[1] z^2 - ... - ar[p-1] z^p
//   MA  theta(z) = 1 + ma[0] z + ma[1] z^2 + ... + ma[q-1] z^q
// Trailing zero coefficients lower the true degree m; only the m roots of the
// true-degree polynomial are considered (a zero-degree polynomial has none).
// Roots come from `polyroot()` (ctsa_replacements/polyroot.c).
//
// int archeck(int p, double *ar)
//   Returns 1 if the AR part is stationary -- no root of phi lies strictly
//   inside the unit circle -- and 0 otherwise. As upstream, a root of modulus
//   exactly 1 does not make the result 0 (there is no tolerance band: roots
//   within rounding of the circle go either way, as they did upstream).
//   p <= 0 or true degree 0 -> 1. `ar` is not modified.
//   Degenerate input -- a non-finite coefficient, or polyroot() failing
//   (allocation failure, divergence) -> 0, i.e. "not stationary"; the ctsa
//   callers then report fit status 10/12 instead of continuing. Upstream
//   returned 1 for NaN/inf coefficients, but only because it compared the
//   moduli of undefined roots (uninitialised memory with upstream polyroot,
//   NaN with ours); the conservative answer is chosen deliberately.
//
// int invertroot(int q, double *ma)
//   Makes the MA part invertible, in place: every root r of theta with
//   |r| < 1 is replaced by its reflection in the unit circle, 1 / conj(r);
//   theta is rebuilt from the roots with constant term 1 and the imaginary
//   parts of the rebuilt coefficients are dropped. Coefficients beyond the
//   true degree (trailing zeros) are left untouched. Returns 1 if the
//   coefficients were rewritten, 0 if they were left as they are: no root
//   inside the circle, q <= 0, true degree 0, or degenerate input (a
//   non-finite coefficient, polyroot() failing, or a rebuilt coefficient
//   that is not finite).
//   Reflecting by 1 / conj(r) rather than 1 / r gives the same polynomial
//   whenever both members of a complex-conjugate pair are reflected, and stays
//   correct when rounding puts only one member of a pair on the inside of the
//   circle (e.g. theta(z) = 1 + z^4): upstream's result there was not
//   theta. See PROVENANCE.md.
//
// Rebuilding: theta(z) = prod_i (1 - w_i z), w_i = 1 / root_i, expanded in
// complex arithmetic. For a reflected root w_i = conj(r_i) exactly, so no
// division is needed for it.

#include <math.h>
#include <stdlib.h>

#include "polyroot.h"
#include "talg.h"

// True degree of 1 +/- c[0] z + ... + c[n-1] z^n: the index of the last
// non-zero coefficient plus one. Returns -1 if a coefficient is not finite.
static int arma_true_degree(int n, const double* c) {
  int i, m = 0;
  for (i = 0; i < n; ++i) {
    if (!isfinite(c[i])) return -1;
    if (c[i] != 0.0) m = i + 1;
  }
  return m;
}

// Roots of 1 + sign * (c[0] z + ... + c[m-1] z^m), m >= 1, c[m-1] != 0, into
// a block of 2m doubles allocated here (real parts, then imaginary parts).
// Returns NULL on allocation or polyroot() failure.
static double* arma_roots(int m, const double* c, double sign) {
  double* coeff = (double*)malloc(sizeof(double) * (size_t)(3 * m + 1));
  double* roots;
  int i;
  if (coeff == NULL) return NULL;
  roots = coeff + (m + 1);
  coeff[0] = 1.0;
  for (i = 0; i < m; ++i) coeff[i + 1] = sign * c[i];
  if (polyroot(coeff, m, roots, roots + m) != 0) {
    free(coeff);
    return NULL;
  }
  // Hand back a block that starts at the roots: move them to the front.
  for (i = 0; i < 2 * m; ++i) coeff[i] = roots[i];
  return coeff;
}

int archeck(int p, double* ar) {
  double* roots;
  int i, m, stationary;

  if (p <= 0) return 1;
  m = arma_true_degree(p, ar);
  if (m < 0) return 0;
  if (m == 0) return 1;

  roots = arma_roots(m, ar, -1.0);
  if (roots == NULL) return 0;
  stationary = 1;
  for (i = 0; i < m; ++i) {
    // `!(x >= 1)` also catches a NaN modulus.
    if (!(hypot(roots[i], roots[m + i]) >= 1.0)) {
      stationary = 0;
      break;
    }
  }
  free(roots);
  return stationary;
}

int invertroot(int q, double* ma) {
  double *roots, *wr, *wi, *cr, *ci;
  int i, k, m, reflected;

  if (q <= 0) return 0;
  m = arma_true_degree(q, ma);
  if (m <= 0) return 0;

  roots = arma_roots(m, ma, 1.0);
  if (roots == NULL) return 0;

  wr = (double*)malloc(sizeof(double) * (size_t)(4 * m + 2));
  if (wr == NULL) {
    free(roots);
    return 0;
  }
  wi = wr + m;
  cr = wi + m;
  ci = cr + (m + 1);

  // w_i = 1 / root_i, with root_i first reflected if it lies inside.
  reflected = 0;
  for (i = 0; i < m; ++i) {
    const double re = roots[i], im = roots[m + i];
    if (!(hypot(re, im) >= 1.0)) {
      if (!isfinite(re) || !isfinite(im)) break;  // undefined root: give up
      // root -> 1 / conj(root), hence w = conj(root).
      wr[i] = re;
      wi[i] = -im;
      reflected = 1;
    } else if (fabs(re) >= fabs(im)) {
      // 1 / (re + i im) by Smith's algorithm (no overflow in re^2 + im^2).
      const double t = im / re, d = re + im * t;
      wr[i] = 1.0 / d;
      wi[i] = -t / d;
    } else {
      const double t = re / im, d = re * t + im;
      wr[i] = t / d;
      wi[i] = -1.0 / d;
    }
  }
  free(roots);
  if (i < m || !reflected) {
    free(wr);
    return 0;
  }

  // c(z) = prod_i (1 - w_i z), built one factor at a time, highest power
  // first so that each step can be done in place.
  cr[0] = 1.0;
  ci[0] = 0.0;
  for (i = 0; i < m; ++i) {
    cr[i + 1] = 0.0;
    ci[i + 1] = 0.0;
    for (k = i + 1; k >= 1; --k) {
      cr[k] -= wr[i] * cr[k - 1] - wi[i] * ci[k - 1];
      ci[k] -= wr[i] * ci[k - 1] + wi[i] * cr[k - 1];
    }
  }

  for (k = 1; k <= m; ++k) {
    if (!isfinite(cr[k])) {
      free(wr);
      return 0;
    }
  }
  for (k = 1; k <= m; ++k) ma[k - 1] = cr[k];
  free(wr);
  return 1;
}
