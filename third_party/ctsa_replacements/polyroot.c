// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// Replacement for ctsa's `polyroot()` (upstream `src/polyroot.c`), which was a
// translation of ACM Algorithm 419 (CPOLY, Jenkins & Traub) distributed under
// the ACM licence for non-commercial use, and has been removed from the
// vendored tree. See third_party/ctsa/PROVENANCE.md, "Clean-room
// replacements". This file was written without reading upstream's code; the
// contract below was taken from the vendored header (`polyroot.h`), from the
// BSD-3-Clause call sites (`archeck()`/`invertroot()` in talg.c, `myarima()` in
// ctsa.c) and from black-box tests of the upstream binary.
//
// Contract
//   coeff   DEGREE + 1 real coefficients in INCREASING powers:
//           p(z) = coeff[0] + coeff[1] z + ... + coeff[DEGREE] z^DEGREE.
//           Not modified.
//   ZEROR/ZEROI  DEGREE entries each: real and imaginary parts of the roots.
//   Returns 0 on success and 1 on failure, as upstream:
//     * DEGREE < 0, or a non-finite coefficient, or p == 0 identically
//       -> 1; every output slot is set to NaN (upstream left the outputs
//       untouched, i.e. uninitialised memory for every caller in ctsa);
//     * coeff[DEGREE] == 0 (the polynomial has a lower degree than stated)
//       -> 1, as upstream; the outputs are nevertheless defined: the roots of
//       the polynomial of its true degree, followed by +inf for each missing
//       root (the limit of those roots as the leading coefficient -> 0).
//       Upstream's `archeck()` ignored the return value and read all DEGREE
//       roots, so with upstream it read uninitialised memory in this case
//       (the current `archeck()`, ctsa_replacements/arma_roots.c, drops
//       trailing zero coefficients first and checks the status);
//     * DEGREE == 0 with a finite non-zero constant -> 0, nothing written.
//   Roots at the origin (low-order zero coefficients) are returned exactly as
//   0 + 0i, first. The order of the remaining roots is unspecified; no caller
//   depends on it (they test moduli or rebuild the polynomial from all roots).
//
// Method
//   Degree 1 and 2: closed form (the quadratic by the cancellation-free
//   formula -- the larger root from -(b + sign(b) sqrt(disc)) / 2, the other
//   from Vieta's x1 x2 = c / a; see N. J. Higham, "Accuracy and Stability of
//   Numerical Algorithms", 2nd ed., SIAM 2002, sec. 1.8).
//   Degree >= 3: the Ehrlich-Aberth simultaneous iteration
//     z_i <- z_i - N_i / (1 - N_i * sum_{j != i} 1 / (z_i - z_j)),
//     N_i = p(z_i) / p'(z_i)
//   (L. W. Ehrlich, "A modified Newton method for polynomials", Comm. ACM 10
//   (1967) 107-108; O. Aberth, "Iteration methods for finding all zeros of a
//   polynomial simultaneously", Math. Comp. 27 (1973) 339-344), applied in
//   Gauss-Seidel fashion, with the starting points, the evaluation of the
//   reversed polynomial for |z| > 1 and the backward-error stopping rule
//   described by D. A. Bini, "Numerical computation of polynomial zeros by
//   means of Aberth's method", Numerical Algorithms 13 (1996) 179-200:
//     * starting points on circles whose radii come from the upper convex
//       hull ("Newton polygon") of the points (i, log|a_i|);
//     * a root stops moving once |p(z)| <= tol * sum_i |a_i| |z|^i, i.e. it is
//       the exact root of a polynomial whose coefficients differ from the
//       given ones by a relative amount of order the unit roundoff.
//   Complex arithmetic is written out by hand: MSVC's C compiler has no C99
//   _Complex.

#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "polyroot.h"

typedef struct {
  double re, im;
} pr_cplx;

static pr_cplx pr_make(double re, double im) {
  pr_cplx z;
  z.re = re;
  z.im = im;
  return z;
}

static pr_cplx pr_add(pr_cplx a, pr_cplx b) {
  return pr_make(a.re + b.re, a.im + b.im);
}

static pr_cplx pr_sub(pr_cplx a, pr_cplx b) {
  return pr_make(a.re - b.re, a.im - b.im);
}

static pr_cplx pr_mul(pr_cplx a, pr_cplx b) {
  return pr_make(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}

static pr_cplx pr_scale(pr_cplx a, double s) {
  return pr_make(a.re * s, a.im * s);
}

// a / b by Smith's algorithm (R. L. Smith, "Algorithm 116: Complex division",
// Comm. ACM 5 (1962) 435), which avoids overflow in |b|^2. b != 0.
static pr_cplx pr_div(pr_cplx a, pr_cplx b) {
  double r, d;
  if (fabs(b.re) >= fabs(b.im)) {
    r = b.im / b.re;
    d = b.re + b.im * r;
    return pr_make((a.re + a.im * r) / d, (a.im - a.re * r) / d);
  }
  r = b.re / b.im;
  d = b.re * r + b.im;
  return pr_make((a.re * r + a.im) / d, (a.im * r - a.re) / d);
}

static double pr_abs(pr_cplx a) { return hypot(a.re, a.im); }

static int pr_is_zero(pr_cplx a) { return a.re == 0.0 && a.im == 0.0; }

// Evaluates, at z, the Newton correction N = p(z) / p'(z) of the polynomial
// b[0] + b[1] z + ... + b[m] z^m (b[0] != 0, b[m] != 0) and reports whether z
// already satisfies the backward-error test. For |z| > 1 the reversed
// polynomial is evaluated at w = 1/z so that nothing overflows:
//   p(z) = z^m q(w),  q(w) = b[m] + b[m-1] w + ... + b[0] w^m,
//   p'(z) / p(z) = m w - w^2 q'(w) / q(w).
// Returns 1 if z is accepted as a root (and then *corr is 0), 0 otherwise.
static int pr_newton(const double* b, int m, pr_cplx z, double tol,
                     pr_cplx* corr) {
  pr_cplx p, dp, w;
  double s, az;
  int i;

  az = pr_abs(z);
  if (az <= 1.0) {
    p = pr_make(b[m], 0.0);
    dp = pr_make(0.0, 0.0);
    s = fabs(b[m]);
    for (i = m - 1; i >= 0; --i) {
      dp = pr_add(pr_mul(dp, z), p);
      p = pr_add(pr_mul(p, z), pr_make(b[i], 0.0));
      s = s * az + fabs(b[i]);
    }
    if (pr_abs(p) <= tol * s) {
      *corr = pr_make(0.0, 0.0);
      return 1;
    }
    if (pr_is_zero(dp)) {
      // Stationary point that is not a root: nudge instead of dividing by 0.
      *corr = pr_make(-(az + 1.0) * 1e-3, (az + 1.0) * 1e-3);
      return 0;
    }
    *corr = pr_div(p, dp);
    return 0;
  }

  w = pr_div(pr_make(1.0, 0.0), z);
  {
    const double aw = 1.0 / az;
    pr_cplx q = pr_make(b[0], 0.0);
    pr_cplx dq = pr_make(0.0, 0.0);
    pr_cplx ratio, den;
    s = fabs(b[0]);
    for (i = 1; i <= m; ++i) {
      dq = pr_add(pr_mul(dq, w), q);
      q = pr_add(pr_mul(q, w), pr_make(b[i], 0.0));
      s = s * aw + fabs(b[i]);
    }
    if (pr_abs(q) <= tol * s) {
      *corr = pr_make(0.0, 0.0);
      return 1;
    }
    ratio = pr_div(dq, q);
    den = pr_sub(pr_scale(w, (double)m), pr_mul(pr_mul(w, w), ratio));
    if (pr_is_zero(den)) {
      *corr = pr_make(-(az + 1.0) * 1e-3, (az + 1.0) * 1e-3);
      return 0;
    }
    *corr = pr_div(pr_make(1.0, 0.0), den);
    return 0;
  }
}

// Starting points (Bini 1996, sec. 3): upper convex hull of (i, log|b_i|),
// i = 0..m, over the non-zero coefficients; each hull edge from k to l gives
// l - k points on the circle of radius (|b_k| / |b_l|)^(1 / (l - k)).
// Returns 0, or -1 on allocation failure.
static int pr_start(const double* b, int m, pr_cplx* z) {
  const double sigma = 0.7;  // angular offset, breaks real-axis symmetry
  const double two_pi = 6.283185307179586476925286766559;
  int* hull = (int*)malloc(sizeof(int) * (size_t)(m + 1));
  int nh = 0, i, t, j, pos = 0;
  if (hull == NULL) return -1;

  for (i = 0; i <= m; ++i) {
    if (b[i] == 0.0) continue;
    // Pop while the last two hull points and i do not make a right turn
    // (i.e. the middle point lies on or below the chord).
    while (nh >= 2) {
      const int i1 = hull[nh - 2], i2 = hull[nh - 1];
      const double y1 = log(fabs(b[i1])), y2 = log(fabs(b[i2]));
      const double y3 = log(fabs(b[i]));
      const double cross =
          (double)(i2 - i1) * (y3 - y1) - (y2 - y1) * (double)(i - i1);
      if (cross >= 0.0) {
        --nh;
      } else {
        break;
      }
    }
    hull[nh++] = i;
  }

  for (t = 0; t + 1 < nh; ++t) {
    const int k = hull[t], l = hull[t + 1], cnt = l - k;
    const double r =
        exp((log(fabs(b[k])) - log(fabs(b[l]))) / (double)cnt);
    for (j = 0; j < cnt; ++j) {
      const double ang =
          two_pi * (double)j / (double)cnt + two_pi * (double)k / (double)m +
          sigma;
      z[pos++] = pr_make(r * cos(ang), r * sin(ang));
    }
  }
  free(hull);
  return pos == m ? 0 : -1;
}

// Roots of b[0] + ... + b[m] z^m, m >= 3, b[0] != 0, b[m] != 0, all finite,
// already scaled so that max|b_i| == 1. Returns 0, or 1 on failure.
static int pr_aberth(const double* b, int m, pr_cplx* z) {
  const int max_sweeps = 1000;
  // Running-error bound of Horner's rule is about 2 m u sum|b_i||z|^i
  // (Higham 2002, eq. 5.3); accept a root once the residual is at that level.
  const double tol = 4.0 * (double)(m + 1) * DBL_EPSILON;
  unsigned char* done;
  int i, j, sweep, remaining;

  if (pr_start(b, m, z) != 0) return 1;
  done = (unsigned char*)calloc((size_t)m, 1);
  if (done == NULL) return 1;

  remaining = m;
  for (sweep = 0; sweep < max_sweeps && remaining > 0; ++sweep) {
    for (i = 0; i < m; ++i) {
      pr_cplx corr, sum, den;
      if (done[i]) continue;
      if (pr_newton(b, m, z[i], tol, &corr)) {
        done[i] = 1;
        --remaining;
        continue;
      }
      sum = pr_make(0.0, 0.0);
      for (j = 0; j < m; ++j) {
        pr_cplx diff;
        if (j == i) continue;
        diff = pr_sub(z[i], z[j]);
        if (pr_is_zero(diff)) continue;  // coincident iterates: skip the pole
        sum = pr_add(sum, pr_div(pr_make(1.0, 0.0), diff));
      }
      den = pr_sub(pr_make(1.0, 0.0), pr_mul(corr, sum));
      if (!pr_is_zero(den)) corr = pr_div(corr, den);
      z[i] = pr_sub(z[i], corr);
      if (!isfinite(z[i].re) || !isfinite(z[i].im)) {
        free(done);
        return 1;
      }
    }
  }
  free(done);
  // Not every root met the backward-error test within the sweep budget. This
  // happens for clusters of (nearly) multiple roots, whose iterates are then
  // as accurate as the conditioning allows; they are returned as they are,
  // which matches what a root finder can deliver for such input.
  return 0;
}

// Roots of b0 + b1 z + b2 z^2 (b0 != 0, b2 != 0, all finite).
//
// The variable is scaled first, z = lambda y with lambda = sqrt|b0 / b2| (the
// geometric mean of the root moduli), which turns the equation into
//   1 + B y + s y^2 = 0,  B = b1 lambda / b0,  s = sign(b0 b2),
// so that nothing under- or overflows unless a root itself is out of range
// (coefficient scaling alone cannot do this: 1e-300 + z + 1e300 z^2 loses
// its constant term to underflow).
static void pr_quadratic(double b0, double b1, double b2, pr_cplx* z) {
  const double lambda = sqrt(fabs(b0)) / sqrt(fabs(b2));
  const double s = ((b0 > 0.0) == (b2 > 0.0)) ? 1.0 : -1.0;
  const double B = copysign(1.0, b0) * (b1 / sqrt(fabs(b0)) / sqrt(fabs(b2)));
  const double aB = fabs(B);
  if (s < 0.0 || aB >= 2.0) {
    // Real roots: disc = B^2 - 4 s >= 0, evaluated without overflowing B^2.
    const double root_disc =
        aB > 1e150 ? aB * sqrt(1.0 - 4.0 * s / aB / aB)
                   : sqrt(B * B - 4.0 * s);
    // q != 0: |q| >= |B| / 2, and B == 0 only with s = -1, where q = -1.
    const double q = -0.5 * (B + copysign(root_disc, B));
    z[0] = pr_make(lambda * (q * s), 0.0);  // q / s
    z[1] = pr_make(lambda * (1.0 / q), 0.0);
  } else {
    // s = +1, |B| < 2: complex pair on the circle |y| = 1.
    const double re = -0.5 * B;
    const double im = 0.5 * sqrt((2.0 - aB) * (2.0 + aB));
    z[0] = pr_make(lambda * re, lambda * im);
    z[1] = pr_make(lambda * re, -lambda * im);
  }
}

int polyroot(double* coeff, int DEGREE, double* ZEROR, double* ZEROI) {
  int i, top, low, m, status, pos;
  double* b = NULL;
  pr_cplx* z = NULL;
  double bmax;

  if (DEGREE < 0 || coeff == NULL) return 1;
  for (i = 0; i <= DEGREE; ++i) {
    if (!isfinite(coeff[i])) goto undefined;
  }

  top = -1;
  for (i = DEGREE; i >= 0; --i) {
    if (coeff[i] != 0.0) {
      top = i;
      break;
    }
  }
  if (top < 0) goto undefined;  // p == 0: every z is a root
  if (DEGREE == 0) return 0;

  low = 0;
  while (coeff[low] == 0.0) ++low;  // terminates: coeff[top] != 0
  m = top - low;

  pos = 0;
  for (i = 0; i < low; ++i, ++pos) {
    ZEROR[pos] = 0.0;
    ZEROI[pos] = 0.0;
  }

  status = 0;
  if (m == 1) {
    ZEROR[pos] = -coeff[low] / coeff[low + 1];
    ZEROI[pos] = 0.0;
    ++pos;
  } else if (m == 2) {
    pr_cplx q[2];
    pr_quadratic(coeff[low], coeff[low + 1], coeff[low + 2], q);
    for (i = 0; i < 2; ++i, ++pos) {
      ZEROR[pos] = q[i].re;
      ZEROI[pos] = q[i].im;
    }
  } else if (m >= 3) {
    b = (double*)malloc(sizeof(double) * (size_t)(m + 1));
    z = (pr_cplx*)malloc(sizeof(pr_cplx) * (size_t)m);
    if (b == NULL || z == NULL) {
      free(b);
      free(z);
      goto undefined;
    }
    bmax = 0.0;
    for (i = 0; i <= m; ++i) bmax = fmax(bmax, fabs(coeff[low + i]));
    for (i = 0; i <= m; ++i) b[i] = coeff[low + i] / bmax;
    if (pr_aberth(b, m, z) != 0) {
      free(b);
      free(z);
      goto undefined;
    }
    for (i = 0; i < m; ++i, ++pos) {
      ZEROR[pos] = z[i].re;
      ZEROI[pos] = z[i].im;
    }
    free(b);
    free(z);
  }

  // Missing roots of a degree-deficient polynomial are at infinity.
  for (; pos < DEGREE; ++pos) {
    ZEROR[pos] = HUGE_VAL;
    ZEROI[pos] = 0.0;
  }
  if (top < DEGREE) status = 1;
  return status;

undefined:
  for (i = 0; i < DEGREE; ++i) {
    ZEROR[i] = NAN;
    ZEROI[i] = NAN;
  }
  return 1;
}
