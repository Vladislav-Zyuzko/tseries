// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// Clean-room replacement for ctsa's LGPL-3.0 `brent_local_min` (upstream
// `src/brent.c`, removed from this tree). Written from the published algorithm
// description only; the upstream body was not consulted. Process and sources:
// third_party/ctsa/PROVENANCE.md, "Clean-room replacements".
//
// Algorithm: R. P. Brent, "Algorithms for Minimization without Derivatives",
// Prentice-Hall, 1973, chapter 5 — golden-section search on a bracketing
// interval, accelerated by successive parabolic interpolation through the three
// best points seen so far, with the safeguards described there (a parabolic
// step is only taken if it lands inside the interval and is shorter than half
// of the step before last; no evaluation closer than `tol` to a previous one or
// to the interval ends).
//
// Only caller in the tree: `fminbnd()` (optimc.c), used by Box-Cox automatic
// lambda selection.

#include <float.h>
#include <math.h>

#include "brent.h"

// Upper bound on iterations. Brent's method has no cap of its own; it is
// guaranteed to terminate for t > 0, but `fminbnd` passes a relative tolerance
// of one machine epsilon, so a cap keeps a pathological objective (NaN, flat
// plateaus) from spinning. Each iteration shrinks the interval by at least the
// golden-section factor every few steps, so 500 is far beyond what a
// well-posed problem needs (~80 golden steps reach 1e-16 of the start width).
#define TSERIES_BRENT_MAX_ITER 500

double brent_local_min(custom_funcuni *funcuni, double a, double b, double t,
                       double eps, double *x) {
  // (3 - sqrt(5)) / 2: the golden-section fraction.
  const double golden = 0.38196601125010515;
  double lo = a < b ? a : b;
  double hi = a < b ? b : a;
  double xm, fxm;  // best point so far and its value
  double w, fw;    // second best
  double v, fv;    // previous value of w
  double d = 0.0;  // current step
  double e = 0.0;  // step before last
  int iter;

  // A relative tolerance below machine precision cannot be honoured; keep the
  // termination test meaningful when t == 0.
  if (!(eps >= DBL_EPSILON)) eps = DBL_EPSILON;
  if (!(t >= 0.0)) t = 0.0;

  xm = w = v = lo + golden * (hi - lo);
  fxm = fw = fv = FUNCUNI_EVAL(funcuni, xm);

  for (iter = 0; iter < TSERIES_BRENT_MAX_ITER; ++iter) {
    const double mid = 0.5 * (lo + hi);
    const double tol = eps * fabs(xm) + t;
    const double tol2 = 2.0 * tol;
    int take_golden = 1;
    double u, fu;

    // Stop when the bracket around xm is within 2*tol.
    if (fabs(xm - mid) <= tol2 - 0.5 * (hi - lo)) break;

    if (fabs(e) > tol) {
      // Parabola through (xm, fxm), (w, fw), (v, fv); step p/q from xm.
      double r = (xm - w) * (fxm - fv);
      double q = (xm - v) * (fxm - fw);
      double p = (xm - v) * q - (xm - w) * r;
      double e_prev = e;
      q = 2.0 * (q - r);
      if (q > 0.0) {
        p = -p;
      } else {
        q = -q;
      }
      e = d;
      if (fabs(p) < fabs(0.5 * q * e_prev) && p > q * (lo - xm) &&
          p < q * (hi - xm)) {
        d = p / q;
        u = xm + d;
        // Do not evaluate too close to the interval ends.
        if (u - lo < tol2 || hi - u < tol2) d = xm < mid ? tol : -tol;
        take_golden = 0;
      }
    }

    if (take_golden) {
      // Golden-section step into the larger of the two sub-intervals.
      e = xm < mid ? hi - xm : lo - xm;
      d = golden * e;
    }

    // Never evaluate closer than tol to the current best point.
    if (fabs(d) >= tol) {
      u = xm + d;
    } else {
      u = d > 0.0 ? xm + tol : xm - tol;
    }
    fu = FUNCUNI_EVAL(funcuni, u);

    if (fu <= fxm) {
      if (u < xm) {
        hi = xm;
      } else {
        lo = xm;
      }
      v = w;
      fv = fw;
      w = xm;
      fw = fxm;
      xm = u;
      fxm = fu;
    } else {
      if (u < xm) {
        lo = u;
      } else {
        hi = u;
      }
      if (fu <= fw || w == xm) {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      } else if (fu <= fv || v == xm || v == w) {
        v = u;
        fv = fu;
      }
    }
  }

  *x = xm;
  return fxm;
}
