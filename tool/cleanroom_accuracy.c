// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// Accuracy harness for tseries' clean-room replacements of ctsa's former
// LGPL units (third_party/ctsa_replacements/). Compares against scipy
// reference values in test/fixtures/cleanroom_reference.txt (regenerate with
// tool/cleanroom_reference.py) and exits non-zero if any gate is exceeded.
//
// The functions under test are internal to the `ctsa` code asset (hidden
// visibility, not part of the Dart API), so they are compiled here directly
// from source rather than reached through FFI.
//
// Build and run from packages/tseries (any C99 compiler; host clang shown):
//
//   clang -O2 -std=c11 -D_CRT_SECURE_NO_WARNINGS \
//     -Ithird_party/ctsa_replacements -Ithird_party/ctsa/src \
//     tool/cleanroom_accuracy.c third_party/ctsa_replacements/*.c \
//     third_party/ctsa/src/dist.c third_party/ctsa/src/pdist.c \
//     third_party/ctsa/src/lls.c third_party/ctsa/src/matrix.c \
//     -o "$TMP/cleanroom_accuracy" -lm      (drop -lm on Windows)
//   "$TMP/cleanroom_accuracy" test/fixtures/cleanroom_reference.txt
//
// Exit status 0 = every gate passed. Last run (2026-10-05, clang 22.1.7,
// Windows x64): erfinv/erfcinv <= 4.7e-16 rel., betainv vs scipy <= 1.7e-12
// rel., tinv <= 4.3e-13 rel., all gates OK.

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <stdlib.h>

#include "brent.h"
#include "lls.h"
#include "pdist.h"

// --- brent test functions (ids match tool/cleanroom_reference.py) ----------
static int g_evals;
static double brent_f(double x, void *params) {
  const int id = *(const int *)params;
  ++g_evals;
  switch (id) {
    case 0: return (x - 2.0) * (x - 2.0);
    case 1: return cos(x);
    case 2: return x * x * x * x - 2.0 * x;
    case 3: return -x * exp(-x);
    case 4: return fabs(x - 0.3);
    case 5: return exp(x) - 3.0 * x;
    case 6: return x;
    case 7: return (x - 1e-3) * (x - 1e-3);
    case 8: return x * log(x);
    default: return NAN;
  }
}

// Backward error of x as a root of the *computed* residual
// r(y) = ibeta(y, a, b) - alpha, expressed as a relative change of x:
// |r(x)| / (pdf(x) * x), with the residual floored at one ulp of alpha (below
// that only ibeta's rounding noise is visible: ctsa's ibeta is a staircase at
// ulp level, with plateaus of ~20 ulps of x observed). Independent of ibeta's
// own error versus the true function, which the scipy comparison measures.
static double backward_error(double x, double alpha, double a, double b) {
  const double ulp_alpha = nextafter(alpha, 1.0) - alpha;
  double r = fabs(ibeta(x, a, b) - alpha);
  const double pdf = exp((a - 1.0) * log(x) + (b - 1.0) * log1p(-x) -
                         beta_log(a, b));
  r = r > ulp_alpha ? r - ulp_alpha : 0.0;
  return r / (pdf * x);
}

static int edge(int ok, const char *what) {
  if (!ok) printf("edge case FAILED: %s\n", what);
  return !ok;
}

typedef struct {
  const char *name;
  double gate;      // max allowed relative error
  double worst;     // worst relative error seen
  char worst_case[160];
  int n;
} stat_t;

static void record(stat_t *s, double got, double ref, const char *what) {
  double err;
  if (isinf(ref) || isinf(got)) {
    err = (got == ref) ? 0.0 : INFINITY;
  } else if (ref == 0.0) {
    err = fabs(got);
  } else {
    err = fabs(got - ref) / fabs(ref);
  }
  if (isnan(got) || isnan(err)) err = INFINITY;
  ++s->n;
  if (err > s->worst) {
    s->worst = err;
    snprintf(s->worst_case, sizeof s->worst_case, "%s got=%.17g ref=%.17g",
             what, got, ref);
  }
}

int main(int argc, char **argv) {
  // Gates: see PROVENANCE.md "Clean-room replacements" for the reasoning.
  stat_t s_erfinv = {"erfinv", 1e-14, 0, "", 0};
  stat_t s_erfcinv = {"erfcinv", 1e-14, 0, "", 0};
  stat_t s_betainv = {"betainv vs scipy, alpha in [1e-12, 1-1e-6]", 1e-10, 0, "", 0};
  stat_t s_betainv_root = {"betainv backward error |r(x)|/(pdf(x) x)", 1e-13, 0, "", 0};
  stat_t s_ibeta = {"ibeta (ctsa, CC0; diagnostic only)", 1.0, 0, "", 0};
  stat_t s_tinv = {"tinv vs |scipy t.ppf| (tinv returns |t|)", 1e-10, 0, "", 0};
  stat_t s_lls = {"lls_svd2 (minfit) vs numpy: |dx|/(max|x| cond eps)", 100, 0, "", 0};
  stat_t s_brent = {"brent_local_min argmin (rel. to max(1,|x*|))", 1e-7, 0, "", 0};
  stat_t *all[] = {&s_erfinv, &s_erfcinv, &s_betainv, &s_betainv_root,
                   &s_ibeta, &s_tinv, &s_lls, &s_brent};
  static char line[65536];
  char kind[32];
  int failed = 0, max_evals = 0;
  size_t i;
  FILE *f;

  if (argc < 2) {
    fprintf(stderr, "usage: %s cleanroom_reference.txt\n", argv[0]);
    return 2;
  }
  f = fopen(argv[1], "r");
  if (!f) {
    perror(argv[1]);
    return 2;
  }
  while (fgets(line, sizeof line, f)) {
    double a1, a2, a3, ref;
    char what[128];
    if (line[0] == '#' || sscanf(line, "%31s", kind) != 1) continue;
    if (!strcmp(kind, "erfinv") &&
        sscanf(line, "%*s %lf %lf", &a1, &ref) == 2) {
      snprintf(what, sizeof what, "x=%.17g", a1);
      record(&s_erfinv, erfinv(a1), ref, what);
    } else if (!strcmp(kind, "erfcinv") &&
               sscanf(line, "%*s %lf %lf", &a1, &ref) == 2) {
      snprintf(what, sizeof what, "c=%.17g", a1);
      record(&s_erfcinv, erfcinv(a1), ref, what);
    } else if (!strcmp(kind, "betainv") &&
               sscanf(line, "%*s %lf %lf %lf %lf", &a1, &a2, &a3, &ref) == 4) {
      const double x = betainv(a1, a2, a3);
      snprintf(what, sizeof what, "alpha=%.17g a=%g b=%g", a1, a2, a3);
      record(&s_betainv, x, ref, what);
      // How well we invert ctsa's own ibeta, independent of ibeta's error.
      record(&s_betainv_root, backward_error(x, a1, a2, a3), 0.0, what);
    } else if (!strcmp(kind, "ibeta") &&
               sscanf(line, "%*s %lf %lf %lf %lf", &a1, &a2, &a3, &ref) == 4) {
      snprintf(what, sizeof what, "x=%.17g a=%g b=%g", a1, a2, a3);
      record(&s_ibeta, ibeta(a1, a2, a3), ref, what);
    } else if (!strcmp(kind, "tinv") &&
               sscanf(line, "%*s %lf %lf %lf", &a1, &a2, &ref) == 3) {
      snprintf(what, sizeof what, "p=%.17g df=%g", a1, a2);
      // ctsa's tinv() returns the two-sided magnitude for df <= 1000 (its
      // `sign` is computed but never applied; regress() relies on that).
      record(&s_tinv, fabs(tinv(a1, (int)a2)), fabs(ref), what);
    } else if (!strcmp(kind, "lls")) {
      // lls M N cond A[M*N] b[M] x[N]. Error is judged in units of
      // cond(A) * eps: a backward-stable solver cannot do better than that.
      char *p = line + 3, *end;
      const int m = (int)strtol(p, &end, 10);
      const int n = (int)strtol(end, &end, 10);
      double *A = malloc(sizeof(double) * m * n), *bv = malloc(sizeof(double) * m);
      double *xr = malloc(sizeof(double) * n), *xo = malloc(sizeof(double) * n);
      double xmax = 0.0;
      int k;
      double cond;
      p = end;
      cond = strtod(p, &p);
      for (k = 0; k < m * n; ++k) A[k] = strtod(p, &p);
      for (k = 0; k < m; ++k) bv[k] = strtod(p, &p);
      for (k = 0; k < n; ++k) {
        xr[k] = strtod(p, &p);
        if (fabs(xr[k]) > xmax) xmax = fabs(xr[k]);
      }
      lls_svd2(A, bv, m, n, xo);
      for (k = 0; k < n; ++k) {
        snprintf(what, sizeof what, "M=%d N=%d x[%d]", m, n, k);
        record(&s_lls, (xo[k] - xr[k]) / (xmax * cond * DBL_EPSILON), 0.0, what);
      }
      free(A);
      free(bv);
      free(xr);
      free(xo);
    } else if (!strcmp(kind, "brent") &&
               sscanf(line, "%*s %lf %lf %lf %lf", &a1, &a2, &a3, &ref) == 4) {
      int id = (int)a1;
      double x;
      custom_funcuni fu = {brent_f, &id};
      g_evals = 0;
      // Same tolerances fminbnd() (optimc.c) passes.
      brent_local_min(&fu, a2, a3, 1e-12, DBL_EPSILON, &x);
      if (g_evals > max_evals) max_evals = g_evals;
      snprintf(what, sizeof what, "f%d on [%g, %g] (%d evals)", id, a2, a3,
               g_evals);
      {
        const double scale = fabs(ref) > 1.0 ? fabs(ref) : 1.0;
        record(&s_brent, (x - ref) / scale, 0.0, what);  // |x - x*| / scale
      }
    }
  }
  fclose(f);

  // Edge cases (domain ends and invalid input).
  {
    int edge_fail = 0;
#define EDGE(cond) edge_fail |= edge(cond, #cond)
    EDGE(erfinv(1.0) == INFINITY && erfinv(-1.0) == -INFINITY);
    EDGE(isnan(erfinv(1.5)) && isnan(erfinv(-1.0000001)) && isnan(erfinv(NAN)));
    EDGE(erfinv(0.0) == 0.0 && erfinv(-0.0) == 0.0 && signbit(erfinv(-0.0)));
    EDGE(erfcinv(0.0) == INFINITY && erfcinv(2.0) == -INFINITY);
    EDGE(erfcinv(1.0) == 0.0 && isnan(erfcinv(-0.1)) && isnan(erfcinv(2.1)));
    EDGE(betainv(0.0, 2.0, 3.0) == 0.0 && betainv(1.0, 2.0, 3.0) == 1.0);
    EDGE(isnan(betainv(0.5, -1.0, 3.0)) && isnan(betainv(0.5, 2.0, 0.0)));
    EDGE(isnan(betainv(NAN, 2.0, 3.0)));
    EDGE(fabs(betainv(0.5, 3.0, 3.0) - 0.5) < 1e-15);  // symmetric shape
    EDGE(r8_max(1.0, 2.0) == 2.0 && r8_max(-1.0, -2.0) == -1.0);
#undef EDGE
    printf("edge cases: %s\n", edge_fail ? "FAIL" : "OK");
    if (edge_fail) failed = 1;
  }

  for (i = 0; i < sizeof all / sizeof all[0]; ++i) {
    const stat_t *s = all[i];
    const int ok = s->worst <= s->gate;
    printf("%-55s n=%4d  worst rel. err %.3g  (gate %.0e) %s\n    worst: %s\n",
           s->name, s->n, s->worst, s->gate, ok ? "OK" : "FAIL",
           s->worst_case);
    if (!ok || s->n == 0) failed = 1;
  }
  printf("brent_local_min: max evaluations per call = %d\n", max_evals);
  return failed;
}
