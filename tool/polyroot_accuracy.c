// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// Accuracy harness for tseries' clean-room `polyroot()`
// (third_party/ctsa_replacements/polyroot.c), which replaced ctsa's ACM
// Algorithm 419 translation. Checks, against numpy reference roots in
// test/fixtures/polyroot_reference.txt (regenerate with
// tool/polyroot_reference.py):
//   * backward error of every computed root (independent of numpy);
//   * forward error against numpy, scaled by each root's condition number;
//   * the verdict of ctsa's archeck() (stationary or not) -- the decision
//     behind fit status codes 10/12 -- on every AR case;
//   * the coefficients produced by ctsa's invertroot() (MA root reflection);
//   * the degenerate inputs documented in polyroot.c.
// Exits non-zero if any gate fails.
//
// archeck()/invertroot() are internal to the `ctsa` code asset, so this
// harness compiles the sources directly rather than going through FFI. Since
// 2026-10-05 both live in ctsa_replacements/arma_roots.c (see also
// tool/arma_roots_accuracy.c). Build and run from packages/tseries (host
// clang shown):
//
//   clang -O2 -std=c11 -D_CRT_SECURE_NO_WARNINGS \
//     -Ithird_party/ctsa_replacements -Ithird_party/ctsa/src \
//     tool/polyroot_accuracy.c third_party/ctsa_replacements/polyroot.c \
//     third_party/ctsa_replacements/arma_roots.c \
//     -o "$TMP/polyroot_accuracy" -lm      (drop -lm on Windows)
//   "$TMP/polyroot_accuracy" test/fixtures/polyroot_reference.txt
//
// Last run: see third_party/ctsa/PROVENANCE.md, "Clean-room replacements".

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "polyroot.h"
#include "talg.h"

#define MAXN 512

static int g_fail = 0;

static void gate(const char *what, double worst, double limit) {
  const int ok = worst <= limit;
  printf("%-46s worst %.3g  gate %.3g  %s\n", what, worst, limit,
         ok ? "OK" : "FAIL");
  if (!ok) g_fail = 1;
}

// |p(z)| / sum |a_i| |z|^i, evaluated on the reversed polynomial for |z| > 1.
static double backward_error(const double *a, int n, double zr, double zi) {
  double pr, pi, s, t, wr, wi, az = hypot(zr, zi);
  int i;
  if (az <= 1.0) {
    pr = a[n];
    pi = 0.0;
    s = fabs(a[n]);
    for (i = n - 1; i >= 0; --i) {
      t = pr * zr - pi * zi + a[i];
      pi = pr * zi + pi * zr;
      pr = t;
      s = s * az + fabs(a[i]);
    }
  } else {
    wr = zr / (az * az);
    wi = -zi / (az * az);
    pr = a[0];
    pi = 0.0;
    s = fabs(a[0]);
    for (i = 1; i <= n; ++i) {
      t = pr * wr - pi * wi + a[i];
      pi = pr * wi + pi * wr;
      pr = t;
      s = s / az + fabs(a[i]);
    }
  }
  return hypot(pr, pi) / s;
}

static int read_doubles(char **cursor, double *out, int count) {
  int i;
  for (i = 0; i < count; ++i) {
    char *end;
    out[i] = strtod(*cursor, &end);
    if (end == *cursor) return 0;
    *cursor = end;
  }
  return 1;
}

static void degenerate_cases(void) {
  double zr[4], zi[4], ar[4];
  int ret, bad = 0;

  {  // zero constant term: exact zero root first
    double c[] = {0.0, 1.0, -0.5};
    ret = polyroot(c, 2, zr, zi);
    bad |= !(ret == 0 && zr[0] == 0.0 && zi[0] == 0.0 &&
             fabs(zr[1] - 2.0) < 1e-15 && zi[1] == 0.0);
  }
  {  // degree deficient: status 1 as upstream, outputs defined
    double c[] = {1.0, -0.5, 0.0, 0.0};
    ret = polyroot(c, 3, zr, zi);
    bad |= !(ret == 1 && fabs(zr[0] - 2.0) < 1e-15 && isinf(zr[1]) &&
             isinf(zr[2]));
  }
  {  // identically zero
    double c[] = {0.0, 0.0, 0.0};
    ret = polyroot(c, 2, zr, zi);
    bad |= !(ret == 1 && isnan(zr[0]) && isnan(zi[1]));
  }
  {  // non-finite coefficient
    double c[] = {1.0, NAN, 0.3};
    ret = polyroot(c, 2, zr, zi);
    bad |= !(ret == 1 && isnan(zr[0]) && isnan(zr[1]));
  }
  {  // degree 0, non-zero constant: success, nothing to write
    double c[] = {1.0};
    zr[0] = 7.0;
    ret = polyroot(c, 0, zr, zi);
    bad |= !(ret == 0 && zr[0] == 7.0);
  }
  {  // quadratic with complex roots +-2i, and with huge/tiny coefficients
    double c[] = {1.0, 0.0, 0.25};
    double d[] = {1e-300, 1.0, 1e300};
    ret = polyroot(c, 2, zr, zi);
    bad |= !(ret == 0 && zr[0] == 0.0 && fabs(fabs(zi[0]) - 2.0) < 1e-15 &&
             zi[1] == -zi[0]);
    ret = polyroot(d, 2, zr, zi);
    bad |= !(ret == 0 && fabs(zr[0] + 5e-301) < 1e-315 &&
             fabs(fabs(zi[0]) - 8.660254037844386e-301) < 1e-314);
  }
  {  // archeck() with a zero last AR coefficient (upstream: read garbage)
    ar[0] = 0.5;
    ar[1] = 0.0;  // 1 - 0.5 z: root 2, stationary
    bad |= archeck(2, ar) != 1;
    ar[0] = 1.5;
    ar[1] = 0.0;  // 1 - 1.5 z: root 2/3, not stationary
    bad |= archeck(2, ar) != 0;
  }
  printf("%-46s %s\n", "degenerate inputs (documented contract)",
         bad ? "FAIL" : "OK");
  if (bad) g_fail = 1;
}

int main(int argc, char **argv) {
  static char line[1 << 20];
  static double a[MAXN + 1], ref[3 * MAXN], zr[MAXN], zi[MAXN], ma[MAXN],
      expect[MAXN];
  static int used[MAXN], order[MAXN];
  FILE *f;
  int npoly = 0, ninv = 0, nverdict = 0, nambig = 0, nfwd = 0;
  double worst_be = 0.0, worst_fwd = 0.0, worst_inv = 0.0;
  char worst_fwd_tag[128] = "", worst_be_tag[128] = "";

  if (argc < 2) {
    fprintf(stderr, "usage: %s polyroot_reference.txt\n", argv[0]);
    return 2;
  }
  f = fopen(argv[1], "r");
  if (f == NULL) {
    perror(argv[1]);
    return 2;
  }

  while (fgets(line, sizeof line, f) != NULL) {
    char kind[16], tag[128];
    int n, consumed, i, j;
    char *cur;
    if (line[0] == '#' || line[0] == '\n') continue;
    if (sscanf(line, "%15s %127s %d%n", kind, tag, &n, &consumed) != 3 ||
        n < 1 || n > MAXN) {
      fprintf(stderr, "bad line: %.60s\n", line);
      return 2;
    }
    cur = line + consumed;

    if (strcmp(kind, "poly") == 0) {
      int ret;
      double minmod_ref = HUGE_VAL;
      if (!read_doubles(&cur, a, n + 1) || !read_doubles(&cur, ref, 3 * n)) {
        fprintf(stderr, "bad poly line %s\n", tag);
        return 2;
      }
      ret = polyroot(a, n, zr, zi);
      if (ret != 0) {
        printf("polyroot failed (%d) on %s\n", ret, tag);
        g_fail = 1;
        continue;
      }
      ++npoly;
      for (i = 0; i < n; ++i) {
        const double be =
            backward_error(a, n, zr[i], zi[i]) / ((double)(n + 1) * DBL_EPSILON);
        if (!(be <= worst_be)) {
          worst_be = be;
          snprintf(worst_be_tag, sizeof worst_be_tag, "%s", tag);
        }
        used[i] = 0;
        order[i] = i;
        minmod_ref = fmin(minmod_ref, hypot(ref[3 * i], ref[3 * i + 1]));
      }
      // Match best-conditioned reference roots first.
      for (i = 1; i < n; ++i) {
        for (j = i; j > 0 && ref[3 * order[j] + 2] < ref[3 * order[j - 1] + 2];
             --j) {
          const int t = order[j];
          order[j] = order[j - 1];
          order[j - 1] = t;
        }
      }
      for (i = 0; i < n; ++i) {
        const int k = order[i];
        const double rr = ref[3 * k], ri = ref[3 * k + 1], cond = ref[3 * k + 2];
        double best = HUGE_VAL;
        int bj = -1;
        for (j = 0; j < n; ++j) {
          double d;
          if (used[j]) continue;
          d = hypot(zr[j] - rr, zi[j] - ri);
          if (d < best) {
            best = d;
            bj = j;
          }
        }
        used[bj] = 1;
        // Forward error in units of (condition number x unit roundoff): both
        // root finders are backward stable, so this should be O(n).
        if (isfinite(cond) && cond < 1e8) {
          const double scaled =
              best / hypot(rr, ri) / (cond * DBL_EPSILON) / (double)n;
          ++nfwd;
          if (scaled > worst_fwd) {
            worst_fwd = scaled;
            snprintf(worst_fwd_tag, sizeof worst_fwd_tag, "%s", tag);
          }
        }
      }
      // archeck() verdict on AR cases (a_0 = 1, ar_k = -a_k).
      if (strncmp(tag, "ar", 2) == 0 || strncmp(tag, "sar", 3) == 0 ||
          strncmp(tag, "near", 4) == 0) {
        static double arc[MAXN];
        int want, got;
        for (i = 0; i < n; ++i) arc[i] = -a[i + 1] / a[0];
        if (fabs(minmod_ref - 1.0) <= 1e-9) {
          ++nambig;
        } else {
          want = minmod_ref >= 1.0;
          got = archeck(n, arc);
          ++nverdict;
          if (want != got) {
            printf("archeck verdict differs on %s: numpy %d, ours %d\n", tag,
                   want, got);
            g_fail = 1;
          }
        }
      }
    } else if (strcmp(kind, "inv") == 0) {
      int ret, any_inside = 0;
      if (!read_doubles(&cur, ma, n) || !read_doubles(&cur, expect, n)) {
        fprintf(stderr, "bad inv line %s\n", tag);
        return 2;
      }
      {
        double c[MAXN + 1];
        c[0] = 1.0;
        for (i = 0; i < n; ++i) c[i + 1] = ma[i];
        if (polyroot(c, n, zr, zi) == 0) {
          for (i = 0; i < n; ++i) any_inside |= hypot(zr[i], zi[i]) < 1.0;
        }
      }
      ret = invertroot(n, ma);
      ++ninv;
      if (ret != any_inside) {
        printf("invertroot status differs on %s\n", tag);
        g_fail = 1;
      }
      for (i = 0; i < n; ++i) {
        const double e = fabs(ma[i] - expect[i]) / fmax(1.0, fabs(expect[i]));
        if (e > worst_inv) worst_inv = e;
      }
    }
  }
  fclose(f);

  printf("polynomials %d, roots compared %d, archeck verdicts %d (%d within "
         "1e-9 of |z| = 1 skipped), invertroot cases %d\n",
         npoly, nfwd, nverdict, nambig, ninv);
  gate("backward error / ((n + 1) eps)", worst_be, 8.0);
  printf("  (worst on %s)\n", worst_be_tag);
  gate("forward error / (cond * eps * n)", worst_fwd, 8.0);
  printf("  (worst on %s)\n", worst_fwd_tag);
  gate("invertroot coefficients vs numpy reflection", worst_inv, 1e-12);
  degenerate_cases();
  printf("%s\n", g_fail ? "SOME GATES FAILED" : "all gates OK");
  return g_fail;
}
