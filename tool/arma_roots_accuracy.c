// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// Differential harness for tseries' clean-room `archeck()` / `invertroot()`
// (third_party/ctsa_replacements/arma_roots.c), which replaced ctsa's ports of
// R's arCheck()/maInvert() (GPL-2+). See third_party/ctsa/PROVENANCE.md,
// "Clean-room replacements", round 3.
//
// A deterministic set of 20 000 polynomials (families below) is generated from
// a fixed seed using only IEEE-exact operations (+ - * / sqrt; build with
// FP contraction off), so every platform sees bit-identical inputs. Each
// polynomial c(z) = 1 + c_1 z + ... + c_n z^n is fed to
//   archeck(n, ar)     with ar_k = -c_k  (phi(z) = c(z)), and
//   invertroot(n, ma)  with ma_k =  c_k  (theta(z) = c(z)).
//
// Two modes.
//
// * Default: compares the current archeck()/invertroot() with the outputs of
//   the upstream code recorded in test/fixtures/arma_roots_reference.txt
//   (every archeck verdict and invertroot status; the invertroot coefficients
//   of every 20th case, %.17g) and runs independent checks on all cases:
//   after invertroot() no root of theta lies inside the unit circle, and for
//   polynomials built from known roots the result matches the reflection
//   computed from those roots. Build from packages/tseries:
//
//     clang -O2 -std=c11 -ffp-contract=off -D_CRT_SECURE_NO_WARNINGS \
//       -Ithird_party/ctsa_replacements -Ithird_party/ctsa/src \
//       tool/arma_roots_accuracy.c third_party/ctsa_replacements/arma_roots.c \
//       third_party/ctsa_replacements/polyroot.c -o "$TMP/arma_roots_accuracy"
//     "$TMP/arma_roots_accuracy" test/fixtures/arma_roots_reference.txt
//
// * -DARMA_ROOTS_RECORD: links, in addition, the upstream functions renamed to
//   old_archeck()/old_invertroot() (compile upstream talg.c from git history,
//   `git show 7a5bf14:packages/tseries/third_party/ctsa/src/talg.c`, with
//   -Darcheck=old_archeck -Dinvertroot=old_invertroot and link it with the
//   rest of the vendored sources), compares old and new on every case with
//   full coefficients, and writes the fixture to stdout.
//
// Divergence policy (what the gates allow):
//   * cases with a root within the noise band of |z| = 1 (1e-6 for simple
//     roots, 1e-3 when the family has multiple roots): the verdict of either
//     version is a rounding accident, the outputs are only counted;
//   * non-finite coefficients, or polyroot() failing on the input: archeck()
//     now returns 0 where upstream returned 1 (documented in arma_roots.c);
//   * everything else: identical archeck verdicts and invertroot statuses;
//     invertroot coefficients within 1e-12 relative to max|c_k| or, where
//     they differ by more, within 4x the change that a +-4 ulp perturbation
//     of the input causes (the problem's own conditioning); no root left
//     inside the circle; never more than 2x further than upstream from the
//     reflection of the exact roots, where those are known.
// Exits non-zero if any gate fails.

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "polyroot.h"
#include "talg.h"

#pragma STDC FP_CONTRACT OFF

#define NCASES 20000
#define MAXN 32
#define COEF_EVERY 20

#ifdef ARMA_ROOTS_RECORD
int old_archeck(int p, double *ar);
int old_invertroot(int q, double *ma);
#endif

typedef struct {
  int n;                 // stated degree (array length)
  double c[MAXN + 1];    // c[0] = 1
  int known;             // roots below are exact (as doubles)
  int m;                 // number of known roots (true degree)
  double rr[MAXN], ri[MAXN];
  double band;           // noise band around |z| = 1
  int degenerate;        // non-finite coefficient
  const char *family;
} arma_case;

// ---- deterministic generator -------------------------------------------

static uint64_t g_state;

static uint64_t next_u64(void) {  // splitmix64
  uint64_t z = (g_state += 0x9E3779B97F4A7C15ull);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  return z ^ (z >> 31);
}

static double unif(void) { return (double)(next_u64() >> 11) * 0x1.0p-53; }
static int irange(int lo, int hi) { return lo + (int)(next_u64() % (uint64_t)(hi - lo + 1)); }

static const double k_pow10[] = {1e-1, 1e-2, 1e-3, 1e-4,  1e-5,  1e-6,  1e-7, 1e-8,
                                 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15};

// A modulus drawn from [lo, hi] on a coarse log scale built from exact ops.
static double modulus_between(double lo, double hi) {
  return lo + (hi - lo) * unif() * unif();
}

// Modulus near 1: 1 +- 10^-k, or exactly 1.
static double modulus_near_one(void) {
  const int k = irange(0, 15);
  if (k == 15) return 1.0;
  return unif() < 0.5 ? 1.0 - k_pow10[k] : 1.0 + k_pow10[k];
}

// Unit complex number from a random point (no trigonometry).
static void unit_dir(double *x, double *y) {
  double a, b, r;
  do {
    a = 2.0 * unif() - 1.0;
    b = 2.0 * unif() - 1.0;
    r = sqrt(a * a + b * b);
  } while (r < 0.05 || r > 1.0);
  *x = a / r;
  *y = b / r;
}

// Adds a real root or a conjugate pair with the given modulus.
static void add_root(arma_case *k, double mod, int pair) {
  if (pair && k->m + 2 <= MAXN) {
    double x, y;
    unit_dir(&x, &y);
    if (fabs(y) < 1e-3) y = 0.5;
    k->rr[k->m] = mod * x;
    k->ri[k->m] = mod * y;
    k->rr[k->m + 1] = mod * x;
    k->ri[k->m + 1] = -mod * y;
    k->m += 2;
  } else {
    k->rr[k->m] = unif() < 0.5 ? -mod : mod;
    k->ri[k->m] = 0.0;
    k->m += 1;
  }
}

// c(z) = prod (1 - z / r_i), real parts.
static void expand_roots(arma_case *k) {
  double cr[MAXN + 1], ci[MAXN + 1];
  int i, j;
  cr[0] = 1.0;
  ci[0] = 0.0;
  for (i = 0; i < k->m; ++i) {
    const double d = k->rr[i] * k->rr[i] + k->ri[i] * k->ri[i];
    const double wr = k->rr[i] / d, wi = -k->ri[i] / d;
    cr[i + 1] = 0.0;
    ci[i + 1] = 0.0;
    for (j = i + 1; j >= 1; --j) {
      cr[j] -= wr * cr[j - 1] - wi * ci[j - 1];
      ci[j] -= wr * ci[j - 1] + wi * cr[j - 1];
    }
  }
  for (i = 0; i <= k->m; ++i) k->c[i] = cr[i];
  for (i = k->m + 1; i <= MAXN; ++i) k->c[i] = 0.0;
  k->n = k->m;
}

static void roots_case(arma_case *k, int deg, double lo, double hi) {
  k->m = 0;
  while (k->m < deg) add_root(k, modulus_between(lo, hi), unif() < 0.6 && k->m + 2 <= deg);
  k->known = 1;
  expand_roots(k);
}

static void make_case(int id, arma_case *k) {
  int i, deg;
  memset(k, 0, sizeof *k);
  g_state = 0x5EEDull * 1000003ull + (uint64_t)id;
  k->band = 1e-6;
  k->c[0] = 1.0;

  if (id < 4000) {  // A: random coefficients
    static const double scales[] = {0.1, 0.3, 0.5, 1.0, 2.0};
    const double s = scales[irange(0, 4)];
    k->family = "random";
    k->n = irange(1, 24);
    for (i = 1; i <= k->n; ++i) k->c[i] = s * (2.0 * unif() - 1.0) / (double)i;
    if (k->c[k->n] == 0.0) k->c[k->n] = 0.25;
  } else if (id < 8000) {  // B: known roots, moduli 0.3 .. 3
    k->family = "roots";
    roots_case(k, irange(1, 24), 0.3, 3.0);
  } else if (id < 12000) {  // C: roots near the unit circle
    k->family = "near-circle";
    deg = irange(1, 20);
    k->m = 0;
    while (k->m < deg) {
      const double mod = unif() < 0.5 ? modulus_near_one() : modulus_between(0.4, 2.5);
      add_root(k, mod, unif() < 0.6 && k->m + 2 <= deg);
    }
    k->known = 1;
    expand_roots(k);
  } else if (id < 14000) {  // D: multiple roots
    const int mult = irange(2, 4);
    double mod, x, y;
    k->family = "multiple";
    k->band = 1e-3;
    k->m = 0;
    mod = unif() < 0.3 ? modulus_near_one() : modulus_between(0.4, 2.5);
    if (unif() < 0.5) {
      unit_dir(&x, &y);
      if (fabs(y) < 1e-3) y = 0.5;
      for (i = 0; i < mult && k->m + 2 <= MAXN; ++i) {
        k->rr[k->m] = mod * x; k->ri[k->m] = mod * y;
        k->rr[k->m + 1] = mod * x; k->ri[k->m + 1] = -mod * y;
        k->m += 2;
      }
    } else {
      const double r = unif() < 0.5 ? -mod : mod;
      for (i = 0; i < mult; ++i) { k->rr[k->m] = r; k->ri[k->m] = 0.0; ++k->m; }
    }
    deg = irange(0, 6);
    while (deg-- > 0) add_root(k, modulus_between(0.4, 2.5), 0);
    k->known = 1;
    expand_roots(k);
  } else if (id < 16000) {  // E: trailing zeros
    const int extra = irange(1, 6);
    k->family = "tail-zeros";
    roots_case(k, irange(1, 16), 0.3, 3.0);
    k->n = k->m + extra;  // c[m+1..n] already 0
  } else if (id < 19000) {  // F: seasonal, sparse
    static const int periods[] = {2, 3, 4, 6, 7, 12, 24};
    const int s = periods[irange(0, 6)];
    const int P = irange(1, 2), p = irange(0, 2);
    double seas[MAXN + 1] = {0}, full[MAXN + 1] = {0};
    double Phi[2], phi[2];
    int a, b;
    k->family = "seasonal";
    for (i = 0; i < 2; ++i) {
      Phi[i] = unif() < 0.25 ? (unif() < 0.5 ? -1.0 : 1.0) * modulus_near_one()
                             : 1.6 * unif() - 0.8;
      phi[i] = 1.8 * unif() - 0.9;
    }
    // (1 + Phi1 z^s + Phi2 z^2s)(1 + phi1 z + phi2 z^2), truncated to MAXN.
    seas[0] = 1.0;
    seas[s] = Phi[0];
    if (P == 2 && 2 * s <= MAXN) seas[2 * s] = Phi[1] * 0.5;
    for (a = 0; a <= MAXN; ++a) {
      if (seas[a] == 0.0) continue;
      for (b = 0; b <= p && a + b <= MAXN; ++b)
        full[a + b] += seas[a] * (b == 0 ? 1.0 : phi[b - 1]);
    }
    k->n = 0;
    for (i = 0; i <= MAXN; ++i) {
      k->c[i] = full[i];
      if (full[i] != 0.0) k->n = i;
    }
    if (k->n == 0) { k->n = 1; k->c[1] = 0.5; }
  } else {  // G: degenerate and extreme inputs
    const int kind = id % 8;
    k->family = "degenerate";
    k->n = irange(1, 8);
    for (i = 1; i <= k->n; ++i) k->c[i] = 2.0 * unif() - 1.0;
    switch (kind) {
      case 0: k->c[irange(1, k->n)] = NAN; k->degenerate = 1; break;
      case 1: k->c[irange(1, k->n)] = (unif() < 0.5 ? -1 : 1) * HUGE_VAL; k->degenerate = 1; break;
      case 2: for (i = 1; i <= k->n; ++i) k->c[i] = 0.0; break;  // degree 0
      case 3: for (i = 1; i <= k->n; ++i) k->c[i] *= 1e150; break;
      case 4: for (i = 1; i <= k->n; ++i) k->c[i] *= 1e-150; break;
      case 5: k->c[k->n] = 1e-300; break;
      case 6: k->c[1] = 1e8; break;
      default: k->c[k->n] = 0.0; if (k->n == 1) k->c[1] = 0.0; break;
    }
  }
}

// ---- analysis -----------------------------------------------------------

// Facts about a case's input polynomial, computed with polyroot().
typedef struct {
  double dist;       // smallest ||z| - 1| over exact and computed roots
  int root_failure;  // polyroot() failed on a finite input
  int m;             // true degree
} case_info;

static void analyse(const arma_case *k, case_info *ci) {
  double zr[MAXN], zi[MAXN];
  int i;
  ci->dist = HUGE_VAL;
  ci->root_failure = 0;
  ci->m = 0;
  if (k->degenerate) return;
  if (k->known) {
    for (i = 0; i < k->m; ++i) ci->dist = fmin(ci->dist, fabs(hypot(k->rr[i], k->ri[i]) - 1.0));
  }
  for (i = k->n; i >= 1; --i) if (k->c[i] != 0.0) { ci->m = i; break; }
  if (ci->m == 0) return;
  if (polyroot((double *)k->c, ci->m, zr, zi) != 0) {
    ci->root_failure = 1;
    return;
  }
  for (i = 0; i < ci->m; ++i) {
    ci->dist = fmin(ci->dist, fabs(hypot(zr[i], zi[i]) - 1.0));
  }
}

// Number of roots of 1 + ma_1 z + ... inside |z| < 1 - slack.
static int roots_inside(int n, const double *ma, double slack) {
  double c[MAXN + 1], zr[MAXN], zi[MAXN];
  int i, m = 0, count = 0;
  c[0] = 1.0;
  for (i = 0; i < n; ++i) { c[i + 1] = ma[i]; if (ma[i] != 0.0) m = i + 1; }
  if (m == 0) return 0;
  if (polyroot(c, m, zr, zi) != 0) return -1;
  for (i = 0; i < m; ++i) count += hypot(zr[i], zi[i]) < 1.0 - slack;
  return count;
}

// Reflection computed from the known (exact) roots.
static void reflect_known(const arma_case *k, double *out) {
  arma_case t = *k;
  int i;
  for (i = 0; i < t.m; ++i) {
    const double d = t.rr[i] * t.rr[i] + t.ri[i] * t.ri[i];
    if (d < 1.0) { t.rr[i] /= d; t.ri[i] /= d; }  // 1 / conj(r) = r / |r|^2
  }
  expand_roots(&t);
  for (i = 0; i < k->n; ++i) out[i] = i < t.m ? t.c[i + 1] : 0.0;
}

// max_k |x_k - ref_k| / max(1, max_k |ref_k|)
static double coef_error(int n, const double *x, const double *ref) {
  double num = 0.0, den = 1.0;
  int i;
  for (i = 0; i < n; ++i) {
    num = fmax(num, fabs(x[i] - ref[i]));
    den = fmax(den, fabs(ref[i]));
  }
  return num / den;
}

// Conditioning of the problem itself: how far invertroot() moves when every
// input coefficient is perturbed by +-4 eps (relative), measured like
// coef_error() against the unperturbed output `base`. Two implementations
// that differ in rounding cannot be expected to agree more closely than this.
static double input_sensitivity(const arma_case *k, const double *base) {
  const uint64_t saved = g_state;
  double worst = 0.0, ma[MAXN];
  int t, i;
  g_state = 0xC0FFEEull;
  for (t = 0; t < 6; ++t) {
    for (i = 0; i < k->n; ++i)
      ma[i] = k->c[i + 1] * (1.0 + ((next_u64() >> 63) ? 4.0 : -4.0) * DBL_EPSILON);
    invertroot(k->n, ma);
    worst = fmax(worst, coef_error(k->n, ma, base));
  }
  g_state = saved;
  return worst;
}

static int g_fail = 0;
static FILE *g_out;  // stdout, or stderr while recording

static void gate(const char *what, double worst, double limit) {
  const int ok = worst <= limit;
  fprintf(g_out, "%-52s worst %.3g  gate %.3g  %s\n", what, worst, limit, ok ? "OK" : "FAIL");
  if (!ok) g_fail = 1;
}

typedef struct {
  int cases, band, degenerate, root_failure;
  int verdict_diff, verdict_diff_band, verdict_diff_degen;
  int status_diff, status_diff_band;
  int coef_cmp, coef_above_1e12;
  double worst_coef, worst_coef_ratio;
  int coef_band_cmp, coef_band_big, old_nonfinite;
  double worst_coef_band;
  int inside_new, inside_band_new, inside_band_old;
  int known_cmp, known_worse;
  int band_known, band_new_closer, band_old_closer;
  double worst_known_new, worst_known_old;
} stats;

// Compares one case: verdicts/statuses (old_v, old_s) and, if old_ma is not
// NULL, coefficients. Runs the independent checks on the new output.
static void check_case(int id, const arma_case *k, int old_v, int old_s,
                       const double *old_ma, stats *st, int verbose) {
  double ar[MAXN], ma[MAXN], ref[MAXN];
  case_info ci;
  int i, v, s, in_band, odd;

  analyse(k, &ci);
  in_band = ci.dist < k->band;
  odd = k->degenerate || ci.root_failure;
  for (i = 0; i < k->n; ++i) { ar[i] = -k->c[i + 1]; ma[i] = k->c[i + 1]; }
  v = archeck(k->n, ar);
  s = invertroot(k->n, ma);
  ++st->cases;
  st->band += in_band;
  st->degenerate += k->degenerate;
  st->root_failure += ci.root_failure;

  if (v != old_v) {
    if (odd) ++st->verdict_diff_degen;
    else if (in_band) ++st->verdict_diff_band;
    else {
      ++st->verdict_diff;
      if (verbose) fprintf(g_out, "  archeck differs: case %d (%s) old %d new %d dist %.3g\n",
                           id, k->family, old_v, v, ci.dist);
    }
  }
  if (s != old_s && !odd) {
    if (in_band) ++st->status_diff_band;
    else {
      ++st->status_diff;
      if (verbose) fprintf(g_out, "  invertroot status differs: case %d (%s) old %d new %d\n",
                           id, k->family, old_s, s);
    }
  }
  if (old_ma != NULL && !odd) {
    for (i = 0; i < k->n && isfinite(old_ma[i]); ++i) {}
    if (i < k->n) {
      ++st->old_nonfinite;
      old_ma = NULL;
      if (verbose) fprintf(g_out, "  upstream invertroot gave non-finite coefficients: case %d (%s), new status %d\n",
                           id, k->family, s);
    }
  }
  if (old_ma != NULL && !odd) {
    const double e = coef_error(k->n, ma, old_ma);
    if (in_band) {
      ++st->coef_band_cmp;
      st->worst_coef_band = fmax(st->worst_coef_band, e);
      if (e > 1e-12) ++st->coef_band_big;
    } else {
      ++st->coef_cmp;
      st->worst_coef = fmax(st->worst_coef, e);
      if (e > 1e-12) {
        // Above the nominal tolerance: must be within the conditioning.
        const double ratio = e / input_sensitivity(k, ma);
        ++st->coef_above_1e12;
        if (verbose && ratio > 1.0)
          fprintf(g_out, "  invertroot coefficients differ beyond conditioning: case %d (%s) rel %.3g, "
                  "%.3g x input sensitivity\n", id, k->family, e, ratio);
        st->worst_coef_ratio = fmax(st->worst_coef_ratio, ratio);
      }
    }
  }
  if (!odd) {
    // Independent: no root left inside the circle after invertroot. In the
    // noise band clustered roots are only resolved to the conditioning, so
    // there it is only counted (for both versions).
    if (in_band) {
      if (old_ma != NULL) {  // counted on the cases compared, for both versions
        st->inside_band_new += roots_inside(k->n, ma, 1e-9) > 0;
        st->inside_band_old += roots_inside(k->n, old_ma, 1e-9) > 0;
      }
      if (old_ma != NULL && k->known) {
        // Which version is closer to the reflection of the exact roots.
        double en, eo;
        reflect_known(k, ref);
        en = coef_error(k->n, ma, ref);
        eo = coef_error(k->n, old_ma, ref);
        ++st->band_known;
        st->band_new_closer += en * 2.0 + 1e-9 < eo;
        st->band_old_closer += eo * 2.0 + 1e-9 < en;
      }
    } else if (roots_inside(k->n, ma, 1e-9) > 0) {
      ++st->inside_new;
      if (verbose) fprintf(g_out, "  root left inside after invertroot: case %d (%s)\n", id, k->family);
    }
    // Against the reflection of the exact roots: the input coefficients are
    // rounded, so for ill-conditioned root sets neither version can match it
    // closely; the gate is that the new code is no worse than upstream.
    if (k->known && !in_band && old_ma != NULL) {
      double en, eo;
      reflect_known(k, ref);
      en = coef_error(k->n, ma, ref);
      eo = coef_error(k->n, old_ma, ref);
      ++st->known_cmp;
      st->worst_known_new = fmax(st->worst_known_new, en);
      st->worst_known_old = fmax(st->worst_known_old, eo);
      if (en > 2.0 * eo + 1e-12) {
        ++st->known_worse;
        if (verbose) fprintf(g_out, "  further from exact reflection than upstream: case %d (%s) new %.3g old %.3g\n",
                             id, k->family, en, eo);
      }
    }
  }
}

static void report(const stats *st, int full_coefs) {
  fprintf(g_out, "cases %d: in noise band %d, non-finite input %d, polyroot failure %d\n",
          st->cases, st->band, st->degenerate, st->root_failure);
  fprintf(g_out, "archeck verdict differences: %d outside band, %d in band (allowed), "
          "%d on non-finite input / polyroot failure (documented)\n",
          st->verdict_diff, st->verdict_diff_band, st->verdict_diff_degen);
  fprintf(g_out, "invertroot status differences: %d outside band, %d in band (allowed)\n",
          st->status_diff, st->status_diff_band);
  if (st->verdict_diff || st->status_diff) g_fail = 1;
  fprintf(g_out, "invertroot coefficients vs upstream, outside band: %d cases%s, "
          "%d differ by > 1e-12 rel. to max|c| (worst %.3g)\n",
          st->coef_cmp, full_coefs ? "" : " (recorded subset)", st->coef_above_1e12,
          st->worst_coef);
  gate("  those > 1e-12: difference / input sensitivity", st->worst_coef_ratio, 4.0);
  fprintf(g_out, "in band: %d compared, %d differ by > 1e-12 (worst %.3g)\n",
          st->coef_band_cmp, st->coef_band_big, st->worst_coef_band);
  fprintf(g_out, "in band, root inside |z| < 1 - 1e-9 after invertroot: new %d, upstream %d\n",
          st->inside_band_new, st->inside_band_old);
  fprintf(g_out, "in band, known roots: %d; clearly closer to exact reflection: new %d, upstream %d\n",
          st->band_known, st->band_new_closer, st->band_old_closer);
  fprintf(g_out, "upstream invertroot outputs with non-finite coefficients: %d (excluded)\n",
          st->old_nonfinite);
  fprintf(g_out, "outside band, root left inside after new invertroot: %d\n", st->inside_new);
  if (st->inside_new) g_fail = 1;
  fprintf(g_out, "known-root cases vs exact reflection: %d; worst new %.3g, upstream %.3g; "
          "new worse than 2 x upstream: %d\n",
          st->known_cmp, st->worst_known_new, st->worst_known_old, st->known_worse);
  if (st->known_worse) g_fail = 1;
}

static void degenerate_contract(void) {
  double x[3];
  int bad = 0;
  x[0] = NAN; bad |= archeck(1, x) != 0;
  x[0] = 0.5; x[1] = HUGE_VAL; bad |= archeck(2, x) != 0;
  x[0] = 0.0; x[1] = 0.0; bad |= archeck(2, x) != 1;
  bad |= archeck(0, x) != 1;
  x[0] = 1.0; bad |= archeck(1, x) != 1;   // root exactly on the circle
  x[0] = 2.0; x[1] = NAN; bad |= invertroot(2, x) != 0 || x[0] != 2.0;
  x[0] = 2.0; bad |= invertroot(0, x) != 0 || x[0] != 2.0;
  x[0] = 0.0; x[1] = 4.0; x[2] = 0.0;
  bad |= invertroot(3, x) != 1 || x[0] != 0.0 || fabs(x[1] - 0.25) > 1e-15 || x[2] != 0.0;
  x[0] = 0.0; x[1] = 0.0; x[2] = 1.0;    // 1 + z^3: all roots on the circle
  invertroot(3, x);
  bad |= fabs(x[0]) > 1e-12 || fabs(x[1]) > 1e-12 || fabs(x[2] - 1.0) > 1e-12;
  fprintf(g_out, "%-50s %s\n", "degenerate inputs (documented contract)", bad ? "FAIL" : "OK");
  if (bad) g_fail = 1;
}

#ifdef ARMA_ROOTS_RECORD

int main(void) {
  static arma_case k;
  stats st;
  char vbuf[101];
  int id, i, nv = 0;
  memset(&st, 0, sizeof st);
  g_out = stderr;
  printf("# arma_roots_reference.txt -- outputs of upstream ctsa archeck()/invertroot()\n"
         "# (tseries HEAD 7a5bf14, talg.c, linked with tseries polyroot) on the case set\n"
         "# of tool/arma_roots_accuracy.c. Lines: 'V <first id> <digits>' with digit =\n"
         "# 2 * archeck + invertroot status per case; 'C <id> <n> <coefficients>' with\n"
         "# the invertroot output of every %dth case.\n", COEF_EVERY);
  for (id = 0; id < NCASES; ++id) {
    double ar[MAXN], ma[MAXN];
    int v, s;
    make_case(id, &k);
    for (i = 0; i < k.n; ++i) { ar[i] = -k.c[i + 1]; ma[i] = k.c[i + 1]; }
    v = old_archeck(k.n, ar);
    s = old_invertroot(k.n, ma);
    check_case(id, &k, v, s, ma, &st, 1);  // live, full coefficients
    vbuf[nv++] = (char)('0' + 2 * v + s);
    if (nv == 100) { vbuf[nv] = 0; printf("V %d %s\n", id - 99, vbuf); nv = 0; }
    if (id % COEF_EVERY == 0) {
      printf("C %d %d", id, k.n);
      for (i = 0; i < k.n; ++i) printf(" %.17g", ma[i]);
      printf("\n");
    }
  }
  report(&st, 1);
  degenerate_contract();
  fprintf(g_out, "%s\n", g_fail ? "SOME GATES FAILED" : "all gates OK");
  return g_fail;
}

#else

int main(int argc, char **argv) {
  static arma_case k;
  static int verdict[NCASES], status[NCASES], has_coef[NCASES];
  static double coefs[NCASES / COEF_EVERY + 1][MAXN];
  static char line[1 << 16];
  stats st;
  FILE *f;
  int id, nrec = 0;

  if (argc < 2) {
    fprintf(stderr, "usage: %s arma_roots_reference.txt\n", argv[0]);
    return 2;
  }
  f = fopen(argv[1], "r");
  if (f == NULL) { perror(argv[1]); return 2; }
  for (id = 0; id < NCASES; ++id) verdict[id] = -1;
  while (fgets(line, sizeof line, f) != NULL) {
    int first, n, consumed, i;
    char digits[128];
    if (line[0] == '#' || line[0] == '\n') continue;
    if (line[0] == 'V' && sscanf(line + 1, "%d %127s", &first, digits) == 2) {
      for (i = 0; digits[i] && first + i < NCASES; ++i) {
        verdict[first + i] = (digits[i] - '0') >> 1;
        status[first + i] = (digits[i] - '0') & 1;
      }
    } else if (line[0] == 'C' && sscanf(line + 1, "%d %d%n", &first, &n, &consumed) == 2 &&
               first >= 0 && first < NCASES && first % COEF_EVERY == 0 && n <= MAXN) {
      char *cur = line + 1 + consumed;
      for (i = 0; i < n; ++i) coefs[first / COEF_EVERY][i] = strtod(cur, &cur);
      has_coef[first] = 1;
      ++nrec;
    } else {
      fprintf(stderr, "bad line: %.60s\n", line);
      return 2;
    }
  }
  fclose(f);

  memset(&st, 0, sizeof st);
  g_out = stdout;
  for (id = 0; id < NCASES; ++id) {
    if (verdict[id] < 0) { fprintf(stderr, "case %d missing from fixture\n", id); return 2; }
    make_case(id, &k);
    check_case(id, &k, verdict[id], status[id],
               has_coef[id] ? coefs[id / COEF_EVERY] : NULL, &st, 1);
  }
  printf("recorded coefficient vectors: %d\n", nrec);
  report(&st, 0);
  degenerate_contract();
  printf("%s\n", g_fail ? "SOME GATES FAILED" : "all gates OK");
  return g_fail;
}

#endif
