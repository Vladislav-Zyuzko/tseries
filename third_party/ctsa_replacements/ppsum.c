// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// Replacement for ctsa's `ppsum()` (upstream `src/talg.c`), which carried
// "Copyright (C) 1997-2000 Adrian Trapletti" (R package tseries, GPL) and has
// been removed from the vendored tree without being read. See
// third_party/ctsa/PROVENANCE.md, "Clean-room replacements". The contract was
// taken from the prototype in talg.h and from black-box tests of the
// upstream binary.
//
// Contract
//   *sum += (2 / n) * sum_{s=1}^{l} w(s, l) * sum_{t=s}^{n-1} u[t] u[t-s],
//   w(s, l) = 1 - s / (l + 1)   (Bartlett weights).
//   It ADDS to *sum: callers put the lag-0 term (1/n) sum u[t]^2 there first.
//   l <= 0 adds nothing; lags s >= n contribute nothing.
//
// This is the autocovariance part of the Newey-West (Bartlett-kernel)
// long-run variance estimator,
//   s^2(l) = (1/n) sum_t u_t^2 + (2/n) sum_{s=1}^{l} w(s, l) sum_{t>s} u_t u_{t-s},
// W. K. Newey and K. D. West, "A simple, positive semi-definite,
// heteroskedasticity and autocorrelation consistent covariance matrix",
// Econometrica 55 (1987) 703-708; used in the same form by the KPSS test,
// D. Kwiatkowski, P. C. B. Phillips, P. Schmidt and Y. Shin, J. Econometrics
// 54 (1992) 159-178, eq. (10), and the Phillips-Perron test.

#include "talg.h"

void ppsum(double* u, int n, int l, double* sum) {
  double total = 0.0;
  int s, t;
  if (u == NULL || sum == NULL || n <= 0 || l <= 0) return;
  for (s = 1; s <= l && s < n; ++s) {
    const double w = 1.0 - (double)s / ((double)l + 1.0);
    double cov = 0.0;
    for (t = s; t < n; ++t) cov += u[t] * u[t - s];
    total += w * cov;
  }
  *sum += 2.0 * total / (double)n;
}
