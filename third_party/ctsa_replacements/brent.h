// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// tseries replacement for ctsa's `brent.h` (see brent_local_min.c and
// third_party/ctsa/PROVENANCE.md, "Clean-room replacements").
//
// ctsa's `neldermead.h` does `#include "brent.h"`; the upstream header was
// removed from the vendored tree, so the include resolves here (this directory
// is on the build's include path). It provides what the vendored code needs
// from that include: the `brent_local_min` prototype and, transitively,
// `secant.h` (which pulls in the optimiser typedefs `neldermead.h` uses).
//
// `brent_zero` (upstream's root finder) is intentionally not declared: nothing
// in ctsa or in the tseries shim calls it, so it was dropped, not replaced.

#ifndef TSERIES_CTSA_REPLACEMENT_BRENT_H_
#define TSERIES_CTSA_REPLACEMENT_BRENT_H_

#include "secant.h"

#ifdef __cplusplus
extern "C" {
#endif

// Minimises `funcuni` on [a, b] (Brent 1973, ch. 5). On return `*x` is the
// abscissa of the minimum; the function value there is returned.
//
// `t` is the absolute and `eps` the relative tolerance: the interval is shrunk
// until the minimum is located to within 2 * (eps * |x| + t).
double brent_local_min(custom_funcuni *funcuni, double a, double b, double t,
                       double eps, double *x);

#ifdef __cplusplus
}
#endif

#endif  // TSERIES_CTSA_REPLACEMENT_BRENT_H_
