<!--
Draft of an issue for https://github.com/rafat/ctsa. NOT SENT — whether and
when to send it is the project owner's decision. Written 2026-10-05 against
the vendored pin (see third_party/ctsa/PROVENANCE.md, "Pin").
-->

# Exact likelihood is wrong for pure AR(1) models: `iupd` stays 0 in the AR(1) special case of `fas154*`

## Summary

For any model whose AR/MA polynomials multiply out to `ip == 1, iq == 0` —
ARIMA(1,d,0) with no seasonal AR/MA terms, for any `d` and any seasonal
differencing `D` — the exact (AS 154) likelihood treats the first
observation as having variance σ² instead of the stationary σ²/(1 − φ²), and
drops the ½·ln(1 − φ²) term. Since this is the objective being maximised, the
estimates of φ, the mean/regression coefficients and σ², the reported
log-likelihood/AIC and (through the estimates) the forecasts are not the exact
MLE ones.

Affected: `arima_exec`, `sarima_exec`, `sarimax_exec` (and `ar_exec`) with
the exact-likelihood methods (MLE, CSS-MLE).

## Cause

`emle.c`, in `fas154()`, `fas154_seas()` and `fas154x_seas()`:

```c
iupd = 0;

if ((ip == 1 && iq == 0) || (ip == 0 && iq == 0)) {
    *V = 1.0;
    *A = 0.0;
    *P = 1.0 / (1.0 - phi[0] * phi[0]);
}
else {
    iupd = 1;
    ifault = starma(ip, iq, phi, theta, A, P, V);
}
```

The special case follows AS 154's own note (`starma` is not suitable for
AR(1); set `V(1) = 1`, `A(1) = 0`, `P(1) = 1/(1 − φ²)`), and the resulting
`P` is — like `starma`'s output — already the *predicted* variance for the
first observation. But `iupd` stays 0, so `karma()` performs a prediction step
before the first observation. With `ir == 1` that step recomputes `P` from
`V`:

```c
if (iupd != 1 || i > 0) {
    ...
    P[ind] = V[ind];          /* P = V = 1: the stationary P0 is lost */
    ...
}
```

so the first observation is filtered with variance 1 (times σ²) instead of
1/(1 − φ²). `forkal()` handles the same case correctly (it sets `iupd = 1`
after initialising `P = 1/(1 − φ²)` for `ir == 1`), which is why forecasts
*at given parameters* are right.

## Proposed fix

Set `iupd = 1` in the special case of all three functions (equivalently,
initialise `iupd = 1` for both branches):

```diff
     if ((ip == 1 && iq == 0) || (ip == 0 && iq == 0)) {
         *V = 1.0;
         *A = 0.0;
         *P = 1.0 / (1.0 - phi[0] * phi[0]);
+        iupd = 1;
     }
```

For `ip == 0 && iq == 0` (white noise after differencing) `phi[0]` is 0, so
`P0 = V = 1` and the result is identical either way.

## Minimal reproduction

AR(1) with a mean on the first 40 values of R's `LakeHuron`, fitted by
CSS-MLE, compared with the exact AR(1) likelihood evaluated in closed form at
ctsa's *own* estimates (so the comparison does not depend on the optimiser).
The closed form uses ctsa's `log(2 * 3.14159)` so the constants match.

```c
#include <math.h>
#include <stdio.h>
#include "ctsa.h"

static double y[] = {580.38, 581.86, 580.97, 580.80, 579.79, 580.39, 580.42,
  580.82, 581.40, 581.32, 581.44, 581.68, 581.17, 580.53, 580.01, 579.91,
  579.14, 579.16, 579.55, 579.67, 578.44, 578.24, 579.10, 579.09, 579.35,
  578.82, 579.32, 579.01, 579.00, 579.80, 579.83, 579.72, 579.89, 580.01,
  579.37, 578.69, 578.19, 578.67, 579.55, 578.92};

int main(void) {
  const int n = sizeof(y) / sizeof(y[0]);
  sarima_object obj = sarima_init(1, 0, 0, 0, 0, 0, 0, n); /* AR(1) + mean */
  sarima_setMethod(obj, 0);                                 /* CSS-MLE */
  sarima_exec(obj, y);
  double phi = obj->phi[0], mu = obj->mean;

  /* exact AR(1) log-likelihood at (phi, mu), sigma^2 profiled out */
  double s = (1 - phi * phi) * (y[0] - mu) * (y[0] - mu);
  for (int t = 1; t < n; ++t) {
    double e = (y[t] - mu) - phi * (y[t - 1] - mu);
    s += e * e;
  }
  double exact = -0.5 * n * (log(2 * 3.14159) + log(s / n) + 1)
                 + 0.5 * log(1 - phi * phi);
  printf("phi = %.6f  mu = %.6f  retval = %d\n", phi, mu, obj->retval);
  printf("ctsa loglik           = %.6f\n", obj->loglik);
  printf("exact AR(1) loglik    = %.6f   (at the same phi, mu)\n", exact);
  printf("ctsa sigma2           = %.6f\n", obj->var);
  printf("exact sigma2 (S/n)    = %.6f\n", s / n);
  sarima_free(obj);
  return 0;
}
```

Output (clang 22, Windows x86_64, `-O2`; other compilers may differ in the
last digits):

| | current | with the fix | statsmodels (exact MLE) |
| --- | --- | --- | --- |
| φ | 0.840884 | 0.817150 | 0.817150 |
| μ | 580.004244 | 579.801616 | 579.801620 |
| ctsa log-likelihood | −32.815643 | −33.132663 | −33.132663 ¹ |
| exact log-likelihood at ctsa's (φ, μ) | −33.263641 | −33.132663 | |
| σ² | 0.302071 | 0.298559 | 0.298559 |

¹ statsmodels' −33.132679 plus the constant ½·n·ln(π/3.14159) that ctsa's
`log(2 * 3.14159)` adds (see the separate note on that constant).

## How the fix was checked

In the copy vendored by the Dart package `tseries` (Sweet Limit project),
2026-10-05:

1. **The likelihood function at fixed parameters.** `fas154`, `fas154_seas`
   and `fas154x_seas` called directly, at statsmodels' optimum and at ctsa's
   optimum for 34 models — AR(1) with and without a mean, with a regressor,
   (1,1,0) with and without drift, (1,0,0)(0,1,0)[12], (1,1,0)(0,1,0)[12],
   (0,0,0)(1,0,0)[1], synthetic AR(1) with φ ∈ {−0.8, 0.05, 0.98, 0.995} and
   n ∈ {30, 200}, plus non-AR(1) neighbours (AR(2), MA(1), (0,1,1),
   (1,0,0)(1,0,0)[12], (1,1,0)(0,1,1)[12]) — against statsmodels 0.15.0
   `SARIMAX(..., concentrate_scale=True).loglike` (stationary initialisation)
   on the same differenced series. Before: off by 2e-7 … 134 in exactly the
   AR(1) cases. After: |Δ| ≤ 6e-12 relative everywhere.
2. **The fitted optimum** of every case through `sarima_exec`/`sarimax_exec`:
   within 1e-4 of statsmodels in log-likelihood, 1.6e-5 in the coefficients,
   1 ± 2.3e-5 in σ².
3. **Forecasts and their standard errors** at the fitted parameters against
   statsmodels `get_forecast`: ≤ 1e-9 relative.
4. **Nothing else moves:** 6 787 fits (11 series, 26 orders, all methods,
   with/without mean, drift, regressor) dumped bit for bit before and after;
   only fits that reach this branch and run the MLE step changed.

## Related observations (separate issues, for reference)

- `ctsa.c` computes the log-likelihood with `log(2 * 3.14159)`, so every
  reported log-likelihood is ½·n·ln(π/3.14159) ≈ 4.2e-7·n too large.
- The log-likelihood and σ² that `as154*` report come from the *last*
  objective evaluation, which happens inside `hessian_fd` at a perturbed
  point, not at the returned estimate (up to ~3e-5 in log-likelihood in our
  tests; which point is last depends on the compiler's evaluation order of
  the four calls in one expression).
- `karma()`'s `*nit != 0` branch (unreachable today: every caller passes
  `nit = 0`) has `for (j = 0; j < iq; ++i)` — `i` instead of `j`.
