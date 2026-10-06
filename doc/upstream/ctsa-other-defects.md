<!--
Drafts of short issues for https://github.com/rafat/ctsa. NOT SENT — whether
and when to send them is the project owner's decision. Written 2026-10-05
against the vendored pin (third_party/ctsa/PROVENANCE.md, "Pin"). Line numbers
are upstream's at that pin, approximately.
-->

# Other ctsa defects found while vendoring (drafts)

One section per issue. The AR(1) likelihood defect has its own draft:
`ctsa-ar1-likelihood.md`.

## 1. `as154x()` / `cssx()`: `rank()` is called with rows and columns swapped

`XX` (the differenced regressors) is stored column-major — column `j` at
`XX + N*j`, as `fas154x_seas()` reads it (`reg[j*N + i]`) and `regress()`
consumes it — but the collinearity check calls `rank(XX, N, ncxreg)`, and
`rank(A, M, N)` reads `A` as a row-major M×N matrix. The rank is therefore
that of a reshaped matrix:

- false negatives: an all-zero regressor, or a duplicated column, passes the
  check; the fit then ends in NaN (retval 15) or returns retval 1 with
  meaningless coefficients;
- false positives: two full-rank regressors whose values are held for two
  steps (e.g. sampled at half the series' rate) are rejected as collinear
  (retval 7) at d = 0 — every reshaped row has the form (a, a).

**Fix:** `rank(XX, ncxreg, N)` in both functions. A column-major N×k buffer is
the row-major k×N transpose, which has the same rank; `rank()` transposes
internally when M < N. Full-rank designs are accepted exactly as before. The
same pattern (`rank(xreg, N, r)`) is in the auto-ARIMA path of `ctsa.c`.

Also in `matrix.c`: `rank_c()` returns early without freeing `U`, `V`, `q`
when the SVD does not converge (seen only for designs spanning more than ~200
orders of magnitude).

## 2. `boxcox()` reads an uninitialised variable

```c
double boxcox(double *x, int N, double *lambda, double *y) {
    ...
    int cc, cn;
    ...
    if (cn == 1) {
        printf("Input vector cannot have negative values. Exiting. \n");
        exit(-1);
    }
```

`cn` is never assigned (the negative-value check was apparently never
written), so the behaviour is undefined; in optimised builds we observed the
process exiting on every call, including on positive data. **Fix:** compute
it, e.g. `cn = 0; for (i = 0; i < N; ++i) if (x[i] <= 0) cn = 1;` (Box-Cox
needs strictly positive data), and preferably return an error code instead of
calling `exit()` from a library.

Minor, same file: `inv_boxcox_eval()` is declared `double` but has no
`return` on its exit path.

## 3. `ar_estimate()`: AIC order selection never compares against order 1

```c
lvar = log(var);
aic = lvar + 2 * (double)(i + 1) / N;
if (i == 0) {
    aic = aic0;      /* overwrites the order-1 AIC with the initial 0.0 */
    p = 1;
}
else {
    if (aic < aic0) {
        aic0 = aic;
        p = i + 1;
    }
}
```

At `i == 0` the assignment goes the wrong way, so `aic0` stays at its initial
`0.0` and every later order is compared with 0 instead of with the order-1 AIC
(the selected order depends on the sign of `log(var)`, i.e. on the units of
the series). **Fix:** `aic0 = aic;` in the `i == 0` branch.

## 4. The Gaussian constant uses `3.14159`

`ctsa.c` (`sarima_exec`, `sarimax_exec`, `arima_exec`, …) computes

```c
obj->loglik = -0.5 * (obj->Nused * (2 * obj->loglik + 1.0 + log(2 * 3.14159)));
```

so every reported log-likelihood (and AIC) is larger than the exact one by
½·n·ln(π/3.14159) ≈ 4.2e-7·n (n = observations after differencing). Harmless
for comparisons on the same series, but it makes the numbers disagree with
R/statsmodels beyond rounding. There are eight such lines in `ctsa.c`
(`arima_exec` ×2, `sarimax_exec` ×3, `sarima_exec` ×2, `ar_exec`); every other
π in the library is already full precision. **Fix:** use a full-precision
constant (as `PIVAL` in `erfunc.h`; `M_PI` needs `_USE_MATH_DEFINES` on MSVC).
The conversion runs after the optimiser, so estimates do not change — we
checked bit for bit: only the log-likelihood (by exactly
−½·n·ln(π/3.14159)) and AIC move.

## 5. `forkal()` leaves forecasts and their MSEs unset for white noise

For `ip == 0 && iq == 0` (e.g. ARIMA(0,1,0), or (0,0,0) with a mean),
`forkal()` sets `ifault = 4` and returns before writing `Y`/`AMSE`; the
`*_predict()` wrappers ignore the return value and add the mean to whatever
`xpred` held, and `amse` is left as the caller's (typically uninitialised)
buffer. The model is perfectly well defined — the h-step MSE of a random walk
is h·σ², of white noise σ² — so the case needs handling (or at least an error
code propagated to the caller).

**Reproduce:** `sarima_init(0, 1, 0, 0, 0, 0, 0, N)`, `sarima_exec`, then
`sarima_predict(obj, x, 6, xpred, amse)` with `xpred` zero-filled and `amse`
fresh from `malloc`: `xpred` stays 0 (not the last observation) and `amse` is
whatever `malloc` returned — under MemorySanitizer every element is reported
uninitialised (`ctsa.c`, `xpred[i] += wmean`, and the caller's first read of
`amse`).

**Fix (what we do on our side):** a white-noise model is the AR(1) model with
φ = 0, which AS 182 already handles (`ir == 1`, `P = 1/(1 − φ²) = 1`, no
`starma`). Calling `forkal(1, iq, id, phi /* phi[0] == 0 */, ...)` when
`ip == iq == 0` gives the exact forecasts and MSEs; alternatively drop the
`ip² + iq² == 0` check in `forkal()` and size `ir` as `max(ip, iq + 1, 1)`.

## 5a. `sarima_predict()` / `sarimax_predict()` read past `obj->res`

Both copy the residuals with

```c
for (i = 0; i < N; ++i) {
    ...
    resid[i] = obj->res[i];
}
```

but `obj->res` has `N - d - s*D` elements (`sarima_init`/`sarimax_init`
allocate `Nused` of them, followed by the `vcov` block of (k + M)² doubles).
Whenever `d + s*D > (k + M)²` the loop reads past the end of the object's heap
block — AddressSanitizer: *heap-buffer-overflow READ, ctsa.c:901
(sarima_predict) / ctsa.c:1008 (sarimax_predict)*. That covers every
differenced white-noise model (k = M = 0) and e.g. (0,0,1)(0,1,0)[12]. The
copied values are dead (`karma()` overwrites `resid[0 .. N - id)` before
`forkal()` reads any of it), so the result is unaffected; the read itself can
fault. **Fix:** copy only `Nused` values (or none).

## 5b. Heap overflow in `fcssx()` for short seasonal series

`xlik_css_init()` allocates `2r + (3 + M)·N` doubles (N = differenced length,
M = regressors incl. intercept); `fcssx()` then zeroes
`obj->x[offset + N + i]` for `i < ncond = p + s*P`. When `p + s*P > (2 + M)·N`
this writes past the block — AddressSanitizer: *heap-buffer-overflow WRITE,
emle.c:2288*, e.g. `sarimax_init(0,0,0, 1,1,0, 12, 0, 0, 15)` with method 0 or
2. Before that point (`N ≤ p + s*P`) the conditional sum of squares is empty
and the optimiser starts from NaN. **Fix:** reject `N <= p + s*P` in
`cssx()` (there is nothing to condition on). `fcss_seas()` has the same
indexing but allocates with the undifferenced length and did not overflow in
our runs.

## 5b'. Heap overflow differencing the regressors under d > 0 and D > 0 (2026-10-06)

In `cssx()` and `as154x()` (`emle.c`), `XX` is allocated with
`N*(r+1)` doubles, N the fully differenced length (`length - d - s*D`). With
regressors and `D > 0`, each column is first seasonally differenced into
`XX + N1*i` with `N1 = length - s*D` (the `diffs(x0+length*i, …, XX+N1*i)`
loop), i.e. up to `N1*ncxreg` doubles. That exceeds the block exactly when
`d*ncxreg > N` — AddressSanitizer: *heap-buffer-overflow WRITE in `diffs()`
(talg.c) from `cssx()`*, e.g. `sarimax_init(0,2,0, 0,1,0, 4, 3, 0, 11)` (three
regressors, 11 points) with any method. **Fix:** allocate `XX` with
`length*(r+1)` doubles (as `x0` already is), or seasonally difference into a
separate scratch buffer.

## 5c. Box-Jenkins (`nlalsms`) memory errors

- p = q = P = Q = 0 with d + D > 0 (nothing to estimate, `k = 0`): the 0×0
  system handed to `linsolve()` (matrix.c:1302) is written out of bounds —
  heap corruption, then a segfault (it does not reach any `exit()`).
- `avaluem()` reads `w[t - i]` for `t < p + s*P` and `e[t2 - i]` for
  `t2 < q + s*Q` (boxjenkins.c:410/421): out of bounds when the expanded lag
  exceeds the differenced length (seasonal model on less than a season).
- `USPE_seasonal()` with s = 0 (every non-seasonal Box-Jenkins fit) allocates
  `cv` with `K = (p + q + 1)·s = 0` elements and reads `cv[0]`
  (boxjenkins.c:208).
- `USPE_seasonal()` on seasonal fits of short series uses uninitialised
  autocovariances (MemorySanitizer, boxjenkins.c:240/290, then
  `linsolve`/`pludecomp`); the results of such fits vary run to run.
  Mechanism: it asks `autocovar()` for K = (P + Q + 1)·s lags and reads
  lags 0, s, …, (P + Q)·s; when K > N, `autocovar()` sets `M = N - 1` and
  fills only lags 0 … N − 2 (lag N − 1, which it could compute, is skipped
  too). So lag (P + Q)·s is never written when N ≤ (P + Q)·s + 1 — exactly
  where we see reports (sweep over s ∈ {4, 7, 12}, P, Q ≤ 1). **Fix:** in
  `autocovar()` clamp to `M = N` (and zero the requested lags it cannot
  compute), or have `nlalsms()` reject N ≤ (P + Q)·s + 1.

## 5d. `bfgs_min()` reads its gradient uninitialised when f(x0) is NaN

`secant.c`: when the first objective value is NaN (`rcode = 15`), the code
still computes the stopping test from `jac[]`, which was never filled
(MemorySanitizer, secant.c:514/519; buffer from secant.c:432). Reached from
`css_seas()` when the differenced series is not longer than `p + s*P` (empty
CSS sum); the subsequent MLE starts from garbage and the fit's outcome
(status and estimates) differs between runs. **Fix:** return 15 before using
`jac[]` when f(x0) is not finite, and/or reject N ≤ p + s·P in `css_seas()`.
(The MA analogue, N ≤ q + s·Q with N > p + s·P, is fine: the CSS sum is not
empty.)

## 6. Reported log-likelihood / σ² are from a perturbed point (compiler-dependent)

`as154()`, `as154_seas()`, `as154x()`, `css()`, `css_seas()` and `cssx()`
copy `obj->ssq` (→ `*var`), `obj->loglik` and the residuals out of the
likelihood object *after* `hessian_fd()`. The objective functions store those
fields as a side effect, so the values belong to the Hessian's last
evaluation, not to the returned estimate `tf`. Which evaluation is last is
unspecified C: `hessian_fd()` computes
`((FUNCPT_EVAL(f, xij) - FUNCPT_EVAL(f, xi)) - (FUNCPT_EVAL(f, xj) - FUNCPT_EVAL(f, x)))`.
GCC and clang happen to evaluate `f(x)` (= `tf`) last, so their numbers are
right; MSVC evaluates `f(xi)` last, i.e. `tf` with its last parameter moved by
one finite-difference step (≈ 6e-6 relative). On MSVC builds we measured up to
5.6e-5 in log-likelihood and 2.3e-5 relative in σ² on ordinary fits (Nile,
lynx, LakeHuron, AirPassengers, USAccDeaths, WWWusage), 1.2e-4 relative in σ²
near an MA unit root, and much more on degenerate fits. AIC follows the
log-likelihood. Re-evaluating the objective at exactly that stencil point
reproduces the reported numbers bit for bit, which pins the mechanism.

**Fix** (six one-line additions, one per function, right before `*var = …`):

```c
FUNCPT_EVAL(&as154_min, tf, pq);   /* as154, as154x, as154_seas */
FUNCPT_EVAL(&fcss_min, tf, pq);    /* css */
FUNCPT_EVAL(&css_min, tf, pq);     /* css_seas, cssx */
```

Coefficients, the Hessian/covariance and forecasts are untouched (the Hessian
is computed before; `forkal()` estimates its own σ²). Afterwards σ² and the
log-likelihood equal the objective evaluated at the estimate, and are the same
on every compiler for the same estimate.

## 7. `karma()`: infinite loop in the (currently unreachable) fast-recursion branch

In the `*nit != 0` branch: `for (j = 0; j < iq; ++i)` increments `i`, not
`j`. Every current caller passes `nit = 0`, so it is unreachable today.
