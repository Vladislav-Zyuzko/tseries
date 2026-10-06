# Vendored `ctsa` — provenance & license notes

This directory contains a **verbatim, pinned** vendoring of the C sources of
[`ctsa`](https://github.com/rafat/ctsa) by Rafat Hussain, a univariate
time-series / ARIMA-SARIMA library in ANSI C.

## Pin

| | |
| --- | --- |
| Upstream | https://github.com/rafat/ctsa |
| Commit | `6e39431fe0832bfecd0ce2e99b044599b0a74418` |
| Commit date | 2026-06-13 |
| Commit subject | "NaN bug solved for ARIMA 0,1,0 model" |
| Vendored on | 2026-07-12 |

This exact revision was chosen deliberately: it is **not older than** the
2026-06-13 fix for the NaN bug in the ARIMA(0,1,0) model, which our
due-diligence flagged as a must-have. Do not downgrade below this commit.

## What was vendored

- `src/**.c`, `src/**.h` — the full C source tree (compiled by `hook/build.dart`),
  **minus** the code removed for licensing reasons on 2026-10-05: the LGPL
  units (`brent.c`, `brent.h`, `erfunc.c`, two sections of `dist.c`), the
  unattributed MIT routines of `lls.c`, `polyroot.c` (ACM Algorithm 419,
  non-commercial licence), `ppsum()` from `talg.c` (from the GPL R package
  tseries), and unreachable code of restrictive or GPL origin (Numerical
  Recipes routines in `matrix.c`, AS 197 in `emle.c`/`talg.c`, R-derived
  parameter transforms in `talg.c`), and — in 0.8.0 — the whole auto-ARIMA
  cluster, a C port of GPL R code (`autoutils.c`, `unitroot.c`, `seastest.c`,
  `stl.c` and their headers deleted; the auto functions of `ctsa.c`/`ctsa.h`
  and `supsmu` of `talg.c`/`talg.h` cut), and — after 0.8.0 was tagged —
  `initest.c`/`initest.h` (Burg translation of unknown licence) with its
  unreachable callers in `ctsa.c` — see "Local modifications",
  "Clean-room replacements" and "Auto-ARIMA removal" below. Dead code is dropped at link time
  (`--gc-sections` / `-dead_strip`, hidden visibility), but note that the
  **source** of everything in this directory is published with the package,
  reachable or not. The tree is compiled **without floating-point
  contraction** (`-ffp-contract=off` for clang/GCC, `/fp:precise` for MSVC;
  2026-10-06) so that no `a * b + c` becomes a fused multiply-add on arm64 —
  a build flag only, no source change (`hook/build.dart`,
  `doc/platform_verification.md`).
- `COPYRIGHT` — upstream BSD-3-Clause copyright (Rafat Hussain, 2014).
- `components_copyright/` — upstream's per-component license texts
  (`LGPL_COPYRIGHT` removed on 2026-10-05 together with the last LGPL code).
- `README.upstream.md` — upstream README (renamed to avoid clashing with ours).

The header ffigen parses is **`src/ctsa.h`**, NOT the repo's `header/ctsa.h`.
Reason: the two are out of sync — `src/ctsa.h` (the one the `.c` files actually
`#include`) declares extra struct fields (e.g. `cssml`) and extra structs that
`header/ctsa.h` lacks. Binding against the stale public header would give wrong
struct field offsets and corrupt memory. `header/ctsa.h` is intentionally NOT
vendored to remove the trap.

## License status

**LGPL: resolved (2026-10-05).** Upstream's SPDX audit (contributed to ctsa in
PR #15 by bavay, merged 2024-10-25, which tags files and individual functions)
marks `src/brent.c`, `src/brent.h`, `src/erfunc.c` and two sections of
`src/dist.c` (`r8_max`, `betainv`) as LGPL-3.0-or-later. All of them are gone
from this tree and replaced by independent BSD-3-Clause implementations; see
"Clean-room replacements". The rest is BSD-3-Clause, MIT (`log1p` in `dist.c`,
GeographicLib, attributed in `THIRD_PARTY_NOTICES`) and CC0/public domain, plus the open
questions listed under "Open licensing questions" below.

History, kept for context: on 2026-07-13 the project owner knowingly deferred
the LGPL question ("do not delete, do not modify, do not rewrite until before a
public release"). On 2026-10-05, ahead of a pub.dev release under BSD-3-Clause,
the owner started the replacement work recorded here; the deferral no longer
applies.

## Local modifications to the vendored tree

The rule is "do not modify the vendored code unless it cannot be avoided". We
have departed from it repeatedly (entries 1-11, 13, 14 and 17; entries 12, 15 and 16 record defects
worked around in the shim without touching vendored files; whole files deleted for
licensing reasons are listed under "Clean-room replacements"). Each departure is listed here so it is **explicit and traceable** rather
than silent, because every local patch is a re-apply-and-re-verify burden on the
next ctsa version bump.

### When patching upstream is allowed

**Only when the problem cannot be solved on our side of the boundary.** We own
the C shim (`native/tseries_ctsa.*`) and the Dart wrapper; anything expressible
there belongs there, and anything avoidable by excluding a source file from
`hook/build.dart` (as was done for the wavelet cluster) belongs there instead.
Patching upstream blurs the "whose bug is this?" line and costs us on every
update — it needs to earn its place.

A worked example of the rule biting: ctsa's `sarima_init` leaves `loglik`/`aic`
uninitialised, so they read back as malloc garbage for any method that never
assigns them. Initialising them to NaN in `*_init` was considered and
**rejected** — the shim is the only reader of those fields, so simply not
reading them unless ctsa assigned them achieves the same result with no
vendored diff. See the gating in `tseries_sarima_fit`.

### 1. `src/boxcox.c` — `inv_boxcox_eval` missing `return`

**Superseded 2026-10-05:** `boxcox.c` was deleted with its last caller (entry
11), so there is nothing left to patch. Kept for the record.

- **What:** added `return lambda;` at the end of `inv_boxcox_eval`.
- **Why:** the function is declared `double` but had no `return` on its exit
  path (`-Wreturn-type` / MSVC C4715) — undefined behaviour if its value were
  ever consumed. It is never called today, but its translation unit is
  compiled (`ctsa.c` `sarimax_wrapper` uses `boxcox_eval`; until 0.8.0 the
  removed `seastest.c` also used `boxcox`), so the warning was live in every
  build.
- **Not fixable on our side:** it is a defect inside a function body we compile;
  no shim or build-hook change can add a missing `return`.
- **Signature/ABI:** unchanged. The value is never consumed.
- **On version bump:** check whether upstream added a `return`. If so, drop this
  patch. If not, re-apply and re-check the warning is gone.

### 2. `src/emle.c` — memory leaks on the `retval` 10/12 early return

Two sites, same defect: an early `return ERR` that skips a live `malloc`.

- **`as154()`** (~line 1476): leaked `x`, `(N-d)` doubles, allocated at the top
  of the function. Added `free(x)` before the return.
- **`as154_seas()`** (~line 2713): leaked `inp2`, `(N - s*D)` doubles. Added
  `free(inp2)` before the return. `x` is not yet allocated at that point, so
  `inp2` is the only live allocation.
- **Why it matters:** `as154_seas` is the path `sarima_exec` actually calls, and
  the app re-fits on chart interaction, so this leaked continuously in a live
  process rather than in a short-lived isolate. Codes 10/12 are reachable in
  practice (see `test/arima_fit_status_test.dart`, which reproduces both).
- **Not fixable on our side:** `x`/`inp2` are function-local allocations that
  never escape; nothing outside `emle.c` can free them.
- **Control flow / return values:** unchanged — only the `free` was added.
- **Note:** `as154x()` (~line 1649) has the same early return but is **not**
  leaking: both of its early returns precede its first `malloc`. It was checked
  and deliberately left alone.
- **On version bump:** re-check all three `checkroots_cerr` call sites in
  `emle.c` for early returns that skip a live allocation.

### 3. `src/emle.c` — collinearity check reads a column-major matrix as row-major

- **What:** in `as154x()` and `cssx()` (the two SARIMAX estimators), the call
  `rank(XX,N,ncxreg)` became `rank(XX,ncxreg,N)` — the two dimension arguments
  swapped. Nothing else.
- **The defect:** `XX` (the differenced regressors) is **column-major** —
  column `j` at `XX + N*j`, exactly as `fas154x_seas` reads it
  (`reg[j*N+i]`) and `regress()` consumes it. `rank(A,M,N)` (matrix.c) reads
  `A` as **row-major** M×N. Upstream therefore computed the rank of a
  *reshaped* matrix, with two consequences, both reproduced on a host bench
  (2026-10-03/04):
  - **false negatives:** an all-zero regressor passed the check and the fit
    ran into NaN (`retval` 15); a duplicated column passed with `retval` 1 and
    garbage coefficients;
  - **false positives:** two full-rank regressors whose values are held for two
    steps (sampled at half the series' rate — routine for resampled exogenous
    data) were rejected as collinear (`retval` 7) at `d = 0`. The
    reshaped rows are all of the form `(a, a)`, i.e. rank 1.
- **Why it could not be fixed on our side:** the false negatives *can* be —
  and are — handled by the shim: `tseries_sarimax_fit` screens every regressor
  itself before ctsa sees it (zero / constant-after-differencing / collinear,
  two-pass Gram-Schmidt on the differenced columns; see `native/tseries_ctsa.h`).
  The **false positives cannot**: the input is legitimate, ctsa rejects it, and
  no transformation of the regressors that keeps the model identical changes a
  sample-and-hold structure (any linear recombination of held columns is still
  held). The only fix is the argument order.
- **Why the swap is correct:** a column-major N×k buffer *is* the row-major
  k×N transpose, whose rank is the same. `rank()` with `M = k < N` transposes
  it internally and runs the SVD on the true N×k design.
- **Behaviour on normal input:** unchanged. A design of full rank was accepted
  before and is accepted now; only the two defective verdicts change. Covered by
  `test/sarimax_test.dart` ("regressors sampled at half rate are not falsely
  rejected"), which fails without this patch.
- **Not patched, on purpose (now moot):** upstream `ctsa.c:1209` had the same
  pattern (`rank(xreg,N,r)`) in the auto-ARIMA path with regressors. It was
  never reachable for us, and in 0.8.0 it left the tree with the auto-ARIMA
  cluster (entry 9).
- **Also known, not patched:** `rank_c()` (matrix.c) returns early without
  freeing `U`, `V`, `q` when its SVD fails to converge (the fit then reports
  `retval` 7, which tseries surfaces as a `TseriesNumericException`). Measured
  on a host fuzz of `rank()` (2026-10-04, 200 000 random N×k designs per row):
  | entries of one design span | SVD failures |
  | --- | --- |
  | ≤ 100 orders of magnitude (incl. 1e150-scale columns beside unit ones) | 0 |
  | ~200 orders | 11 |
  | ~300 orders | 31 |
  | ~600 orders (1e±300) | 4 832 |
  The shim rejects any value beyond ±1e150 (`TSERIES_MAX_ABS_VALUE`), which
  removes the overflow-driven failures but not, in principle, a design mixing
  1e-150 with 1e150. Such data is not a measurement of anything, the leak is a
  few KB per such call, and the fix would be a second patch to a function we
  otherwise never touch — so it is documented, not patched. Revisit if a real
  input ever lands there.
- **On version bump:** check the call sites; if upstream fixed the layout,
  drop this patch and keep the regression test.

### 4. `src/dist.c` — LGPL sections removed (2026-10-05)

- **What:** the file's tail from the `// SPDX-License-Identifier:
  LGPL-3.0-or-later` marker before `r8_max()` to the end of the file (upstream
  lines 760-1034: `r8_max()` and `betainv()`, Remark AS R19 / Algorithm AS 109)
  was cut and replaced by a `tseries patch (license)` comment. Nothing else in
  the file changed; the prototypes in `dist.h` are unchanged.
- **Why not on our side:** the file is mixed-license. Excluding the whole
  translation unit would also drop its BSD/CC0 functions (`ibeta`,
  `gamma_log`, ...) that the ARIMA path needs, and a second definition of the
  same symbols elsewhere would not link. The section boundaries were located by
  grepping for the SPDX markers and function headers, without reading bodies.
- **Replacement:** `third_party/ctsa_replacements/betainv.c`.
- **On version bump:** re-apply the cut at the LGPL markers; check whether
  upstream added callers of `betainv`/`r8_max`.

### 5. `src/lls.c` — unattributed MIT routines removed (2026-10-05)

- **What:** `svd_gr()`, `svd_gr2()` and `minfit()` (upstream lines 478-1340,
  three `SPDX-License-Identifier: MIT` sections) were cut and replaced by a
  `tseries patch (license)` comment. `lls.h` is unchanged.
- **Why:** MIT requires reproducing the copyright notice, and none exists. The
  MIT tags were added by the 2024 SPDX audit (ctsa PR #15, by bavay), not by
  an author; upstream's `components_copyright/MIT_COPYRIGHT` is an unfilled
  template (`<YEAR> <COPYRIGHT HOLDER>`); the file's history (initial commit by
  Rafat Hussain, 2014-10-30; later commits also his, plus the SPDX commit)
  names no other author. The code is a C translation of the EISPACK Fortran
  translation of Golub & Reinsch's ALGOL `SVD`/`MINFIT` (Num. Math. 14,
  403-420, 1970; EISPACK itself is distributed by netlib without a copyright
  notice) — the same lineage and the same header comment as ctsa's
  BSD-3-Clause `svd()` in `matrix.c`, which suggests, but does not prove, the
  same translator. No holder could be established, so the code was removed
  rather than attributed.
- **Reachability:** `svd_gr`/`svd_gr2` have no callers anywhere. `minfit` is
  called only by `lls_svd2()`, reached only from `regress()` when its
  least-squares method is `"svd"`; the default is `"qr"` and `setLLSMethod()`
  has no caller in ctsa or the shim — it was linked, never executed.
- **Replacement:** `third_party/ctsa_replacements/minfit.c`, a thin
  BSD-3-Clause wrapper over `svd()` from `matrix.c` with the same documented
  contract. Verified against `numpy.linalg.lstsq` through `lls_svd2()`.
- **On version bump:** re-apply the cut; re-check that nothing new calls the
  removed functions.

### 6. `src/talg.c` — `ppsum()` removed (2026-10-05)

- **What:** upstream lines 1171-1194 (`ppsum()`, whose body carries
  "Copyright (C) 1997-2000 Adrian Trapletti" — R package tseries, GPL) were cut
  and replaced by a `tseries patch (license)` comment. `talg.h` is unchanged.
- **Why not on our side:** `talg.c` also holds BSD code on the fit path
  (`diff`, `poly`, ...), so the unit cannot be
  excluded; a second definition elsewhere would not link.
- **How it was located:** the function's first line and its closing brace
  were found by a line map that printed only line numbers and the classes
  "blank / text / closing brace" — the body was not displayed.
- **Replacement:** `third_party/ctsa_replacements/ppsum.c`.
- **On version bump:** re-apply the cut. (Its upstream callers were in
  `unitroot.c`, deleted in 0.8.0 — see entry 9; nothing in the tree calls
  `ppsum` any more, the replacement is kept for a future caller.)

### 7. Unreachable code of restrictive origin removed (2026-10-05)

Nothing in ctsa or tseries calls any of the following (the gc-sections link
already dropped all of it from the binary), but the package publishes the
**source**, so it was cut and replaced by `tseries patch (license)` comments.
Prototypes in the headers are left as they are (no definitions, no callers).
Boundaries were found with the same numbers-only line map as in entry 6.

| File | Removed | Origin |
| --- | --- | --- |
| `src/matrix.c` | `tred2`, `pythag`, `tqli`, `eigensystem` (upstream 1499-1663) | "Modified version of Numerical recipes tred2"; tqli/pythag the NR routines. The NR licence forbids redistributing its source. |
| `src/emle.c` | `flikam`, `fas197`, `as197` (upstream 3045-3346) | Applied Statistics AS 197 (Mélard 1984), RSS "no fee" terms. Not read. |
| `src/talg.c` | `twacf` (upstream 456-614) | autocovariance routine used only by `flikam` (AS 197). Not read. |
| `src/talg.c` | `artrans`, `arinvtrans`, `gradtrans` (616-740), `transall`, `invtransall` (955-1003) | C versions of `partrans`/`invpartrans`/`ARIMA_Gradtrans` of R `stats/src/arima.c` (GPL-2+): `partrans` and `ARIMA_Gradtrans` compared line by line (same loops, same `w1[100]`/`w2[100]`/`w3[100]` work arrays). |

- **On version bump:** re-apply; re-check that nothing new calls them (the
  MSVC link fails on an unresolved external even from unreachable code, so a
  new caller shows up at build time).

### 8. `src/talg.c` — `archeck()` and `invertroot()` removed (2026-10-05)

- **What:** upstream lines 467-678 of the current tree (`archeck()`, 467-531,
  and `invertroot()`, 533-678) were cut and replaced by a
  `tseries patch (license)` comment. `talg.h` is unchanged.
- **Why:** both are C versions of R's `arCheck()`/`maInvert()`
  (`stats/R/arima.R`, GPL-2+) — see "License audit" — and run on every fit
  (`checkroots_cerr()` in `emle.c` after CSS, and `invertroot()` directly in
  the MLE paths).
- **Why not on our side:** same as entry 6 — `talg.c` holds other BSD code on
  the fit path, and a second definition would not link.
- **How it was located:** by a grep for function headers only and a line map
  printing only line numbers and the classes "blank / top-level text / closing
  brace"; the bodies were not displayed (the cut itself was a `sed` line-range
  delete). Before the cut, the old object was compiled with
  `-Darcheck=old_archeck -Dinvertroot=old_invertroot` for black-box tests.
- **Replacement:** `third_party/ctsa_replacements/arma_roots.c` (clean-room
  round 3 below).
- **On version bump:** re-apply the cut; check that the callers in `emle.c`
  (`checkroots`, `checkroots_cerr`, the `invertroot(q, tf + p)` calls) still
  use the same conventions — re-run `tool/arma_roots_accuracy.c` in record mode
  against the new upstream.

### 9. Auto-ARIMA cluster removed (2026-10-05, released in 0.8.0)

- **What:** deleted with `git rm`: `src/autoutils.c/.h`, `src/unitroot.c/.h`,
  `src/seastest.c/.h`, `src/stl.c/.h`. Cut and replaced by
  `tseries patch (license)` comments (line numbers are of the tree at
  `fb62800`):

  | File | Removed |
  | --- | --- |
  | `src/ctsa.c` | `auto_arima_init` (229-347), `auto_arima_exec` (428-545), `arima2`, `myarima`, `search_arima`, `set_results`, `newmodel`, `auto_arima1` (671-2182), `auto_arima_setMethod` (2420-2433), `auto_arima_setOptMethod` (2661-2707), `auto_arima_setApproximation` .. `auto_arima_setVerbose` (2733-2828), `auto_arima_predict` (3042-3161), `auto_arima_summary` (3700-3833), `myarima_summary`, `aa_ret_summary` (3965-4103), `myarima_free`, `aa_ret_free`, `auto_arima_free` (4391-4408) |
  | `src/ctsa.h` | `auto_arima_object`/`struct auto_arima_set` and `auto_arima_init` (145-209), `myarima_object`/`struct myarima_set`/`myarima` (227-240), `aa_ret_object`/`struct aa_ret_set`/`auto_arima1` (242-252), `search_arima` (254-255), and the prototypes of every removed function; `#include "autoutils.h"` replaced by `#include "errors.h"` (the only header it pulled in transitively that still exists) |
  | `src/talg.c` | `supsmu`, `supsmu_`, `smooth_` and their f2c common blocks (651-1138, end of file) |
  | `src/talg.h` | `supsmu` prototype; `#include "stl.h"` replaced by `#include "polyroot.h"` (what it pulled in transitively; `ctsa.c` needs `polyroot`) |

  `hook/build.dart` no longer lists the four deleted units; the shim's
  `tseries_auto_arima_forecast` and the Dart `autoArimaForecast` are gone.
- **Why:** licence — see "License audit" and "Auto-ARIMA removal" below. The
  user decided (2026-10-05) to keep ctsa for fixed-order fits, publish under
  BSD-3-Clause, and write order selection in Dart (`autoArima`).
- **Why not on our side:** excluding units from the build is not enough —
  the package publishes the source tree — and `ctsa.c`/`talg.c` also hold the
  fixed-order code, so the clusters had to be cut out of them.
- **Kept on purpose:** `sarimax_wrapper*` (`ctsa.c`) and `boxcox.c`. Neither
  is reachable from the shim any more (the clang gc-sections link drops both),
  but neither is in the audited auto-ARIMA set, and the vendored-code policy
  is to remove only what has to go. See "Open licensing questions" item 3.
  (Both were removed later the same day — entry 11.)
- **Effect on the fixed-order path: none.** Verified 2026-10-05: (1) a dump of
  every double returned by 635 fixed-order fits (`fitSarima` all three
  methods, `fitSarimax` all three methods x mean x drift, the four
  statsmodels-reference SARIMAX cases, the health check) on 6 series and 9
  orders is byte-identical before and after (SHA-256 equal, the 48 expected
  exceptions identical); (2) on arm64, armv7 and x86_64 (NDK 29, -O3,
  gc-sections) every function that survives in the library has the same size
  before and after — only functions reachable solely from the auto cluster
  disappeared; (3) all package tests pass with unchanged tolerances.
- **On version bump:** re-apply the cuts (the function names above are the
  map); re-check with a gc-sections link and `--no-undefined` that no kept
  code calls into the removed set, and that MSVC links (it fails on an
  unresolved external even from unreachable code).

### 10. `src/emle.c` — exact likelihood of pure AR(1) error models (2026-10-05)

- **What:** in the three AS 154 objective functions `fas154()`,
  `fas154_seas()` and `fas154x_seas()`, the special case
  `(ip == 1 && iq == 0) || (ip == 0 && iq == 0)` now sets `iupd = 1` (one
  added line each, marked `tseries patch (AR(1) likelihood)`). `ip`/`iq` are
  the AR/MA lengths **after** the seasonal polynomials are multiplied out
  (`p + s·P`, `q + s·Q`).
- **Defect:** that branch sets the initial state variance to the stationary
  `P0 = 1/(1 − φ²)` (as AS 154's own note for AR(1) prescribes, because
  `starma()` refuses `ir == 1`) but left `iupd = 0`. `karma()` with
  `iupd == 0` runs a prediction step before the first observation; for
  `ir == 1` that step rebuilds `P` from `V = 1`, so `P0` was thrown away. The
  first observation got variance σ² instead of σ²/(1 − φ²) and the likelihood
  lost ½·ln(1 − φ²). The objective was wrong, so φ̂, σ² and the
  log-likelihood were all off. `P0` from this branch is already a one-step
  prediction variance (as `starma()`'s output is, and the general branch uses
  `iupd = 1` for exactly that reason), so skipping the first prediction is the
  correct fix. `forkal()` (forecasts) already did this right (`ir == 1` →
  `P0 = 1/(1 − φ²)`, `iupd = 1`).
- **Who was affected:** every MLE / CSS-MLE fit (`fitSarima`, `fitSarimax`
  with `cssMle`/`mle`, `arimaForecast`, `sarimaxForecast`) whose error model
  multiplies out to `ip == 1, iq == 0`: ARIMA(1,d,0) with no seasonal AR/MA —
  any d, any seasonal differencing D (e.g. (1,1,0)(0,1,0)[12]), with or
  without mean, drift or regressors — and the odd (0,d,0)(1,D,0)[1]. CSS-only
  and Box-Jenkins fits, and fits where CSS stopped first (retval 4/10/12), never
  reached it. Size on the reference matrix (28 in-branch cases): log-likelihood
  off by 2e-7 … 134 — up to 2.6 with an estimated mean or after
  differencing, large for an undifferenced AR(1) without a mean and φ near 1
  (Nile: 16.4; synthetic φ = 0.995, n = 200: 134); σ² too large by up to
  ×3.9; φ̂ off by up to 0.07. Sweet Limit's ladder (2,1,1)/(1,1,1)/(0,1,1)
  never reaches the branch.
- **`ip == 0 && iq == 0`** (white noise after differencing) shares the branch:
  there φ[0] = 0, so `P0 = V = 1` and both `iupd` values give the same
  numbers; verified bit-identical (below).
- **Not fixable on our side:** the objective function is internal to the
  optimiser loop in `as154_seas()`/`as154x()`; neither the shim nor the
  wrapper can change what it evaluates. Recomputing the log-likelihood in Dart
  was rejected (Arina, 2026-10-05): φ̂ and the forecast would come from one
  objective and the log-likelihood from another.
- **Audit of the AS 154/182 branches that depend on ip/iq** (each has a test
  case; "matrix" = `test/fixtures/ar1_reference.dart`, "dump" = the 6 787-fit
  bit dump below):

  | Where | Condition | Status | Exercised by |
  | --- | --- | --- | --- |
  | `fas154*` special branch | `ip == 1 && iq == 0` | **defect, fixed here** | 28 matrix cases: AR(1) ±μ, +regressor, (1,1,0)±drift, (1,0,0)(0,1,0)[12] ±drift, (1,1,0)(0,1,0)[12], (0,0,0)(1,0,0)[1] |
  | `fas154*` special branch | `ip == 0 && iq == 0` | correct either way (P0 = V = 1) | dump: (0,0,0), (0,1,0), (0,1,0)(0,1,0)[12]; random-walk test |
  | `fas154*` else → `starma`, `iupd = 1` | everything else | correct | matrix: AR(2), MA(1) ±μ, (0,1,1), (1,0,0)(1,0,0)[12], (1,1,0)(0,1,1)[12]; dump: 17 more orders |
  | `starma` | `ip == 0` (MA: back-substitution) / `ip != 0` (AS 75 system) | correct | MA(1), (0,1,1) / AR(2), ARMA, seasonal AR |
  | `starma` | `ir == 1` → ifault 8; `ip = iq = 0` → ifault 4 | unreachable from `fas154*` (special branch first) | — |
  | `karma` | `ir != 1` (state shift and update) | correct for `ir == 1` (scalar state) | all AR(1) cases |
  | `karma` | `*nit != 0` (fast recursions) | **unreachable** (every caller passes `nit = 0`); **latent defect:** its MA loop increments `i` instead of `j` (`for (j = 0; j < iq; ++i)`) and would never end. Not touched. | — |
  | `forkal` | `ir == 1` → P0 = 1/(1 − φ²), `iupd = 1` | correct | AR(1) forecasts/SEs vs closed form and statsmodels |
  | `forkal` | `ip² + iq² == 0` → ifault 4, returns before writing forecasts/MSEs | **defect, NOT fixed** (reported 2026-10-05): every (0,d,0)(0,D,0) forecast standard error is uninitialised memory | dump: 660 fits whose SEs differ run to run |
  | `forkal` | `id != 0`, `id >= 2`, `id != 1` (differencing) | correct | d = 1, d = 2, D = 1 (id = 13) |
  | `as154`, `as154_seas` | `p+q(+P+Q) == 0 && d(+D) == 0` with CSS-ML → return after CSS | correct (white noise: CSS = exact) | (0,0,0)+μ |
  | `as154*`, `as154x` | CSS retval ≠ 1 (`as154x`) or CSS roots rejected (10/12) → return before MLE | by design | 28 in-branch dump fits stayed bit-identical for this reason (retval 4/10/12) |

- **Verification (2026-10-05).**
  1. *Function at fixed parameters:* a C harness calling `fas154`,
     `fas154_seas` and `fas154x_seas` directly (all three, before and after),
     at statsmodels' optimum and at tseries' optimum of every matrix case,
     against statsmodels `SARIMAX(...).loglike` (stationary initialisation,
     σ² concentrated) on the same series: after the patch |Δll| ≤ 5.8e-12·
     max(1, |ll|) in all 68 points, identical in the three functions; before,
     off by 2e-7 … 134 in exactly the in-branch cases and unchanged elsewhere.
  2. *Optimum* (`fitSarima`, `fitSarimax` cssMle and mle; 89 fits): ll (minus
     the π offset) ≥ statsmodels' − 3.4e-5 (max shortfall, gate 1e-4; the
     shortfall is the reporting offset described below, not the optimum); ARMA
     coefficients within 1.6e-5 (gate 1e-4); σ²/σ²_sm within 1 ± 2.3e-5;
     retval 1 everywhere.
  3. *Forecasts and SEs* (h = 1..10) at tseries' parameters against
     statsmodels `get_forecast` (exact diffuse initialisation for the
     differenced part): ≤ 8.6e-10 relative, ≤ 1e-15 for the AR(1) cases.
  4. *Locality:* 6 787 fits (11 series × 26 orders × every method, mean,
     drift and regressor combination, plus the Sweet Limit ladder through
     `arimaForecast`/`sarimaxForecast`), every returned double dumped as its
     64-bit pattern before and after. All 4 928 fits that neither reach the
     branch nor are white noise are identical (SHA-256 of the dump equal).
     1 105 in-branch fits changed; 28 in-branch fits did not, all with
     retval 4/10/12 (no MLE step). The 726 white-noise fits are identical
     except the forecast SEs, which are uninitialised memory before and after
     (see the `forkal` row above). Library: on arm64/armv7/x86_64 (NDK 29,
     -O3, gc-sections) only `fas154_seas` and `fas154x_seas` changed size.
  5. *Tests:* `test/ar1_likelihood_test.dart` (closed-form exact AR(1)
     likelihood, σ² and ψ-weight forecast/SE at the returned parameters, and
     statsmodels' optimum from `test/fixtures/ar1_reference.dart`, generated
     by `tool/ar1_likelihood_reference.py`). Against the unpatched library 139
     of its tests fail, in all 28 in-branch cases and in none of the 6 others.
- **Not changed here:** `ctsa.c`'s `log(2 * 3.14159)` (every log-likelihood
  is still +½·T·ln(π/3.14159) high; separate decision). Also unchanged at the
  time and worth knowing: ctsa reported the log-likelihood and σ² of the
  *last* objective evaluation, which is inside its finite-difference Hessian,
  i.e. a step (≈ 6e-6 relative) away from the returned estimate — observed up
  to 3.4e-5 in log-likelihood and 2.3e-5 relative in σ² (MSVC build). Fixed
  later by entry 13.
- **On version bump:** check whether upstream sets `iupd = 1` in the three
  special cases (proposed upstream: `doc/upstream/ctsa-ar1-likelihood.md`);
  if not, re-apply and re-run `test/ar1_likelihood_test.dart` — it fails
  without the patch.

### 11. `src/ctsa.c`, `src/ctsa.h` — `sarimax_wrapper*` removed; `boxcox.c` deleted (2026-10-05)

- **What:** cut from `ctsa.c` and replaced by `tseries patch (license)`
  comments (line numbers of the tree at `1f9c9c4`): `sarimax_wrapper()`
  (346-434), `sarimax_wrapper_predict()` (1145-1265),
  `sarimax_wrapper_summary()` (1683-1811), `sarimax_wrapper_free()`
  (2096-2099). From `ctsa.h`: `sarimax_wrapper_object`,
  `struct sarimax_wrapper_set` and the four prototypes. `git rm`
  `src/boxcox.c`, `src/boxcox.h` — their only remaining user was
  `sarimax_wrapper()` (`boxcox_eval`). `talg.h`'s `#include "boxcox.h"` is
  replaced by `#include "optimc.h"` and `#include "stats.h"`, what it pulled in
  transitively. `hook/build.dart` no longer lists `boxcox.c`.
- **Why:** licence hygiene for the published source tree. The group was never
  audited (its drift rule resembles R forecast's `Arima()`; "Open licensing
  questions" item 3) and tseries never called it: the shim uses
  `sarimax_init`/`exec`/`predict`/`free` directly, and the clang gc-sections
  link already dropped it. `boxcox.c` also carried known undefined behaviour
  (`boxcox()` reads an uninitialised variable; entry 1's missing `return`).
- **Kept:** `fminbnd()` in `optimc.c` (now without a caller in the tree) and
  therefore `ctsa_replacements/brent_local_min.c`/`brent.h`: MSVC resolves
  references from unreachable code, so `fminbnd` still needs
  `brent_local_min`.
- **Effect on the fixed-order path: none.** The 6 787-fit dump after this
  change equals the one taken after entry 10 for every fit except the
  white-noise forecast SEs (uninitialised, see entry 10); per-function sizes
  on arm64/armv7/x86_64 unchanged; links with `-Wl,--no-undefined` on all
  three ABIs; exports still the four `tseries_*` functions; the MSVC build
  links (Windows DLL 470 528 B, unchanged).
- **On version bump:** re-apply the cuts by name; re-check with
  `--no-undefined` and an MSVC link.

### 12. Shim only — memory-safety defects worked around on our side (2026-10-05)

**No vendored file changed.** Recorded here because the shim now replaces the
role of two vendored functions and refuses inputs on which vendored code is
unsafe; on a version bump each item must be re-checked against upstream.

Found by building ctsa + shim on the host with AddressSanitizer (Windows clang
22; Linux clang 14 with ASan+UBSan) and MemorySanitizer (Linux) and running
every order p,q ∈ 0..2, P,Q ∈ 0..1, d ∈ 0..2, D ∈ 0..1, s ∈ {1, 12}, all
methods, horizons 0/1/6/24, 0–2 regressors ± intercept, n ∈ {8, 13, 15, 20,
26, 30, 61, 150}.

- **`sarima_predict()` / `sarimax_predict()` (`ctsa.c`) — no longer called.**
  1. They copy `resid[i] = obj->res[i]` for `i < N`, but `obj->res` holds only
     `N - d - s*D` values: a heap **over-read** of up to `d + s*D` doubles
     past the fit object whenever that exceeds the (coefficients)² of the
     `vcov` block that follows — every differenced white-noise order, and
     e.g. (0,0,1)(0,1,0)[12]. The values are never used (`karma()` overwrites
     `resid[0 .. N-id)` before `forkal()` reads it), so it is not a numerical
     error — but it is a crash: in a Dart stress run of white-noise
     forecasts the VM died with an access violation whose address resolves
     to this load in `sarimax_predict()` (ctsa.dll RVA 0x11CC0, the
     `resid[i] = obj->res[i]` right after the regressor loop): the read ran
     off the end of committed heap memory. 1 run in 6 of ~32 000 fits.
  2. `forkal()` (AS 182) returns `ifault = 4` for `ip == iq == 0` **without
     writing** the forecasts or MSEs; the wrappers ignore the code. Every
     white-noise order (p = q = P = Q = 0: (0,d,0)(0,D,0), incl. the random
     walk) therefore returned **uninitialised forecast MSEs** (random
     standard errors, sometimes NaN → `TseriesNumericException`) and, for
     d + D > 0, a point forecast equal to whatever the caller's buffer held
     (0 from Dart's zeroed arena; the mean was added correctly for d = D = 0).

  The shim's `shim_forecast()` repeats their preparation statement for
  statement (same floating-point order) with a zero-filled residual buffer,
  checks `forkal()`'s return code, and passes a white-noise model as the AR(1)
  model with φ = 0 — the same state-space form AS 182 already uses for every
  `ir == 1` model (`P = 1/(1 − φ²) = 1`, no `starma`), so point forecasts,
  MSEs and σ² (forkal's own mean squared one-step residual) come from the same
  algorithm as for every other order. Verified: forecasts/SEs equal the
  closed form (σ²·Σψ², ψ the weights of 1/((1 − B)^d (1 − B^s)^D)) to ≤ 4e-15
  relative, and statsmodels to the same where nothing but σ² is estimated
  (to 5e-10 when ctsa's numerically optimised mean enters;
  `test/white_noise_test.dart`).
  Non-white-noise forecasts are bit-identical to before.
- **`cssx()`/`fcssx()` (`emle.c`) heap overflow — refused.** The SARIMAX CSS
  step zeroes `obj->x[offset + N + i]` for `i < p + s*P` in a buffer of
  `2r + (3 + M)·N` doubles (N differenced, M regressors incl. intercept):
  once `p + s*P > (2 + M)·N` this **writes past the heap block** (e.g.
  (0,0,0)(1,1,0)[12] on 15 points). Before that point the CSS sum is empty
  and the fit fails with NaN. The shim now returns `TSERIES_ERR_BAD_SIZE`
  (→ `ArgumentError`) for `fitSarimax` with `cssMle`/`css` whenever
  `n − d − s·D ≤ p + s·P`. `mle` (no CSS step) is not restricted.
- **Box-Jenkins (`nlalsms`, `boxjenkins.c`) — three guards.**
  - p = q = P = Q = 0 with d + D > 0: nothing to estimate (no mean after
    differencing); `nlalsms` solves a 0×0 system and `linsolve()` (`matrix.c`)
    **writes out of bounds** → heap corruption, then a crash (not `exit()`
    as previously recorded). Refused with `TSERIES_ERR_BAD_PARAM` → a
    dedicated `ArgumentError`. White noise *with* a mean (d = D = 0) still
    fits.
  - `p + s·P` or `q + s·Q` larger than the differenced length: `avaluem()`
    **reads past** the series/residual buffers. Refused (`BAD_SIZE`).
  - s = 0 (non-seasonal): `USPE_seasonal()` reads `cv[0]` of a zero-length
    buffer on every non-seasonal Box-Jenkins fit. The shim passes s = 1
    instead when P = D = Q = 0 (the period is otherwise unused); results
    unchanged.
- **Left as is (no memory error; recorded):** `karma()`'s `*nit != 0` branch
  increments `i` instead of `j` in `for (j = 0; j < iq; ++i)` — an infinite
  loop, unreachable because every caller passes `nit = 0`. Two
  uninitialised-value reads remain (MemorySanitizer, not AddressSanitizer):
  `bfgs_min()` (`secant.c`) reads its gradient buffer uninitialised when the
  first objective value is NaN — reached by SARIMA MLE/CSS fits whose
  differenced series is not longer than p + s·P (CSS has no residual), which
  makes those fits' outcome **non-deterministic** run to run; and
  `USPE_seasonal()` reads uninitialised autocovariances on seasonal
  Box-Jenkins fits of short series, same consequence. Neither is reachable
  from `arimaForecast`'s length floor for non-seasonal orders; both are
  reported in `doc/upstream/ctsa-other-defects.md`. Both are refused since
  entry 15.
- **Effect:** every fit outside white noise and outside the refused inputs is
  bit-identical (old vs new dump of 238 000 fit records, n ∈ {8, 13, 15, 20, 30, 61, 150}; the only
  other differences fall in the two non-deterministic families above and
  occur equally between two runs of the *same* binary). The sanitizer runs
  are clean (no AddressSanitizer or UBSan report) on everything not refused.
- **On version bump:** re-check whether upstream fixed `forkal`'s ip = iq = 0
  case, the `res` copy, `fcssx`'s conditioning buffer and `nlalsms`; if so,
  the shim may call `*_predict()` again (re-run `test/white_noise_test.dart`
  and the dump).

### 13. `src/emle.c` — σ², log-likelihood and residuals evaluated at the estimate (2026-10-05)

- **What:** one added statement, `FUNCPT_EVAL(&<objective>, tf, pq);`, in each
  of the six estimators that compute a finite-difference Hessian, placed after
  `hessian_fd()` and the covariance inversion and immediately before
  `*var = …` (marked `tseries patch (estimate-point likelihood)`).
- **Audit — every function that reports σ²/log-likelihood after a Hessian:**

  | Function | Objective re-evaluated | Reaches the caller as | Reachable from tseries |
  | --- | --- | --- | --- |
  | `as154_seas()` | `fas154_seas` | `fitSarima` MLE: σ², `logLikelihood`, `aic` | yes |
  | `css_seas()` | `fcss_seas` | `fitSarima` CSS: σ², `logLikelihood`; `fitSarima` MLE when CSS ends the fit (retval 10/12, and the no-parameter mean-only model): σ² | yes |
  | `as154x()` | `fas154x_seas` | `fitSarimax` cssMle/mle: σ², `logLikelihood`, `aic` | yes |
  | `cssx()` | `fcssx` | `fitSarimax` css: σ², `logLikelihood`; cssMle when CSS ends the fit (retval ≠ 1 from CSS, 10/12): σ² | yes |
  | `as154()` | `fas154` | ctsa's `arima_exec`/`ar_exec` | no (not called by the shim) |
  | `css()` | `fcss` | ctsa's `arima_exec`, `as154()`'s CSS start | no |

  No other function in `emle.c` (or elsewhere in the vendored tree) reads
  `ssq`/`loglik` after a Hessian: `nlalsm()`/`nlalsms()` (Box-Jenkins) do not
  report a log-likelihood, and the forecaster `forkal()` computes its own σ²
  from its own one-step residuals (`emle.c`, "Calculate M.L.E. of sigma
  squared") — it does not receive the fit's σ².
- **Defect:** the objective functions store `ssq`, `loglik` and the residuals
  in the likelihood object as a side effect, and the six functions copy them
  out after `hessian_fd()`. The last evaluation was therefore one of the
  Hessian's finite-difference points, not the estimate `tf`. Which one is
  decided by the order in which the compiler evaluates the four calls in
  `((f(xij) - f(xi)) - (f(xj) - f(x)))`, which C leaves unspecified:
  - **clang and GCC** (Android NDK, iOS/macOS, Linux; checked with Windows
    clang 22, Linux clang 14 and GCC) evaluate `f(x)` last — `x` *is* `tf` —
    so the reported numbers were already exact; the patch is a bit-for-bit
    no-op there (whole dumps identical before and after).
  - **MSVC** (Windows desktop builds) evaluates `f(xi)` last: the point
    `tf + δ·e_last`, δ = ε^(1/3)·max(|tf_last|, 1/dx_last), ε^(1/3) ≈ 6.06e-6
    (`dx` = 1, or 10 × the OLS standard error for a regression column).
    Proven, not inferred: re-evaluating the objective at exactly that point
    reproduces the unpatched MSVC σ² and log-likelihood bit for bit in
    720/720 Sweet Limit ladder fits, 140/154 reference fits (the rest have a
    regression coefficient last) and 3 488/3 552 fits of the broad fixture
    matrix (the rest are boundary fits: MA ≈ ±1, or Lake Huron's level fitted
    without a mean).
- **Size of the change (MSVC only; nothing else moves).** Reference set
  (24 comparability cases + the 34-case AR(1) matrix, all fit paths, 154
  fits): |Δll| ≤ 5.6e-5, |Δσ²|/σ² ≤ 2.3e-5. Sweet Limit ladder
  ((2,1,1)/(1,1,1)/(0,1,1), MLE, windows of 61–144 points of CGM-like data,
  720 fits): |Δll| ≤ 4.2e-5, |Δσ²|/σ² ≤ 1.2e-4 (θ ≈ −0.9994), ≤ 1.8e-3 on
  fits that CSS ends with retval 10 (σ² only; their log-likelihood is not
  reported). Broad matrix (7 series × 26 orders × every method, mean, drift
  and regressor combination, 4 150 fits): 3 438 changed, up to 62 in ll and
  72 % in σ² on degenerate fits (a level series fitted with no mean, CSS
  with regressors) — still the Hessian neighbour (above), the objective is
  just extremely steep there. Zero-parameter fits: bit-identical.
- **Not fixable on our side:** the values are produced inside `as154_seas()`/
  `as154x()` from their private likelihood object; the shim only sees the
  copies. Recomputing the likelihood outside would duplicate ctsa's
  objective (and its data preparation) in the shim.
- **Verification (2026-10-05).**
  1. *Locality* — every output double dumped as its bit pattern, unpatched vs
     patched, MSVC: the Dart dump (6 787 fits through `fitSarima`,
     `fitSarimax`, `arimaForecast`, `sarimaxForecast`; 11 series × 26
     orders) and the C dumps above (5 024 fits): status, retval,
     coefficients, mean, covariance and coefficient SEs, point forecasts and
     **forecast SEs** bit-identical; only σ², log-likelihood and AIC change
     (6 314 of 6 787 Dart fits). Synthetic order × length matrix (MSVC,
     n ∈ {8, 15, 30, 61, 150}, d ≤ 2, D ≤ 1, s ∈ {1, 12}, p, q ≤ 2,
     P, Q ≤ 1, both fit paths, every method, h ∈ {0, 1, 6, 24}; 170 100
     fits, unpatched run twice): 36 942 identical, 132 666 changed in
     σ²/ll/AIC only, 492 differ elsewhere — all of them in the two
     families that already differ between two runs of the *unpatched*
     binary (SARIMA MLE/CSS with Nd ≤ p + s·P, seasonal Box-Jenkins on
     short series; entry 12), none outside.
  2. *Consistency* — an independent harness rebuilds each likelihood object
     from the series and evaluates the same objective function at the
     returned estimate: σ² and log-likelihood equal the fit's **bit for bit**
     in every fit with retval 1 (5 024 of 5 024 checked, MSVC and clang);
     unpatched MSVC: equal in 3 % of them.
  3. *Oracle* — log-likelihood and σ² at the returned estimate vs an exact
     Gaussian likelihood by Cholesky factorisation of the ARMA covariance
     matrix: ≤ 2.0e-13 and ≤ 6.0e-13 relative (154 reference fits; gates
     1e-8 / 1e-10); vs statsmodels `loglike` at the same parameters ≤ 1.3e-11
     (statsmodels' own σ² is ≤ 3.3e-10 off the Cholesky value on MA models,
     ctsa's ≤ 4e-16). Unpatched MSVC: 5.4e-7 / 2.3e-5.
  4. *Across compilers* — patched MSVC, Windows clang, Linux clang and Linux
     GCC: wherever the estimate is bit-identical, σ² and the log-likelihood
     are bit-identical too (unpatched: up to 2.3e-5 / 1.8e-3 apart). The
     estimates themselves differ between Windows and Linux in ~2 % of fits
     (up to 4e-5 in a coefficient, 2.6e-6 in ll) — different `libm`, not this
     defect.
  5. *Tests* tightened to what is now achieved: `ar1_likelihood_test.dart`
     (closed-form ll 1e-4 → 1e-10 relative, σ² 1e-4 → 1e-12),
     `white_noise_test.dart` (ll 1e-4 → 1e-8, σ² 1e-5 → 1e-10 with a mean),
     `fit_path_comparability_test.dart` (the two fit paths' ll: `< 1e-6` →
     identical; drift 1e-7 → 1e-10; Nile AR(1) 1e-8 → 1e-12 relative).
  6. *The 64 boundary fits not explained by the stencil model* (the 3 552 −
     3 488 of the broad matrix above; checked 2026-10-06 on the current tree,
     i.e. this patch + entries 14–15, MSVC): **consistency** — σ² and
     log-likelihood re-evaluated by the independent harness at the returned
     estimate equal the fit's bit for bit in 64/64; **oracle** — against the
     exact Gaussian likelihood by Cholesky factorisation in 60-digit decimal
     arithmetic (autocovariances from the ARMA linear system, not from a
     truncated ψ-sum) the worst relative error is 1.5e-13 in ll and 4.8e-13
     in σ² (gates 1e-8 and 1e-10); the one CSS fit against an independent
     conditional sum of squares, 1.3e-16 / 2.1e-16. A double-precision
     Cholesky oracle is NOT adequate for 11 of them (Lake Huron fitted without
     a mean: an AR root at φ ≈ 1 − 8e-7, so the covariance matrix is
     ill-conditioned) — it shows up to 7e-10 in σ², which the decimal oracle
     attributes entirely to the oracle. Coefficients, mean, σ² and covariance
     are bit-identical to the patch-A build. Per fit (series | order | method
     | constant/regressors: mean, x = the sin regressor, drift | rel. |Δll| |
     rel. |Δσ²|):

  | Series | Order | Method | Terms | ll | σ² |
  | --- | --- | --- | --- | --- | --- |
  | airlog | (0,0,1)(0,0,0)[0] | cssMle | drift | 0.00e+00 | 3.03e-16 |
  | airlog | (0,0,1)(0,0,0)[0] | cssMle | x+drift | 3.64e-16 | 6.07e-16 |
  | airlog | (0,0,1)(0,0,0)[0] | mle | x | 3.21e-16 | 3.45e-16 |
  | airlog | (0,0,1)(0,0,0)[0] | mle | drift | 2.42e-16 | 0.00e+00 |
  | airlog | (0,0,1)(0,0,0)[0] | mle | x+drift | 1.21e-16 | 4.55e-16 |
  | airlog | (0,0,2)(0,0,0)[0] | mle | mean+x | 1.45e-16 | 1.12e-15 |
  | airlog | (0,0,2)(0,0,0)[0] | mle | x | 2.16e-16 | 8.51e-16 |
  | airlog | (2,1,1)(0,0,0)[0] | cssMle | mean+drift | 0.00e+00 | 2.17e-16 |
  | airlog | (2,1,1)(0,0,0)[0] | cssMle | drift | 0.00e+00 | 2.17e-16 |
  | huron | (1,0,0)(0,0,0)[0] | cssMle | x | 3.56e-15 | 1.72e-13 |
  | huron | (1,0,0)(0,0,0)[0] | cssMle | drift | 1.09e-15 | 2.07e-13 |
  | huron | (1,0,0)(0,0,0)[0] | cssMle | x+drift | 4.05e-15 | 3.49e-13 |
  | huron | (1,0,0)(0,0,0)[0] | mle | mean+x+drift | 0.00e+00 | 2.37e-16 |
  | huron | (1,0,0)(0,0,0)[0] | mle | x | 4.55e-15 | 1.07e-13 |
  | huron | (2,0,0)(0,0,0)[0] | cssMle | mean+x+drift | 2.91e-16 | 2.61e-16 |
  | huron | (4,0,0)(0,0,0)[0] | cssMle | mean+x+drift | 1.46e-16 | 1.31e-16 |
  | huron | (0,0,1)(0,0,0)[0] | cssMle | x | 0.00e+00 | 5.16e-16 |
  | huron | (0,0,1)(0,0,0)[0] | cssMle | drift | 0.00e+00 | 3.42e-16 |
  | huron | (0,0,1)(0,0,0)[0] | mle | mean+x+drift | 1.33e-16 | 0.00e+00 |
  | huron | (0,0,1)(0,0,0)[0] | mle | x | 1.63e-16 | 1.72e-16 |
  | huron | (0,0,2)(0,0,0)[0] | cssMle | mean+x+drift | 0.00e+00 | 2.52e-16 |
  | huron | (0,0,2)(0,0,0)[0] | cssMle | x | 4.30e-15 | 1.51e-15 |
  | huron | (0,0,2)(0,0,0)[0] | mle | mean+x+drift | 1.43e-16 | 3.77e-16 |
  | huron | (0,0,2)(0,0,0)[0] | mle | x | 2.33e-15 | 1.01e-15 |
  | huron | (0,0,2)(0,0,0)[0] | mle | drift | 5.01e-15 | 5.60e-15 |
  | huron | (0,0,2)(0,0,0)[0] | mle | x+drift | 3.42e-15 | 3.23e-15 |
  | huron | (1,0,1)(0,0,0)[0] | cssMle | mean+x+drift | 1.45e-16 | 1.30e-16 |
  | huron | (1,0,1)(0,0,0)[0] | cssMle | — | 1.85e-15 | 1.48e-14 |
  | huron | (1,0,1)(0,0,0)[0] | mle | x+drift | 1.49e-15 | 2.37e-13 |
  | huron | (1,0,1)(0,0,0)[0] | css | — | 1.30e-16 | 2.05e-16 |
  | huron | (2,1,1)(0,0,0)[0] | cssMle | mean+drift | 2.79e-16 | 2.41e-16 |
  | huron | (2,1,1)(0,0,0)[0] | cssMle | drift | 2.79e-16 | 2.41e-16 |
  | huron | (0,1,1)(0,1,1)[12] | mle | mean+x | 2.67e-16 | 4.16e-16 |
  | huron | (0,1,1)(0,1,1)[12] | mle | x | 2.67e-16 | 4.16e-16 |
  | huron | (1,0,0)(1,0,0)[12] | cssMle | mean+x+drift | 0.00e+00 | 0.00e+00 |
  | huron | (1,0,0)(1,0,0)[12] | cssMle | x+drift | 1.49e-13 | 4.79e-13 |
  | huron | (1,0,0)(1,0,0)[12] | mle | — | 7.95e-15 | 3.37e-13 |
  | huron | (0,0,0)(1,0,0)[1] | cssMle | x | 3.56e-15 | 1.72e-13 |
  | huron | (0,0,0)(1,0,0)[1] | cssMle | drift | 1.09e-15 | 2.07e-13 |
  | huron | (0,0,0)(1,0,0)[1] | cssMle | x+drift | 4.05e-15 | 3.49e-13 |
  | huron | (0,0,0)(1,0,0)[1] | mle | mean+x+drift | 0.00e+00 | 2.37e-16 |
  | huron | (0,0,0)(1,0,0)[1] | mle | x | 4.55e-15 | 1.07e-13 |
  | huron | (0,0,0)(1,0,0)[12] | cssMle | drift | 6.59e-16 | 3.07e-13 |
  | lynx | (2,1,1)(0,0,0)[0] | cssMle | mean+drift | 1.22e-16 | 3.02e-16 |
  | lynx | (2,1,1)(0,0,0)[0] | cssMle | drift | 1.22e-16 | 3.02e-16 |
  | lynx | (2,1,1)(0,0,0)[0] | mle | mean+x | 1.22e-16 | 0.00e+00 |
  | lynx | (2,1,1)(0,0,0)[0] | mle | mean+x+drift | 1.22e-16 | 1.51e-16 |
  | lynx | (2,1,1)(0,0,0)[0] | mle | x | 1.22e-16 | 0.00e+00 |
  | lynx | (2,1,1)(0,0,0)[0] | mle | x+drift | 1.22e-16 | 1.51e-16 |
  | lynx | (0,1,1)(0,1,1)[12] | mle | mean+x | 0.00e+00 | 5.93e-16 |
  | lynx | (0,1,1)(0,1,1)[12] | mle | x | 0.00e+00 | 5.93e-16 |
  | nile | (2,0,2)(0,0,0)[0] | cssMle | mean+x+drift | 1.80e-16 | 2.93e-15 |
  | nile | (0,2,1)(0,0,0)[0] | mle | mean+x | 1.77e-16 | 2.58e-16 |
  | nile | (0,2,1)(0,0,0)[0] | mle | x | 1.77e-16 | 2.58e-16 |
  | usacc | (0,0,1)(0,0,0)[0] | cssMle | x | 3.21e-16 | 1.87e-16 |
  | usacc | (0,0,2)(0,0,0)[0] | cssMle | x | 0.00e+00 | 1.80e-16 |
  | www | (0,0,1)(0,0,0)[0] | cssMle | drift | 0.00e+00 | 2.23e-16 |
  | www | (0,0,1)(0,0,0)[0] | cssMle | x+drift | 0.00e+00 | 2.24e-16 |
  | www | (0,0,1)(0,0,0)[0] | mle | mean+x | 2.60e-16 | 1.64e-16 |
  | www | (0,0,1)(0,0,0)[0] | mle | mean+x+drift | 0.00e+00 | 4.07e-16 |
  | www | (0,0,1)(0,0,0)[0] | mle | x+drift | 1.16e-16 | 5.61e-16 |
  | www | (0,0,2)(0,0,0)[0] | cssMle | x | 5.59e-16 | 2.48e-15 |
  | www | (0,0,2)(0,0,0)[0] | cssMle | x+drift | 6.61e-16 | 2.83e-14 |
  | www | (0,0,2)(0,0,0)[0] | mle | mean+x+drift | 3.06e-16 | 3.98e-15 |

- **On version bump:** check whether upstream re-evaluates the objective at
  the estimate after `hessian_fd()`; if not, re-apply the six lines and re-run
  the consistency harness or the tightened tests on an **MSVC** build (on
  clang/GCC the defect is invisible).

### 14. `src/ctsa.c` — exact π in the Gaussian log-likelihood (2026-10-05)

- **What:** the literal `3.14159` in `log(2 * 3.14159)` replaced by
  `3.14159265358979323846` (π to full double precision, the value of `PIVAL`
  in `erfunc.h`), each line marked `tseries patch (pi)`, plus one explanatory
  block comment above `arima_exec()`.
- **All approximations of π in the vendored tree and the shim** (grep for
  `3.14159`, `3.1415`, `6.28318`, and for `ln 2π`/`√(2π)`-derived constants
  `1.8378`, `0.9189`, `2.5066`, `0.39894`):

  | Where | Literal | Feeds logLik/AIC | Action |
  | --- | --- | --- | --- |
  | `ctsa.c` `arima_exec()` MLE, CSS | `log(2 * 3.14159)` | yes (`arima_*`, not used by tseries) | replaced |
  | `ctsa.c` `sarimax_exec()` CSS-MLE, MLE, CSS | `log(2 * 3.14159)` | yes — `fitSarimax` | replaced |
  | `ctsa.c` `sarima_exec()` MLE, CSS | `log(2 * 3.14159)` | yes — `fitSarima`, `arimaForecast` | replaced |
  | `ctsa.c` `ar_exec()` | `log(2 * 3.14159)` | yes (`ar_*`, not used by tseries) | replaced |
  | `dist.c` (2×), `pdist.c` | `3.1415926535897932384626434`, `spi = 0.9189385332046727417803297` | no (distributions) | already exact |
  | `erfunc.h` `PIVAL`, `hsfft.h` `PI2`, `hsfft.c`, `real.c` | full precision | no | already exact |
  | `regression.c` (`pi = 3.141592653589793`, `log(2 * pi)`) | full precision | regression log-likelihood, not reported by tseries | already exact |
  | `ctsa_replacements/polyroot.c` `two_pi`, `erfinv.c` (2.5066… coefficient) | full precision | no | already exact |
  | `native/tseries_ctsa.c` (shim) | none | — | — |

  So the eight lines above were the only inexact π, all in the conversion
  `loglik = -0.5 * (N * (2 * value + 1 + log(2π)))` that `ctsa.c` applies to
  the objective value **after** the optimiser has returned.
- **Defect:** every reported log-likelihood (and AIC = −2·ll + 2k) was too
  high by ½·T·ln(π/3.14159) ≈ 4.22e-7·T (T = N − d − s·D): 6.1e-5 at T = 144,
  8.4e-4 at T = 2000. Irrelevant for comparisons at one T, but it kept
  absolute values from agreeing with R/statsmodels beyond 1e-5 and every test
  and reference script had to subtract it.
- **Not fixable on our side** without two log-likelihoods for one model
  (Arina, 2026-10-05: correcting in Dart was rejected for exactly that
  reason); the shim could subtract the offset itself, but that would hard-code
  ctsa's wrong constant into our code instead of removing it.
- **Verification (2026-10-05), base = the tree after entry 13.**
  1. *Nothing but ll/AIC moves:* bit dumps before/after — the Dart dump
     (6 787 fits, MSVC) and the C harness (reference set 154 fits and Sweet
     Limit ladder 720 fits, MSVC and clang): status, retval, coefficients,
     mean, σ², covariance, forecasts and forecast SEs bit-identical; 6 405
     Dart fits changed only in ll/AIC (every fit that reports an ll). The
     constant enters no objective function, so no estimate can move — and
     none did.
  2. *Size of the shift:* new ll − old ll = −½·T·ln(π/3.14159) to within
     9.5e-14 absolute, i.e. ≤ 7e-15 relative to the log-likelihood (that is
     one rounding of the ll itself; relative to the shift, which is 1e-5…1e-4,
     the same rounding is ≤ 2.4e-9). 871 log-likelihoods checked.
  3. *Oracle without any correction:* log-likelihood at the returned estimate
     vs exact Cholesky likelihood ≤ 2.0e-13 relative, vs statsmodels
     `loglike` ≤ 1.3e-11 (154 reference fits).
  4. *Tests/scripts:* the π correction was removed from
     `ar1_likelihood_test.dart`, `fit_path_comparability_test.dart`,
     `white_noise_test.dart` and the docstring of
     `tool/ar1_likelihood_reference.py`; they compare raw values now.
- **On version bump:** if upstream still uses `3.14159`, re-apply (eight
  lines in `ctsa.c`); `fit_path_comparability_test.dart` (random walk vs
  closed form, 1e-9) fails without it.

### 15. Shim + Dart only — short series that made ctsa non-deterministic are refused (2026-10-05)

No vendored change. Entry 12 left two uninitialised-value reads open
(MemorySanitizer, not memory-unsafe) that make results differ from run to run
— including fits reported as successful. Both are now refused before ctsa
runs: `tseries_sarima_fit()` returns `TSERIES_ERR_BAD_SIZE`, and `fitSarima`
throws the new typed `SeriesTooShortError` (an `ArgumentError`) even before
the native call, from the same rule (`sarimaMinimumLength()` in
`lib/src/domain/internal/arima_support.dart`). `arimaForecast` skips such a
rung with the note "series length n below structural minimum M for <method>:
<rule>" (like its numeric floor), and raises `ArgumentError` only if no rung is
long enough. With Nd = n − d − s·D:

- **SARIMA MLE and CSS: refused when Nd ≤ p + s·P.** Mechanism: both start with
  `css_seas()`, which conditions on the first p + s·P differenced observations;
  none is left, the CSS objective is `0.5·log(0/0)` = NaN at the start point,
  and `bfgs_min()` (`secant.c`) then reads its never-filled gradient buffer.
  MemorySanitizer reports (`inclu2`, `lnsrch`, `grad_calc`, `jrotate`, …)
  appear only in such configurations (128-configuration check at s = 4, every
  report inside Nd ≤ p + s·P). Same rule as the SARIMAX CSS guard of entry 12.
  **MA analogue checked and not adopted:** with p + s·P < Nd ≤ q + s·Q the CSS
  sum is not empty; 1 920 such fits (s ∈ {4, 7, 12}, n = 5…16, p, q, Q ≤ 1,
  d, D ≤ 1, MLE and CSS) gave no MemorySanitizer report and bit-identical
  results in two runs, so no guard was added.
- **Seasonal Box-Jenkins (P + Q > 0): refused when Nd < (P + Q + 1)·s + 2.**
  Mechanism (derived, then measured): `USPE_seasonal()` takes the seasonal
  start values from the autocovariances at lags 0, s, …, (P + Q)·s, computed
  by `autocovar(…, K = (P + Q + 1)·s)`; when K > Nd that routine clamps to
  `M = N − 1` and fills only lags 0 … Nd − 2. So lag (P + Q)·s is
  uninitialised exactly when **Nd ≤ (P + Q)·s + 1** (first use:
  `boxjenkins.c` 240/290, then `linsolve`/`pludecomp`/`nlalsms`). MemorySanitizer
  sweep, unguarded, s ∈ {4, 7, 12}, p, q ≤ 2, P, Q ≤ 1, d, D ≤ 1, every Nd up to
  (P + Q + 1)·s + p + q + 4 (432 configurations): reports in 244, **all** at
  Nd ≤ (P + Q)·s + 1, and for every (P + Q, s) the largest reported Nd equals
  that bound exactly. The guard adds one season of margin: Nd ≥
  (P + Q + 1)·s + 2. Non-seasonal Box-Jenkins is unaffected (its start values
  need lags ≤ p + q, always computed under the general floor).
- **After the guards** (A + π + guards build): MemorySanitizer over the
  whole supported matrix for both entry points — n ∈ {6, 8, …, 61} (16
  lengths), s ∈ {4, 7, 12}, p, q ≤ 2, P, Q ≤ 1, d ≤ 2, D ≤ 1, every method,
  SARIMAX with and without a mean; 93 312 fits (31 104 per period), **zero
  reports**. Determinism on the same matrix (85 536 fits, two runs of one
  binary): unguarded, 1 232 records differed between runs — 720 SARIMA
  MLE/CSS with Nd ≤ p + s·P, 512 seasonal Box-Jenkins below the bound,
  nothing else — and all 1 232 are refused now; guarded, the 65 514 accepted
  fits are bit-identical run to run. Same on the s ∈ {1, 12} matrix of entry
  12 (n ∈ {8, 15, 30, 61, 150}, MSVC; 494 records differ between two unguarded runs,
  0 guarded). The bound was also checked with P, Q ≤ 2 at s = 4 (P + Q up
  to 4): reports again exactly up to Nd = (P + Q)·s + 1.
- **Sweet Limit:** its ladder (2,1,1)/(1,1,1)/(0,1,1), MLE, has structural
  minima 4/3/2 against numeric floors 31/21/11 and app windows of 61–144
  points; nothing it fits is refused (the 720-fit ladder dump is
  bit-identical with and without the guards, and so are the 4 150-fit
  fixture matrix's statuses).
- **On version bump:** if upstream fixes `bfgs_min()`'s NaN start or
  `autocovar()`'s clamp, the guards can be relaxed; re-run the sweep first.

### 16. Shim only — regressor buffer overflow under d > 0 and D > 0 refused (2026-10-06)

**No vendored file changed.**

- **Defect (`emle.c`, `cssx()` and `as154x()`, the same block in both):** `XX`
  is allocated with `N·(r + 1)` doubles, N the fully differenced length
  `n − d − s·D`. With regressors and D > 0 each column is first seasonally
  differenced into `XX + N1·i`, `N1 = n − s·D`, i < ncxreg — up to `N1·ncxreg`
  doubles. That writes past the heap block exactly when **d·ncxreg > N**
  (ncxreg = regressors in use; no intercept column exists when d + D > 0).
  AddressSanitizer (Windows clang 22): *heap-buffer-overflow WRITE in `diffs()`
  (talg.c) from `cssx()` (emle.c:2441)*, e.g. 3 regressors under
  (0,2,0)(0,1,0)[4] on 11 points, every method. Reachable from `fitSarimax`
  with ≥ 3 regressors, d = 2, D ≥ 1 and a short series (the drift column
  alone cannot reach it: drift requires d + D ≤ 1).
- **Workaround:** after screening, `tseries_sarimax_fit()` returns
  `TSERIES_ERR_BAD_SIZE` (→ `ArgumentError` naming the rule) when
  `D > 0 && r_used > 0 && d·r_used > n − d − s·D`. The boundary itself
  (d·r = N) is safe and still fitted (ASan-clean, n = 12 in the same example).
  Test: `test/short_series_guard_test.dart`, group "fitSarimax: regressors
  under d > 0 and D > 0".
- **Why not in the vendored code:** the fix is a one-line allocation size, but
  the inputs that reach it are fully decided by the caller's arguments, so a
  refusal in the shim closes it with no vendored diff (policy above).
- **Found by** a sanitizer matrix that, unlike entries 12/15, includes
  regressors and drift: `tseries_sarimax_fit()` with exactly the Dart
  wrapper's buffer sizes, s ∈ {4, 7, 12}, p, q ≤ 5, P, Q ≤ 2, d, D ≤ 1 (d ≤ 2
  in the targeted run), mean on/off, 0–2 regressors (+ drift), all three
  methods, horizon 2. Completed slices (2026-10-06, guard included): Windows
  clang ASan+UBSan n ∈ {6, …, 10} — 382 968 fits; Linux clang-14 ASan+UBSan
  n ∈ {6, 7, 8} — 260 496 fits; MemorySanitizer n ∈ {6, 7, 8} — 241 056
  fits; earlier unsliced partial runs up to n = 14 (s = 4) / 12 (s = 7) /
  10 (s = 12), ~1.04 M fits. **Zero reports** anywhere except the overflow
  above. Longer series, sampled (p, q ∈ {0, 2, 5}, P, Q ∈ {0, 2}, d, D ≤ 1,
  mean on/off, 0 or 2 regressors ± drift, cssMle and mle; 1 584 fits per
  (s, n), 20-minute budget per sanitizer): Windows ASan+UBSan complete for
  s = 4, n ∈ {12, 16, 24, 40} and s = 12, n = 12 (7 920 fits), s = 12,
  n = 16 partial; MemorySanitizer complete for s = 4, n ∈ {12, 16}, n = 24
  partial. Zero reports. Not covered within the budget: s = 12 at n ≥ 24
  (ASan) / n ≥ 12 (MSan), s = 4 n = 40 (MSan) — seasonal MLE at s = 12,
  P = 2 costs seconds per fit even unsanitized (s = 12, n = 40: 408 s for
  1 584 fits).
- **Effect on accepted fits:** none — the condition implies fewer than ~2
  observations per regressor per order of differencing; the 6 787-fit Dart bit
  dump does not reach it.
- **On version bump:** if upstream sizes `XX` with the undifferenced length,
  drop the guard (re-run the targeted ASan case first).

### 17. `src/initest.c`, `src/initest.h` deleted; their callers cut from `src/ctsa.c`, `src/ctsa.h` (2026-10-06)

- **What:** `initest.c` and `initest.h` deleted (Yule-Walker `ywalg`/`ywalg2`,
  Burg `burgalg`, Hannan-Rissanen `hralg`/`hrstep2`, innovations `innalg`/
  `ma_inn`, `pacf_yw`/`pacf_burg`/`pacf_mle`). Cut from `ctsa.c` and replaced
  by `tseries patch (license)` comments — every function that called them:
  `ar_exec()`, `ar()`, `ar_estimate()`, `model_estimate()`, `pacf_opt()`,
  `pacf()`, `yw()`, `burg()`, `hr()` (no other caller in the tree; nothing
  else called these nine). Their prototypes are cut from `ctsa.h`. `emle.h`,
  `errors.h` and `pred.h` included `initest.h` only for what it included
  itself; each now includes `boxjenkins.h` directly (one line each, marked).
  `hook/build.dart` no longer lists `initest.c`. `ar_init()`, `ar_predict()`,
  `ar_summary()` and `ar_free()` stay (they reference nothing that was
  removed); without `ar_exec()` the `ar_object` API is inert, which does not
  matter to tseries, whose shim never used it.
- **Why:** licence hygiene for the published source tree ("Open licensing
  questions" item 4). The Burg routine is marked upstream as a "C translation
  of Cedrick Collomb's C++ Implementation", and no licence for that original
  was found. The file was never reachable from the shim, but pub.dev
  publishes source, reachable or not. Removing the callers as well is
  required, not cosmetic: MSVC fails the link on an unresolved external even
  when the referencing function is itself unreachable.
- **Effect on the fixed-order path: none.** Stripped `libctsa.so` (NDK 29,
  the hook's flags plus `-Wl,--no-undefined`) is **byte-identical** before and
  after on arm64-v8a (183 320 B), armeabi-v7a (137 392 B) and x86_64
  (239 408 B); the unstripped gc-sections builds before the change already
  contained none of the 19 symbols. The MSVC build links (Windows DLL
  468 480 B, unchanged) and the full test suite passes.
- **On version bump:** delete `initest.*` again, cut the nine callers by name,
  re-point the three `#include "initest.h"` lines; re-check with
  `--no-undefined` and an MSVC link.

### Note (historical, superseded 2026-10-05) — LGPL code on the fixed-order hot path

SARIMAX estimation runs ctsa's `regress()` on every fit (OLS starting values
for the regression coefficients), which calls `tinv()` → `betainv()` from the
**LGPL** section of `dist.c`. No new LGPL translation unit entered the build —
`dist.c` was already compiled and used by auto-ARIMA — but code under that
license is now executed on the fixed-order SARIMAX path, not only by the
order search. This sharpens the existing deferred-LGPL decision above rather
than changing it: a clean-room replacement of `betainv`/`tinv` would now also
be needed for `fitSarimax`. (`tinv` additionally calls `exit(1)` for
`df <= 0`; the shim's degrees-of-freedom guard makes that unreachable.)

### How to see the local diff

```
git log -p -- packages/tseries/third_party/ctsa/src/ packages/tsforge/third_party/ctsa/src/
```

(The package was named `tsforge` through 0.7.0 and lived at `packages/tsforge/`;
history before the rename is under that path. `git log --follow` only follows a
single file, so pass both directories as above, or `--follow` one file.)

Each patch is marked in-source with a `tseries patch` comment (written
`tsforge patch` before the rename) stating what and why, so it is visible when
diffing against a fresh upstream checkout.

## Clean-room replacements (2026-10-05)

`third_party/ctsa_replacements/` holds tseries' own implementations of the
functions removed above, with upstream's names and signatures so no vendored
caller changes. All files: `SPDX-License-Identifier: BSD-3-Clause`,
Copyright (c) 2026, Vladislav Zyuzko. Compiled by `hook/build.dart`
(`_replacementSources`); the directory is also on the include path because it
supplies the `brent.h` that `neldermead.h` includes.

### Process

Written on 2026-10-05 by Vasily (the project's FFI engineer — an AI agent
working for the owner) in a fresh session, under this discipline:

- The bodies of the LGPL code were **not opened**: `src/brent.c`,
  `src/erfunc.c` and the LGPL sections of `src/dist.c` were never read or
  printed. Signatures came from the headers (`brent.h`, `erfunc.h`, `dist.h`,
  `optimc.h`) and from call sites in BSD code (`optimc.c`, `pdist.c`). Section
  boundaries in `dist.c` were found by grepping for SPDX markers and function
  headers only.
- Disclosure, for completeness: two greps over the tree for caller names
  printed a few single lines from inside those bodies — two comment lines of
  `dist.c`'s `betainv` ("initial approximation", "Modified Newton-Raphson
  method"), two lines of `betainv` calling `r8_max`, and one line of
  `erfunc.c`'s `erfcinv` (`erfinv( 1.0 - x )`). Nothing from them was used; the
  replacements follow the published sources below. (The identity
  erfcinv(c) = erfinv(1 - c) is textbook; `erfinv.c` uses it only where
  1 - c is exact.)
- The files were then deleted with `git rm` (and the `dist.c` sections cut), so
  the published package contains no LGPL code.

### Units

| Removed upstream code | License | Replacement | Method source |
| --- | --- | --- | --- |
| `brent_local_min` (`brent.c`) | LGPL-3.0+ | `brent_local_min.c` | R. P. Brent, *Algorithms for Minimization without Derivatives*, Prentice-Hall 1973, ch. 5 (golden section + parabolic interpolation) |
| `brent_zero` (`brent.c`) | LGPL-3.0+ | **dropped** — no caller in ctsa or the shim | — |
| `brent.h` | LGPL-3.0+ | own `brent.h` (declares `brent_local_min`, includes `secant.h`) | — |
| `erf`, `erfc` (`erfunc.c`) | LGPL-3.0+ | **none** — callers now get C99 `<math.h>` from libm; upstream's same-named definitions had shadowed libm | ISO C99 7.12.8 |
| `erf__`, `erfc1__`, `erfcx` (`erfunc.c`) | LGPL-3.0+ | **dropped** — used only inside `erfunc.c` (prototypes in `erfunc.h` remain, unused) | — |
| `erfinv`, `erfcinv` (`erfunc.c`) | LGPL-3.0+ | `erfinv.c` | P. J. Acklam's rational approximation of the normal quantile (2003) as start; Halley iteration on libm `erf`/`erfc` |
| `r8_max` (`dist.c`) | LGPL-3.0+ | `betainv.c` | trivial (larger of two doubles) |
| `betainv` (`dist.c`, AS 109) | LGPL-3.0+ | `betainv.c` | start from Abramowitz & Stegun 26.5.22 or the tail expansions of I_x(a,b); bracketed Halley/bisection on ctsa's CC0 `ibeta`; bounded ulp polish |
| `minfit` (`lls.c`) | MIT, no holder | `minfit.c` | wrapper over ctsa's BSD `svd()` |
| `svd_gr`, `svd_gr2` (`lls.c`) | MIT, no holder | **dropped** — no callers | — |
| `polyroot` (`polyroot.c`, ACM Alg. 419 CPOLY) | ACM non-commercial | `polyroot.c` | Ehrlich (1967) / Aberth (1973) simultaneous iteration with Bini's (1996) Newton-polygon start and backward-error stop; closed form for degree <= 2 |
| `cpoly`, `cpolyroot` (`polyroot.c`) | ACM non-commercial | **dropped** — called only from inside `polyroot.c` (prototypes in `polyroot.h` remain, unused) | — |
| `ppsum` (`talg.c`) | GPL (R tseries, © Trapletti) | `ppsum.c` | Newey & West (1987); Kwiatkowski et al. (1992) eq. (10) |
| `archeck`, `invertroot` (`talg.c`) | C versions of R `arCheck`/`maInvert` (GPL-2+) | `arma_roots.c` | stationarity / invertibility by the roots of φ(z), θ(z) and reflection of MA roots in the unit circle: Box, Jenkins, Reinsel & Ljung (2015) sec. 3.2-3.4; Brockwell & Davis (1991) sec. 3.1, 4.4 |

Call paths: `erfinv` <- `normalinv` <- `tinv_appx`; `erfcinv` <- `gammainv`
<- `chiinv` <- `regress()`; `betainv` <- `tinv`/`finv`, `tinv` <- `regress()`.
`regress()` runs on every SARIMAX fit (OLS starting values; until 0.8.0 also
inside the removed auto-ARIMA unit-root/seasonality tests). `brent_local_min`
<- `fminbnd` (`optimc.c`), whose only caller, Box-Cox automatic lambda
(`boxcox()`), was deleted with `boxcox.c` (Local modifications 11); it stays
because `optimc.c`'s `fminbnd` references it at link time (MSVC links even
unreachable references). Since 0.8.0 `ppsum` has
no caller at all (its only user was `unitroot.c`); the replacement is kept,
unreachable, and dropped by the linker.

### Accuracy

`tool/cleanroom_reference.py` (scipy 1.18.0, numpy 2.5.1) writes
`test/fixtures/cleanroom_reference.txt`; `tool/cleanroom_accuracy.c` checks the
replacements against it and exits non-zero if a gate fails. Result on
2026-10-05 (host clang 22.1.7, Windows x64):

| Check | n | Worst | Gate |
| --- | --- | --- | --- |
| `erfinv` vs `scipy.special.erfinv`, x from 1e-300 to 1-2^-53, both signs | 111 | 3.5e-16 rel. | 1e-14 |
| `erfcinv` vs `scipy.special.erfcinv`, c from 1e-300 to 2-1e-10 | 66 | 4.7e-16 rel. | 1e-14 |
| `betainv` vs `scipy.special.betaincinv`, alpha 1e-12 ... 1-1e-6, 36 shapes, a,b in [0.5, 500] | 612 | 1.7e-12 rel. | 1e-10 |
| `betainv` backward error against ctsa's own `ibeta` | 612 | 3.3e-15 (rel. in x) | 1e-13 |
| `tinv` (via `betainv`) vs `scipy.stats.t.ppf`, df 1 ... 1000 | 196 | 4.2e-13 rel. | 1e-10 |
| `minfit` via `lls_svd2` vs `numpy.linalg.lstsq` | 50 | 2.5 cond(A) eps | 100 cond(A) eps |
| `brent_local_min` vs analytic argmin, 9 functions, `fminbnd`'s tolerances | 9 | 5.4e-9 | 1e-7 |
| domain ends, NaN, invalid shapes | 10 | exact | exact |

The `betainv`-vs-scipy figure is dominated by ctsa's own `ibeta` (up to 6e-9
off scipy's `betainc` in extreme corners such as x -> 1 with a = 500); the
inversion itself is at rounding level (backward error 3e-15). `erfinv` and
`erfcinv` reach full double precision because each Halley step is taken
against libm's `erf`/`erfc` on whichever side avoids cancellation.
`brent_local_min` finds a smooth minimum to ~sqrt(eps), the limit for any
derivative-free minimiser.

Regression check: all 164 package tests — including the golden fits against R
`stats::arima` and statsmodels SARIMAX, and the native leak test — pass with
unchanged tolerances (no test file was modified).

Upstream behaviour noticed while testing (unchanged): `tinv()` returns the
two-sided magnitude |t| for df <= 1000 (its `sign` is computed but never
applied; `regress()` relies on that), but a signed value via `tinv_appx()` for
df > 1000.

### Second round (2026-10-05): `polyroot`, `ppsum`

Written by Vasily in a fresh session, same discipline:

- `polyroot.c` was not read beyond its first 29 lines. Lines 1-10 are the
  header (the ACM notice); lines 12-29, printed by that read, are the bodies
  of two trivial helpers, `mcon()` (machine constants: `macheps()`, `DBL_MAX`,
  `DBL_MIN`, `FLT_RADIX`) and the first half of `CMOD()` (a scaled complex
  modulus). Nothing from them was used: the replacement takes moduli with libm
  `hypot` and needs no machine constant beyond `DBL_EPSILON`.
- The contract of `polyroot()` was established from `polyroot.h`, from its
  BSD-3-Clause callers (`archeck()`/`invertroot()` in `talg.c`, `myarima()`
  in `ctsa.c`) and by **black-box tests of the upstream binary** (coefficient
  order, return codes, what is written for zero leading/trailing
  coefficients, NaN/Inf, degree 0). `ppsum()`'s contract likewise, from
  `talg.h` and black-box tests (it *adds* to `*sum`).
- Disclosure: greps for licence markers and caller names printed a few single
  lines from inside restricted bodies — the AS citation comments at the top of
  `starma()` and `forkal()` (`emle.c`), and call lines inside `unitroot.c`
  (`ppsum(res,N,l,&s2);`, `cumsum(...)`, `interpolate_linear(...)`). Nothing
  from them was used.
- Disclosure, relevant for future work: in the same session the author read
  the BSD-tagged `archeck()`/`invertroot()` of `talg.c` (they call
  `polyroot()`) and R's `arCheck()`/`maInvert()` (to check their origin, see
  "Open licensing questions"). The author is therefore **not** clean-room for
  a rewrite of those two functions.
- `polyroot.c` deleted (`git rm`); `ppsum` cut from `talg.c` (entry 6 above).

**Accuracy and behaviour.** `tool/polyroot_reference.py` (numpy 2.5.1;
references Newton-polished in 50-digit decimal arithmetic on the exact double
coefficients) writes `test/fixtures/polyroot_reference.txt`;
`tool/polyroot_accuracy.c` gates the replacement. 207 polynomials — random
stationary and non-stationary AR(3..12), roots within 1e-2 ... 1e-7 of the
unit circle on both sides, seasonal products (1 - φz)(1 - Φz^s) for
s = 4 ... 288, multiple and clustered roots, random coefficients up to degree
89, moduli spanning 1e-3 ... 1e3 — 3962 roots, plus 37 `invertroot()` cases.
Result on 2026-10-05 (host clang 22.1.7, Windows x64), with the upstream code
run through the same harness for comparison (compiled from git, not read):

| Check | Replacement | Upstream CPOLY |
| --- | --- | --- |
| backward error / ((n+1) eps), worst | 4.6 | 1.0e11 (s = 96, Φ = 0.99) |
| forward error / (cond · eps · n), worst | 5.2 | 1.0e11 |
| `archeck()` verdict equals numpy's | 171/171 | 166/166 (dies before the s = 288 cases) |
| `invertroot()` coefficients vs numpy reflection | 5.6e-14 | 1.3e-14 |
| degree 120 ... 289 | OK | **process dies** (exit 127 / SIGSEGV) |
| documented degenerate inputs | OK | outputs left uninitialised |

Per category (backward / forward, same units, replacement vs upstream): random
AR 3.7/4.1 vs 44/49; near the unit circle 0.6/0.8 vs 12/16; seasonal
s <= 24 4.6/4.7 vs 117/122; s = 48, 96 3.6/3.6 vs 1e11/1e11; multiple roots
3.5/2.6 vs 14/11.

Upstream behaviour found by the black-box tests (not by reading the code): it
kills the process for degree >= ~120 (reachable from auto-ARIMA's `myarima()`
with seasonal periods of that size), and for a zero leading coefficient or a
non-finite coefficient it returns 1 without writing the outputs — `archeck()`
ignores the status and then read uninitialised memory. The replacement
defines those outputs (see the contract in `ctsa_replacements/polyroot.c`);
the callers' decisions on valid input are unchanged.

`ppsum` against the upstream binary (compiled from git, not read): 20 000
random cases (n 1 ... 400, l -2 ... 57, magnitudes 1e-3 ... 1e3) —
**bit-identical**.

All 164 package tests pass with unchanged tolerances, including
`arima_fit_status_test.dart` (codes 10/12 via `archeck()`) and the golden fits.

### Third round (2026-10-05): `archeck`, `invertroot`

The author of round 2 had read both upstream bodies and R's versions, so this
round was done by a **fresh** Vasily session that had read neither, working
from a written specification (Nexus' brief: textbook stationarity and MA-root
reflection rules, Box-Jenkins / Brockwell-Davis) plus black-box tests.

- **Not read:** the bodies of `archeck()`/`invertroot()` in `talg.c`; R's
  `arima.R`/`arima.c` (neither locally nor online); the bodies in
  `unitroot.c`, `stl.c`, `supsmu` in `talg.c`, `starma`/`karma`/`forkal` in
  `emle.c`. From `talg.c` only these were displayed: lines 1-10 (file
  header), the patch comments at 455-466 and 679-683, the two signature lines
  and the two closing braces of the cut range, and the function-header grep.
- **Read (allowed):** `talg.h` (signatures), the BSD callers in `emle.c`
  (`checkroots`, `checkroots_cerr`, and the `invertroot(q, tf + p)` call
  sites with a few lines of context — none inside AS 154/182 code), our
  `polyroot.c`.
- **Contract established by black-box tests** of the upstream object
  (compiled from the tree at `7a5bf14` with the functions renamed, linked
  with our `polyroot`, called as external code): AR convention
  φ(z) = 1 − Σ ar_k z^k and MA convention θ(z) = 1 + Σ ma_k z^k (decided with
  AR(2) (0.7, 0.4) vs (−0.7, −0.4) and MA(2) (2.5, 1) → (1, 0.25)); `archeck`
  returns 1 = stationary, 0 = not; no tolerance band — a root of modulus
  exactly 1 counts as stationary, 1 + 1e-15 flips the verdict;
  `invertroot` returns 1 if it rewrote the coefficients and 0 otherwise,
  leaves trailing zero coefficients in place and writes nothing beyond `q`;
  NaN/inf input: `archeck` 1 (an accident of comparing NaN moduli),
  `invertroot` unchanged with status 0; `p, q <= 0`: 1 / 0.
- **Disclosure — what the author saw about the upstream functions without
  reading them** (none of it was needed; all of it is consistent with the
  spec): this file's "Open licensing questions" item 3 (that upstream has a
  `q0 == 1` special case and rebuilds the polynomial root by root — read
  *after* `arma_roots.c` was written); the comment in our `polyroot.c` that
  upstream `archeck()` ignored polyroot's status and read all `DEGREE` roots;
  the expectations encoded in our `tool/polyroot_accuracy.c` (`archeck`
  verdict = smallest root modulus ≥ 1; `invertroot` status = some root
  inside); the agent memory note that the two functions mirror R's.
- **Deliberate differences** from upstream, both on inputs where upstream's
  result was an accident:
  1. A non-finite coefficient, or a `polyroot()` failure on finite input
     (e.g. a leading coefficient of 1e-300, whose root overflows), makes
     `archeck()` return 0 (non-stationary → fit status 10/12) instead of 1.
     `invertroot()` leaves such input unchanged, as upstream did.
  2. A root r inside the circle is replaced by 1/conj(r). Black-box tests show
     upstream's result equals a reflection by 1/r: identical whenever both
     roots of a conjugate pair are reflected, but when rounding puts only one
     of them inside (roots on the unit circle — e.g. θ(z) = 1 + z⁴, a seasonal
     MA coefficient of exactly ±1) upstream returned a different polynomial
     (1 + z⁴ → coefficients ≈ (0, −1, −1.414, 0)). The replacement returns
     ≈ θ unchanged there.
  Also: a rebuilt coefficient that is not finite leaves the input unchanged
  (status 0) — upstream returned NaN coefficients for 1 + 1e300 (z + z² + z³).

**Differential test.** `tool/arma_roots_accuracy.c` generates 20 000
polynomials from a fixed seed with IEEE-exact operations only — random
coefficients (degree 1-24), products of known roots with moduli 0.3-3, roots
at 1 ± 10^-k (k = 1-15) and exactly 1, double to quadruple roots, trailing
zeros (1-6), sparse seasonal products (1 + Φ₁z^s + Φ₂z^2s)(1 + φ₁z + φ₂z²) for
s = 2-24, and NaN / inf / ±1e±150 / 1e-300 / degree-0 inputs — and feeds each
to both functions (AR with ar = −c, MA with ma = c). Record mode
(`-DARMA_ROOTS_RECORD`, upstream linked in) compared full outputs and wrote
`test/fixtures/arma_roots_reference.txt` (every verdict and status, the
coefficients of every 20th case); default mode replays it. Result
(host clang 22.1.7, Windows x64):

| Check | Result |
| --- | --- |
| `archeck` verdict, finite input, root > 1e-6 from \|z\| = 1 (1e-3 for multiple roots) | 15 470 / 15 470 identical |
| `archeck` verdict, within that band | 4 254 / 4 254 identical (allowed to differ) |
| `archeck` on NaN/inf (250) and polyroot failure (26) | 0 vs upstream's 1 — deliberate (1.) |
| `invertroot` status, all finite input | identical |
| `invertroot` coefficients outside the band | 15 384 within 1e-12 of max\|c\|; 86 between 1e-12 and 1.65e-9, each within 1.82x of the change caused by a ±4-ulp perturbation of the input (ill-conditioned clusters, degree up to 24) |
| roots left inside \|z\| < 1 − 1e-9 after `invertroot`, outside the band | 0 |
| vs the reflection of the exact roots (8 080 cases) | never worse than 2x upstream's error |
| in the band: coefficients differing by > 1e-12 | 1 064 / 4 254, up to 3.9 — upstream's single-root-of-a-pair defect (2.); a root left inside after the call: replacement 639, upstream 1 098; closer to the exact reflection: replacement 339 cases, upstream 0 |

`tool/polyroot_accuracy.c` (numpy reflection reference) also passes with the
replacement: 171/171 verdicts, coefficients within 5.6e-14.

**Fit level.** 256 `fitSarima` runs (orders p, q ≤ 3, P, Q ≤ 1, MLE and CSS,
on log AirPassengers s = 12 and a synthetic 400-point series s = 24) with the
old and new library: 231 bit-identical, all 256 with the same status; the 13
differing AirPassengers fits agree to ~1e-10 in log-likelihood. 2 of the 12
differing synthetic fits moved visibly (log-likelihood ±0.09, forecasts up to
0.26 %): both have a seasonal AR of 0.9995 and seasonal MA roots within 0.2 % of
|z| = 1, where the
likelihood is flat and an ulp-level change in `invertroot` steers the optimiser
elsewhere — one ended better, one worse than before. All 164 package tests pass
with unchanged tolerances (the same "non-stationary" messages are printed as
before).

## Auto-ARIMA removal (2026-10-05)

Process, recorded because the replacement (`autoArima`, Dart) must be written
from papers only and the removal must not have leaked the removed code into
it:

- The bodies of the removed code were **not opened**: `autoutils.c`,
  `unitroot.c`, `seastest.c`, `stl.c`, the auto cluster of `ctsa.c` and
  `supsmu` in `talg.c`. Boundaries were found by a grep for function
  headers (printing only the line number and the function name), by a call
  graph built from the *names* of identifiers followed by `(` inside a line
  range (no code shown), and by the numbers-only line map of entry 6. Struct
  declarations in `ctsa.h` were inspected as identifier tokens only.
- **Incidental exposure, recorded:** (1) a `git diff` of `ctsa.c` displayed
  about 18 lines of the **kept** `sarimax_wrapper()` (its drift / order set-up)
  because the diff aligned on them; that function is not part of the removed
  set. (2) The token view of `ctsa.h` showed the field names of the removed
  structs `auto_arima_set`, `myarima_set` and `aa_ret_set` (e.g. `idrift`,
  `ic`, `aicc`). No statement of the removed function bodies was displayed.
  Earlier sessions (the 2026-10-05 licence audit) had inspected parts of the
  auto cluster to establish its origin — see "License audit" — so the author
  of these edits is not a clean-room implementer of `autoArima` either; the
  Dart implementation is to be written from the methodology specification
  (Hyndman & Khandakar 2008; Kwiatkowski et al. 1992; FPP3) only.
- The deleted files were removed with `git rm`; history keeps them, the
  published package does not.

## License audit (2026-10-05)

Engineering assessment, not legal advice. "Reach": (a) executed by every
fixed-order fit/forecast (`fitSarima`/`fitSarimax`/`arimaForecast` — what
Sweet Limit uses); (b) only by `autoArimaForecast` (**removed in 0.8.0** — all
rows marked b are gone from the tree); (c) linked but not
executed by the test suite; (d) dropped by the linker (source only).
Established with gc-sections links per export set (arm64, -O0 and -O3) and a
coverage-instrumented build running the full test suite.

| Unit | Origin / licence | OK for BSD + commercial? | Reach |
| --- | --- | --- | --- |
| `emle.c` `starma`, `karma` (AS 154), `forkal` (AS 182), `inclu2`, `regres` (AS 75) | RSS: distribution "provided that no fee is charged" | **unclear** | a |
| `emle.c` rest (CSS, likelihood wrappers, `checkroots*`) | BSD (Rafat Hussain); `checkroots` prints R `arima()`'s messages | yes / low risk | a |
| `talg.c` `archeck`, `invertroot` | BSD-tagged; line-for-line R `arCheck`/`maInvert` (GPL-2+) | replaced (clean-room round 3) | — |
| `neldermead.c` `nel_min` (AS 47) | RSS "no fee" | unclear | c |
| `ctsa.c` `auto_arima*`, `search_arima`, `myarima`, `newmodel`; `autoutils.c`; `seastest.c` | BSD-tagged; C versions of R package **forecast** (GPL-3) | **no** (likely derivative) — **removed 0.8.0** | b |
| `unitroot.c` | R tseries/urca (GPL), per its header | **no** — **removed 0.8.0** | b |
| `stl.c` | netlib STL (no notice; AT&T, 1990) as modified in R by Ripley/Maechler (GPL-2+) | **no / unclear** — **removed 0.8.0** | b |
| `talg.c` `supsmu`, `smooth_`, `supsmu_` | Friedman 1984 "all rights reserved", via R `ppr.f` | **no** — **removed 0.8.0** | b (c in tests) |
| `ctsa.c` `sarimax_wrapper*`; `boxcox.c` | BSD-tagged; `sarimax_wrapper*` never audited (its drift handling resembles R forecast's `Arima()`) | **removed 2026-10-05** (Local modifications 11) | — |
| `lnsrchmp.c` `cvsrch`, `cstep` | MINPACK (Moré-Thuente); Argonne licence, BSD-style + acknowledgment | yes (notice in THIRD_PARTY_NOTICES) | c |
| `nls.c` | MINPACK | yes (notice in THIRD_PARTY_NOTICES) | d |
| `matrix.c` `svd` | EISPACK translation (netlib, no notice) | yes | a (SARIMAX) |
| `matrix.c` `svd_sort` | csa, Pavel Sakov / CSIRO, BSD-style | yes (notice in THIRD_PARTY_NOTICES) | a |
| `dist.c` CC0 parts (`gamma_log`, `pgamma`, `ibeta`, ...) | W. J. Cody's specfun (netlib) / Abramowitz & Stegun formulas | yes, low risk | a |
| `dist.c` `log1p` | MIT, Karney | yes (notice in THIRD_PARTY_NOTICES) | a |
| `initest.c` Burg (whole file) | "C translation of Cedrick Collomb's C++" — no licence found | **removed 2026-10-06** (Local modifications 17) | — |
| `wavefilt.c`, `wavelib.c`, `wtmath.c`, `waveletarima.c` | BSD, Rafat Hussain + Holger Nahrstaedt | yes (notice in THIRD_PARTY_NOTICES); not compiled | — |
| everything else (`boxjenkins`, `conjgrad`, `conv`, `errors`, `filter`, `hsfft`, `lls`, `newtonmin`, `optimc`, `pdist`, `pred`, `real`, `regression`, `secant`, `spectrum`, `stats`, rest of `ctsa.c`/`matrix.c`/`talg.c`) | BSD-3-Clause, Rafat Hussain; no other attribution found | yes | a-d |
| `ctsa_replacements/*`, `native/*` | BSD-3-Clause, Vladislav Zyuzko | yes | a-c |

## Open licensing questions (unresolved; engineering notes, not legal advice)

Resolved on 2026-10-05: ACM Algorithm 419 (`polyroot.c` replaced), `ppsum`
(replaced), the unreachable Numerical Recipes, AS 197 and R-derived
transform code (removed), and `talg.c` `archeck()`/`invertroot()` — C
versions of R's `arCheck()`/`maInvert()` (GPL-2+), on every fit — replaced
clean-room (round 3). Still from R, and left as is: the two short error
messages that `emle.c` `checkroots*()` (BSD, Rafat Hussain) print, which match
R `arima()`'s ("non-stationary AR part", "non-stationary seasonal AR part") —
plain factual strings. What remains:

1. **Applied Statistics algorithms (Royal Statistical Society)** — on the
   path of **every** fit and forecast. `emle.c`: `starma`/`karma` (AS 154 —
   exact ARMA likelihood by Kalman filter), `forkal` (AS 182 — forecasts) and
   `inclu2`/`regres` (AS 75, Gentleman 1974, used by `starma`; not marked
   upstream — identified from `forkal`'s comment "REGRES from AS 75");
   `neldermead.c`: `nel_min` (AS 47 — linked, used only for optimiser method
   0, which tseries never selects). StatLib's apstat index: "The Royal
   Statistical Society holds the copyright to these routines, but has given its
   permission for their distribution provided that no fee is charged."
   (https://lib.stat.cmu.edu/apstat/index). Upstream tags these sections CC0.
   For comparison, R's `stats` ships `starma.c` — "Code in this file based on
   Applied Statistics algorithms AS154/182 (C) Royal Statistical Society 1980,
   1982" — under GPL-2+. Open point: whether a paid app that embeds the
   binary "charges a fee" for the routines.
2. ~~**Auto-ARIMA is a C version of GPL code.**~~ **Resolved in 0.8.0 by
   removal** (Local modifications entry 9). For the record: besides the files
   whose comments said so — `unitroot.c` (R tseries/urca), `stl.c` (R's
   modified STL), `talg.c` `supsmu` (Friedman, "all rights reserved") —
   inspection had shown that ctsa's order search (`ctsa.c`: `auto_arima1`,
   `search_arima`, `myarima`, `newmodel`), `autoutils.c` (`ndiffs`,
   `nsdiffs`, `is_constant`) and `seastest.c` (`OCSBtest`,
   `calcOCSBCritVal`, `SHtest`, `mstl`) follow the R package **forecast**
   (GPL-3): same function names, the OCSB critical-value formula with
   identical constants, `myarima`'s `minroot < 1.01` test, and a verbatim
   warning string (https://github.com/robjhyndman/forecast; GPL-3). Order
   selection is now `autoArima`, written in Dart from the published
   methodology.
3. ~~**`ctsa.c` `sarimax_wrapper*` — not audited.**~~ **Resolved 2026-10-05
   by removal** (Local modifications 11), together with `boxcox.c`, its last
   user. It had been unreachable since 0.8.0 (its only caller was the removed
   `arima2`); a diff had incidentally shown its drift set-up ("no drift when
   d + D > 1"), the same rule as R forecast's `Arima()` — a common modelling
   rule, not evidence of copying, but it was never audited and tseries never
   needed it.
4. ~~`initest.c` (Burg, "C translation of Cedrick Collomb's C++
   Implementation") — no licence found for the original; unreachable, source
   only.~~ **Resolved 2026-10-06 by removal** (Local modifications 17):
   `initest.c`/`initest.h` deleted together with their unreachable callers
   in `ctsa.c`.
5. Not problems, recorded for completeness: `matrix.c` `svd_sort` is from
   Pavel Sakov's `csa` (Copyright 2000-2008 Pavel Sakov and CSIRO,
   BSD-style; notice in `THIRD_PARTY_NOTICES`); `lnsrchmp.c` and `nls.c` are MINPACK
   translations (Argonne licence, BSD-style with an acknowledgment clause;
   notice added to `THIRD_PARTY_NOTICES` on 2026-10-05); the wavelet files carry a second
   BSD holder, Holger Nahrstaedt (notice added to `THIRD_PARTY_NOTICES`).
