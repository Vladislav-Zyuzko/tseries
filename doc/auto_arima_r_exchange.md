# autoArima synthetic study — exchange with the external R run

Gate C of the autoArima specification (§9.C, §13, §14) compares our
decisions with R (`forecast`, `urca`) on the **same** synthetic series. R is
run as an external tool by someone other than the implementer; only its
numeric output comes back. This file fixes the formats in both directions.

All files are comma-separated, UTF-8, with one header line. Numbers are
written with 17 significant digits (round-trip exact). `t` is 1-based.

## 1. Export (ours → R)

Written by `tool/auto_arima_synthetic_study.dart --export <dir>` (the heavy
test `test/auto_arima_synthetic_test.dart` with `AUTO_ARIMA_C_EXPORT=<dir>`).
The study may be split over several runs (`--models`), each with its own
directory; the comparison tool accepts several `--export` directories.

| file | columns | content |
|---|---|---|
| `index.csv` | `series_id,model,n,rep,seed,period,d_star,D_star,d_true,D_true,holdout` | one row per series |
| `series.csv` | `series_id,model,n,rep,seed,period,d_star,D_star,t,value` | **training** part, long format — what R must be run on |
| `holdout.csv` | `series_id,t,value` | the held-out tail (h = 10, or 2s for a seasonal series); R does not need it |
| `ours.csv` | `series_id,d_star,D_star,step_p,step_q,step_P,step_Q,step_c,step_aicc,ex_p,…,ex_aicc,mase_step,no_model` | our stepwise and exhaustive choices |
| `ours_kpss.csv` | `series_id,step,T,lag,eta,rejects` | our KPSS steps (default lag rule `short`) |

`series_id` is `<model>_n<n>_r<rep>`; `model` ∈ `wn, ar1, ma1, arma11, ar2,
ima, ari, arima111, imadrift, airline, ar098`; `period` is 12 for `airline`,
1 otherwise. `d_star`/`D_star` are **our** differencing orders.

## 2. R results (R → us)

Three files in one directory, plus `PROVENANCE.txt` (R, forecast and urca
versions, the exact commands, the date). Per series x (a `ts` with
frequency = `period`, training part only):

| file | columns | how |
|---|---|---|
| `r_orders.csv` | `series_id,mode,p,d,q,P,D,Q,has_mean,has_drift,aicc_R` | `mode = step`: `auto.arima(x, d = d_star, D = D_star, …, stepwise = TRUE, approximation = FALSE)`; `mode = ex`: the same with `stepwise = FALSE, max.order = 5`; `mode = auto`: `auto.arima(x, test = "kpss", stepwise = TRUE, approximation = FALSE, …)` with R's own d, D. `has_mean`/`has_drift` are 0/1 (or TRUE/FALSE). `aicc_R` is reference only. The full argument list is spec §13. |
| `r_kpss.csv` | `series_id,step,T,lag,eta_R,ndiffs_R` | z = x after `D_star` seasonal differences; for step k = 1, 2, … on z differenced k − 1 times: `ur.kpss(z_k, type = "mu", use.lag = lag)@teststat` with **our** lag (`ours_kpss.csv`, column `lag`; rule floor(3·√T/13)) — one row per step that we ran; `ndiffs_R = ndiffs(z, test = "kpss", type = "level", alpha = 0.05)` repeated on every row. |
| `r_seasonal.csv` | `series_id,F_S_STL,nsdiffs_R` | period > 1 only: F_S of `mstl(x)` (max(0, 1 − var(remainder)/var(seasonal + remainder))) and `nsdiffs(x)` with its defaults. |

## 3. Comparison

```
dart run tool/auto_arima_r_compare.dart --export <dir> [--export <dir> …] \
    --r <R results dir> --out report.md
```

Per (model, n) cell it reports:

* **C-d′** — η compared / violations (|η − η_R| > 1e-6·max(1, η) at equal T
  and lag), lag/T mismatches, d̂ = ndiffs_R, and mismatches split into
  "at a tie" (some step with |η − 0.463| ≤ 1e-6) and defects. Any violation
  or defect → `DEFECTS: FOUND`.
* **C-S′** — share (Wilson 95 %) of series where AICc(our stepwise) >
  AICc_ours(R's stepwise order at our d*, D*) + 1e-3; R's order is refitted by
  our fit path, so both AICc are ours (spec §9.B convention). Gate ≤ 5 %,
  upper bound ≤ 8 %. Each case is listed under `details`.
* **MASE against the reference** — median (bootstrap 95 % CI) of
  MASE(our stepwise) / MASE(R `auto` order, its own d and D, fitted by our
  `fitSarimax`) on the held-out tail. Gate ≤ 1.05.
* R's own step = ex share (from `aicc_R`), and D* = nsdiffs_R.

`--selftest <dir>` writes R-format files built from our own results and
compares against them: every C-d′ count must be clean (plumbing check).

## 4. Gate D′ — rolling-origin windows of the reference series

Written by `tool/auto_arima_rolling_windows.dart --export <dir>` (Windows of
spec §9.D/§15.4: the six reference series, the last 10 origins that leave a
full horizon; h = s for a seasonal series, 10 otherwise; orders re-chosen at
every origin).

Same files as §1 (`index.csv`, `series.csv`, `holdout.csv`, `ours.csv`,
`ours_kpss.csv`) with `series_id = <series>_o<k>`, k = 1..10 (k = 10 is the
longest window), `model` = the series name, `n` = the window length,
`rep` = k, `seed` = 0, empty `d_true`/`D_true`; plus `ours_windows.csv`:

| column | meaning |
|---|---|
| `series_id, series, k, window_n` | the window |
| `cycles, D_source, F_S, threshold` | our seasonal-difference decision (§15.1) |
| `d_star, D_star` | our differencing on the window |
| `mase_step, mase_ex` | MASE of our stepwise / exhaustive choice on the holdout |
| `mase_step_d1` | our stepwise with `d: 1` forced (price of d = 0, wwwusage) |
| `mase_step_sD1` | our stepwise with `seasonalD: 1` (seasonal series) |
| `mase_oldR_step, mase_oldR_ex` | old gate D: the full-series R order fitted at our (d*, D*) of the window |

**What R returns for the windows** — the same three files as §2, with one
difference in `r_orders.csv`: on windows R chooses its **own** d and D
(D′ compares procedures, spec §15.2):

* `mode = step`: `auto.arima(x, stepwise = TRUE, approximation = FALSE,
  test = "kpss", max.p = 5, max.q = 5, max.P = 2, max.Q = 2, ic = "aicc",
  seasonal = (s > 1))`;
* `mode = ex`: the same with `stepwise = FALSE, max.order = 5`;
* `mode = auto` may be omitted (it equals `step` here);
* `r_kpss.csv`: `ur.kpss(type = "mu", use.lag = our lag)` on our steps and
  `ndiffs(z, test = "kpss", type = "level", alpha = 0.05)` on our base z
  (C-d′ on the windows);
* `r_seasonal.csv`: optional for the windows.

**Comparison:**

```
dart run tool/auto_arima_r_compare.dart --windows --export <cexp_D> \
    --r <R results for the windows> --out d_prime.md
```

`--windows` groups the C-d′ table per series and appends gate D′: a
per-window table (window length, cycles, our D source, (d, D) ours and R's
step/ex, MASE(ours)/MASE(o_R) for step and ex, the old-gate-D ratio,
MASE(ours)/MASE(ours with d = 1), MASE(ours with seasonalD: 1)/MASE(o_R
step)) and per series the medians with the verdict "≤ 1.05" (one-sided;
values below 0.95 are not failures, spec §15.2).
