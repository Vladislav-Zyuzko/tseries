# Platform verification

How the cross-platform behaviour promised in the README ("Determinism",
under *Automatic order selection*) was checked, on which binaries, and with
what result. Measured on 2026-10-06 with tseries 0.8.0, Dart 3.10.9 /
Flutter 3.38.10.

## Build setting under test

The native core (ctsa + shim + clean-room replacements) is compiled **without
floating-point contraction** (`hook/build.dart`, `_floatingPointFlagsFor`):

| Toolchain | Flag | Default without it |
| --- | --- | --- |
| clang (Android NDK 29, Linux, macOS/iOS) | `-ffp-contract=off` | `on` — fuses `a * b + c` within an expression wherever the target has FMA (every arm64) |
| GCC (accepted for completeness; not the Linux default of `native_toolchain_c`) | `-ffp-contract=off` | `fast` — fuses across expressions too |
| MSVC (Windows) | `/fp:precise` (explicit) | `/fp:precise`, which does not contract since Visual Studio 2022; `/fp:contract` is never passed |

`-ffast-math`, `-Ofast`, `-march=native`/`-mcpu=native` and `/fp:fast` are
not used anywhere in the build. Why contraction matters: a fused multiply-add
rounds once where `a * b + c` rounds twice; the optimiser amplifies the
difference, and near the root-admissibility boundary that is enough to change
which candidate autoArima accepts (see E2/E3 below).

## Gates

* **E0 (static).** `tool/fma_count.dart` counts fused multiply-add
  instructions in a built binary: aarch64 `fmadd/fmsub/fnmadd/fnmsub` and
  vector `fmla/fmls` (plus the SVE forms), x86 `vfmadd*/vfmsub*/vfnmadd*/
  vfnmsub*` (FMA3/FMA4), 32-bit arm `vfma/vfms/vfnma/vfnms` (`vmla/vmls` are
  not fused and not counted). Expected: **0** on every ABI.
* **E2 (numbers).** On fits accepted on both platforms: coefficients within
  1e-3·max(1,|c|), log-likelihood within 1e-4.
* **E3 (decisions).** The same (d, D) and KPSS statistic, the same search
  path and the same selected model.

The E2/E3 set is 166 autoArima runs: set B (6 reference series × stepwise and
exhaustive = 12 runs), the boundary set F (44 runs: short seasonal series at
the length limits of the search) and subset C of the synthetic benchmark (110 runs),
4 063 candidate fits in total. Windows is the reference; every platform runs
the same driver and series (Android: a throwaway probe app, results streamed
through logcat; Linux: the host driver in Docker).

## Results

### E0 — fused instructions

| Binary | Size | Instructions | Fused |
| --- | ---: | ---: | ---: |
| Android arm64-v8a `libctsa.so` (stripped, from the release APK) | 184 768 B | 41 329 | **0** |
| Android x86_64 `libctsa.so` (stripped) | 240 392 B | 52 689 | **0** |
| Android armeabi-v7a `libctsa.so` (stripped) | 137 672 B | 32 783 | **0** |
| Linux x64 `libctsa.so` (Docker `dart:3.10.9`, clang 19.1.7, unstripped hook output) | 268 192 B | 51 956 | **0** |
| Windows x64: the 32 object files of the hook's sources (MSVC 19.44, `/O2 /fp:precise`) | — | 117 280 | **0** |
| *Control:* Android arm64-v8a built **without** the flag | 183 712 B | 41 065 | 497 (`fmadd` 325, `fmla` 116, `fmsub` 32, `fnmsub` 14, `fmls` 10) |

`ctsa.dll` as a whole (468 480 B) shows 163 fused instructions: the DLL links
the C runtime statically, and the UCRT's math functions ship FMA3 variants
(`exp_fma`, `log_fma`, `erfc_fma`, … in `libucrt.lib`) that the runtime
selects by CPUID. They are library code, like bionic's or glibc's `libm` on
the other platforms (which are shared libraries and so are not inside our
binary). Our own code has none; that is what E0 checks on Windows
(`--hook-objects`, below). A consequence worth knowing: on Windows the
elementary functions may take a different code path on a CPU without FMA3
(pre-2013); that was not measured.

Size cost of the flag: arm64 +1 056 B (+0.6 %); x86_64 +608 B (+0.3 %) —
x86_64 had no FMA to begin with (the NDK's x86_64 baseline has none), but
LLVM's `contract` flag also steers instruction selection and vectorisation,
so the code changes slightly. The numbers do not: the x86_64 probe run with
the flag produced a journal and fit table **byte-identical** to the run
without it (166/166 runs, 4 063/4 063 fits).

### E2/E3 — against Windows

| Platform | (d, D), KPSS η | Same model | Same search path | Verdict flips (not selected) | Accepted on both: outside E2 | Bit-identical fits vs Windows |
| --- | --- | --- | --- | --- | --- | --- |
| Android x86_64 (emulator, native) | identical, η bit-identical | 166/166 | 166/166 | 8 | 0 of 3 343 | 2 683 of 4 063 |
| Android arm64-v8a (emulator, ARM translation), with the flag | identical, η bit-identical | 166/166 | 166/166 | 8 (the same 8) | 0 of 3 343 | 2 696 of 4 063 |
| Linux x64 (Docker, glibc, clang 19), with the flag | identical, η bit-identical | 166/166 | 166/166 | 8 (the same 8) | 0 of 3 343 | 2 696 of 4 063 |
| *Before:* Android arm64-v8a without the flag (FMA) | identical | 165/166 | 165/166 | 15 fits | 7 fits with \|Δll\| > 1e-4 and 12 coefficients outside, of 3 328 | — |

Android arm64 (bionic) with the flag and Linux x64 (glibc) give
**byte-identical** journals and fit tables; Android x86_64 differs from them
in 99 of 4 063 fits, with no verdict flip and no E2 violation between them.

The arm64 binary verified here is the hook's own build: the hook-built
`libctsa.so` is byte-identical (SHA-256 `12b31e75…c4caee2`) to the one used for
the 166-run arm64 experiment, so that run counts for the hook build.

### The 8 verdict flips

All 8 are on candidates that were not selected on any platform, and the
selected model and its criterion are the same everywhere. They are the same
8 on Android x86_64, Android arm64 and Linux; all are Windows (UCRT) versus the
other two math libraries, so they are the residual libm effect, not FMA.
Six distinct fits (two appear in both the stepwise and the exhaustive search).

| Run | Candidate | Windows | Android / Linux | Class |
| --- | --- | --- | --- | --- |
| C `wn_n200_r3` exhaustive #26 | ARIMA(2,0,2) | rejected, fit status 10 (AR root 0.99985) | rejected, roots (AR 1.0000001, MA 1.004) | boundary status |
| C `arma11_n200_r0` exhaustive #28 | ARIMA(2,0,3) | rejected, roots (AR 1.00001) | rejected, fit status 10 (AR root 0.99985) | boundary status |
| C `ar2_n200_r3` exhaustive #26 | ARIMA(2,0,2) | rejected, roots (AR 1.0015) | rejected, fit status 10 (AR root 0.99998) | boundary status |
| C `ari_n200_r0` stepwise #15 and exhaustive #24 | ARIMA(2,1,1) | accepted, ll −302.028, MA root 1.63, AICc 612.26 | rejected, MA root 1.0000003 (AR 0.994), ll −301.846 | boundary MA optimum |
| C `imadrift_n200_r0` exhaustive #38 | ARIMA(4,1,1) | accepted, ll −276.979, MA root 1.32, AICc 566.39 | rejected, MA root 1.002 (AR 1.00007), ll −276.721 | boundary MA optimum |
| C `imadrift_n200_r1` stepwise #9 and exhaustive #17 | ARIMA(1,1,2)+c | accepted, ll −293.448, MA root 1.14, AICc 597.21 | rejected, MA root 1.0000002, ll −293.263 | boundary MA optimum |

Classes:

* **Root at 1.01** (a root modulus straddling the `rootMargin` of 1.01): **0**
  of 8.
* **Boundary status** (3): the optimiser ends at an AR root on the unit
  circle; on one platform ctsa reports non-convergence (status 10), on the
  other it converges there and the root check rejects it. Rejected on both.
* **Boundary MA optimum** (5 flips, 3 fits): Windows stops at an interior
  local optimum (accepted), the other platforms reach the higher likelihood on
  the unit circle (|MA root| = 1), which the root check correctly rejects —
  the case described under "root check vs. boundary optima" in the README.
  Where accepted, the candidate's AICc is 3.9–8.2 above the selected model's,
  so it was never close to being chosen.

### Run time

| Platform | Total for the 166 runs |
| --- | ---: |
| Windows x64 (JIT, host) | 65.1 s |
| Linux x64 (Docker, JIT) | 27.0 s |
| Android x86_64 (emulator, AOT, native) | 36.7 s without the flag, 29.8 s with it (single runs) |
| Android arm64 (emulator, AOT, ARM translation) | 601.7 s with FMA, 281.0 s without (single runs) |

The cost of the flag on arm64 is **not measurable on the emulator**: the
code runs under binary translation, the two runs differ by a factor of two in
the "wrong" direction, which says more about the translator's state than
about FMA. A measurement needs a physical arm64 device. Expected order of
magnitude: small — the fused instructions were 1.2 % of the arm64 code.

## Not verified

* iOS and macOS (Apple arm64 is affected by contraction exactly like Android
  arm64; the flag is in place, the binaries have not been built here).
* A physical Android device; armeabi-v7a at run time (built, E0 = 0).
* Windows on a CPU without FMA3 (UCRT math path).

## How to re-run

E0, from `packages/tseries` (the tool is in the repository, not in the
published archive; `llvm-objdump` is taken from `--objdump`, `$LLVM_OBJDUMP`,
the newest Android NDK, or `PATH`):

```
# Android: libraries from a release APK (unzip lib/<abi>/libctsa.so)
dart run tool/fma_count.dart lib/arm64-v8a/libctsa.so lib/x86_64/libctsa.so lib/armeabi-v7a/libctsa.so

# Linux / macOS: the hook output (.dart_tool/lib/ or hooks_runner/shared/tseries/build/<hash>/)
dart run tool/fma_count.dart .dart_tool/lib/libctsa.so

# Windows: our object files, not the DLL (see above)
dart run tool/fma_count.dart --hook-objects .dart_tool/hooks_runner/shared/tseries/build/<hash>
```

Exit code 0 means no fused instruction (`E0 OK`); `--verbose` lists them by
symbol. Do not run it inside the package while `dart test` is running there:
`dart run` re-runs the build hook, and both copy `.dart_tool/lib/ctsa.dll`.

E2/E3 need the 166-run driver and a reference journal from Windows; the
driver and the comparison scripts are development material outside this
package. A binary-level shortcut is valid: when a rebuilt library is
byte-identical to a verified one, the verification carries over.
