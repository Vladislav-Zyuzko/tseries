// SPDX-License-Identifier: BSD-3-Clause
//
// Build hook for tseries's native core.
//
// Compiles the vendored `ctsa` C library (third_party/ctsa) together with our
// thin C shim (native/tseries_ctsa.c) into a single dynamic **code asset**
// named `ctsa` (asset id `package:tseries/ctsa`). The Dart FFI layer looks the
// asset up by that id via `@Native` / `@DefaultAsset`.
//
// Uses `package:native_toolchain_c`'s CBuilder, per the current Dart build
// hooks / code assets flow (`dart.dev/tools/hooks`). No OS-specific build files.

import 'package:code_assets/code_assets.dart';
import 'package:hooks/hooks.dart';
import 'package:native_toolchain_c/native_toolchain_c.dart';

/// Vendored ctsa translation units (see third_party/ctsa/PROVENANCE.md).
///
/// This is the ARIMA/SARIMA path only. The wavelet cluster
/// (wavefilt/wavelib/wtmath/waveletarima) is deliberately excluded: it is a
/// self-contained island — no non-wavelet source references its symbols and no
/// header reachable from `ctsa.h` (what our shim includes) pulls in a wavelet
/// header — so dropping it keeps the link clean while shrinking the binary and
/// narrowing the vendored-license surface. Re-add these four .c files (and the
/// `_USE_MATH_DEFINES` define below) together if a wavelet API is ever exposed.
///
/// The auto-ARIMA cluster (`autoutils.c`, `unitroot.c`, `seastest.c`, `stl.c`,
/// plus its functions in `ctsa.c` and `supsmu` in `talg.c`) is not here because
/// it no longer exists in the tree: it was a C port of GPL-licensed R code
/// (forecast/tseries/urca/stats) and was deleted for licensing reasons in
/// 0.8.0. Order selection is now `autoArima`, in pure Dart.
///
/// `boxcox.c` is gone too: its last caller was ctsa's unaudited
/// `sarimax_wrapper*` group (removed from `ctsa.c`; PROVENANCE.md, Local
/// modifications 11), which tseries never called.
///
/// `initest.c` (Burg / Yule-Walker / Hannan-Rissanen estimators, a translation
/// of C++ code whose licence was never found) is gone as well, together with
/// its unreachable callers in `ctsa.c` (PROVENANCE.md, Local modifications 17).
const _ctsaSources = <String>[
  'third_party/ctsa/src/boxjenkins.c',
  'third_party/ctsa/src/conjgrad.c',
  'third_party/ctsa/src/conv.c',
  'third_party/ctsa/src/ctsa.c',
  'third_party/ctsa/src/dist.c',
  'third_party/ctsa/src/emle.c',
  'third_party/ctsa/src/errors.c',
  'third_party/ctsa/src/filter.c',
  'third_party/ctsa/src/hsfft.c',
  'third_party/ctsa/src/lls.c',
  'third_party/ctsa/src/lnsrchmp.c',
  'third_party/ctsa/src/matrix.c',
  'third_party/ctsa/src/neldermead.c',
  'third_party/ctsa/src/newtonmin.c',
  'third_party/ctsa/src/nls.c',
  'third_party/ctsa/src/optimc.c',
  'third_party/ctsa/src/pdist.c',
  'third_party/ctsa/src/pred.c',
  'third_party/ctsa/src/real.c',
  'third_party/ctsa/src/regression.c',
  'third_party/ctsa/src/secant.c',
  'third_party/ctsa/src/spectrum.c',
  'third_party/ctsa/src/stats.c',
  'third_party/ctsa/src/talg.c',
];

/// tseries' own BSD-3-Clause implementations of ctsa functions whose upstream
/// code was removed from the vendored tree for licensing reasons: the
/// LGPL-3.0-or-later units (`brent.c`, `erfunc.c`, two sections of `dist.c`),
/// `minfit()` from `lls.c` (tagged MIT with no named copyright holder),
/// `polyroot()` (`polyroot.c`, ACM Algorithm 419 under ACM's non-commercial
/// licence), `ppsum()` from `talg.c` (from the GPL R package tseries) and
/// `archeck()`/`invertroot()` from `talg.c` (ports of R's GPL-2+
/// `arCheck()`/`maInvert()`).
/// Same names and signatures as upstream, so the vendored callers are
/// untouched. See third_party/ctsa/PROVENANCE.md, "Clean-room replacements".
///
/// The directory is also on the include path: it supplies the `brent.h` that
/// the vendored `neldermead.h` includes.
const _replacementSources = <String>[
  'third_party/ctsa_replacements/arma_roots.c',
  'third_party/ctsa_replacements/betainv.c',
  'third_party/ctsa_replacements/brent_local_min.c',
  'third_party/ctsa_replacements/erfinv.c',
  'third_party/ctsa_replacements/minfit.c',
  'third_party/ctsa_replacements/polyroot.c',
  'third_party/ctsa_replacements/ppsum.c',
];

/// Libraries the ctsa core must be linked against, for a given target OS.
///
/// ctsa is dense floating-point code: it references `log`, `exp`, `pow`,
/// `log10` and `sin` from the C math library. Where those symbols live is a
/// per-OS decision, and getting it wrong does NOT fail the build — it produces a
/// shared library that links and packages happily and then fails at `dlopen`
/// on-device with `cannot locate symbol "log"`. That was a real, shipped bug
/// (the Windows host build and its tests could never have caught it), so this
/// mapping is explicit rather than left to the toolchain's defaults.
///
/// * **Android / Linux** — libm is a *separate* shared library and is NOT linked
///   by default; it must be requested explicitly, or every math symbol stays
///   undefined at load time. This is the case that broke.
/// * **macOS / iOS** — the math functions live in libSystem, which is always
///   linked; `-lm` resolves to a stub kept for compatibility, so asking for it
///   is harmless and keeps this list honest about the dependency.
/// * **Windows** — math lives in the CRT (linked automatically) and there is no
///   `m.lib`. `libraries` is emitted as `<name>.lib` for MSVC, so naming `m`
///   here would break the host build outright.
List<String> _mathLibrariesFor(OS targetOS) =>
    targetOS == OS.windows ? const [] : const ['m'];

/// Floating-point flags that make the native core compute the same thing on
/// every platform: **no floating-point contraction** (no FMA).
///
/// Why: autoArima compares candidate models by likelihood and by a root
/// admissibility margin. A fused multiply-add rounds once where `a * b + c`
/// rounds twice, so a contracted build drifts from an uncontracted one by a
/// few ulps per operation, and the optimiser amplifies that. On Android arm64
/// clang contracts by default (`-ffp-contract=on`): the build had 371
/// `fmadd`/`fmsub` instructions, only ~10% of fits were bit-identical to
/// x86_64 and 1 of 166 autoArima runs selected a different model (E2/E3,
/// doc/platform_verification.md). With contraction off, arm64 matched x86_64
/// on 166/166 runs. Cross-platform determinism of model selection is part of
/// the package's contract (README, "Determinism"), so this is not optional.
///
/// * **clang (Android, Linux, macOS, iOS)** — default is `on` (contraction
///   within an expression); `-ffp-contract=off` disables it.
/// * **GCC (a Linux host may resolve to it)** — default is `fast` (contraction
///   across expressions too); it accepts the same `-ffp-contract=off`.
/// * **MSVC (Windows)** — `/fp:precise` is the default and, since Visual
///   Studio 2022, does not contract (`/fp:contract` is a separate opt-in that we
///   never pass; older MSVC could contract under `/fp:precise`, but only with
///   `/arch:AVX2` on x64, which we do not pass either). `/fp:precise` is still
///   passed explicitly so that a `/fp:fast` coming from the `CL` environment
///   variable (prepended by cl.exe) cannot silently override it.
///
/// Never add `-ffast-math`, `-Ofast`, `-march=native`/`-mcpu=native`, `-mfma`
/// or `/fp:fast` here: each re-enables contraction or reassociation, and
/// `native` additionally makes the binary depend on the build machine.
/// `tool/fma_count.dart` (gate E0) checks the built binaries for fused
/// instructions; the expected count is 0 on every ABI.
List<String> _floatingPointFlagsFor(OS targetOS) => targetOS == OS.windows
    ? const ['/fp:precise']
    : const ['-ffp-contract=off'];

/// Compiler/linker flags that keep the library's ABI down to the shim and let
/// the linker drop every vendored function nothing reachable calls.
///
/// * `-fvisibility=hidden` — only symbols marked `TSERIES_EXPORT` (the
///   `tseries_*` shim entry points, `visibility("default")`) are exported.
///   Without it all ~450 ctsa functions were exported, including generic names
///   (`mean`, `var`, `filter`, `gamma`, `log1p`, …) that could interpose with,
///   or be interposed by, other libraries in the host process.
/// * `-ffunction-sections -fdata-sections` + linker garbage collection —
///   `--gc-sections` for ELF (Android/Linux), `-dead_strip` for Mach-O
///   (macOS/iOS). Once nothing but the shim is exported, unreachable ctsa code
///   (spectra, unused optimisers, ...) is removed.
/// * **Windows: no flags here** (the `/fp:precise` above is the only one).
///   Exports there are explicit (`__declspec(dllexport)`
///   on the shim) so visibility is already right, `/O2` implies function-level
///   linking (`/Gy`), and the linker's `/OPT:REF` is on by default for non-debug
///   links. The flags here are clang/GCC syntax; passing them on Windows would
///   break `cl.exe`, and `CBuilder` (0.17 to 0.19) offers no slot after `/link`.
List<String> _sizeAndVisibilityFlagsFor(OS targetOS) => switch (targetOS) {
  OS.windows => const [],
  OS.macOS || OS.iOS => const [
    '-fvisibility=hidden',
    '-ffunction-sections',
    '-fdata-sections',
    '-Wl,-dead_strip',
  ],
  _ => const [
    '-fvisibility=hidden',
    '-ffunction-sections',
    '-fdata-sections',
    '-Wl,--gc-sections',
  ],
};

void main(List<String> args) async {
  await build(args, (input, output) async {
    if (!input.config.buildCodeAssets) return;

    final builder = CBuilder.library(
      name: 'ctsa',
      // Asset id becomes `package:<packageName>/ctsa`. Must match the
      // `@DefaultAsset('package:tseries/ctsa')` in the Dart FFI layer.
      assetName: 'ctsa',
      sources: [
        'native/tseries_ctsa.c',
        ..._ctsaSources,
        ..._replacementSources,
      ],
      // Link the C math library where it is a separate library (see above).
      libraries: _mathLibrariesFor(input.config.code.targetOS),
      flags: [
        ..._floatingPointFlagsFor(input.config.code.targetOS),
        ..._sizeAndVisibilityFlagsFor(input.config.code.targetOS),
      ],
      // Header search paths: our shim header, the vendored ctsa headers (src/,
      // which is the real header the .c files compile against) and the
      // clean-room replacements (provides `brent.h`).
      includes: const [
        'native',
        'third_party/ctsa/src',
        'third_party/ctsa_replacements',
      ],
      defines: const {
        // NOTE: `_USE_MATH_DEFINES` (which exposes M_PI/M_SQRT2 in MSVC's
        // <math.h>) is intentionally NOT defined here: only the excluded wavelet
        // cluster used those constants. The ARIMA path uses ctsa's own PI2/PIVAL
        // macros. Restore it alongside the wavelet sources if they come back —
        // defining it here now would redefine the copy wavefilt.h sets itself.
        //
        // Silence MSVC's deprecation of strcpy/sprintf used by ctsa. These are
        // warnings, not correctness issues; keep the build output clean.
        '_CRT_SECURE_NO_WARNINGS': null,
      },
      // ctsa uses C99/C11 features (flexible array members, // comments).
      std: 'c11',
    );

    await builder.run(input: input, output: output);
  });
}
