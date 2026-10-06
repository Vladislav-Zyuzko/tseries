# `hook/` — native build hook

`build.dart` compiles the vendored ctsa C sources (`third_party/ctsa/`), the
clean-room replacements for removed upstream units
(`third_party/ctsa_replacements/`) and the shim (`native/tseries_ctsa.c`) into a
single `ctsa` code asset, using `package:native_toolchain_c`. No OS-specific
build files are required.

On ELF/Mach-O targets it builds with hidden visibility and linker dead-code
elimination (`_sizeAndVisibilityFlagsFor`): only the `tseries_*` shim functions
are exported, and unreachable ctsa code is dropped. Windows needs no flags
(exports are explicit `__declspec(dllexport)`; `/O2` + default `/OPT:REF`).

The asset id is `package:tseries/ctsa`, which must match the
`@DefaultAsset('package:tseries/ctsa')` in the generated bindings
(`lib/src/ffi/bindings/`).

Build hooks / code assets are **stable** as of Dart 3.10 (Flutter 3.38); this
package targets `sdk: ^3.10.0`, so the native build runs on the standard stable
SDK with no `--enable-experiment` flag and no beta channel.

Android (`arm64-v8a`, `armeabi-v7a`, `x86_64`) is built and verified on-device;
`native_toolchain_c` handles the cross-compilation. iOS (device arm64 +
simulator, xcframework) is the next target.

## Linking is not a build-time error — verify per ABI

The one trap this hook exists to remember: **a link mistake here does not fail the
build.** ctsa needs libm (`log`, `exp`, `pow`, `log10`, `sin`), which is a separate
shared library on Android/Linux but lives in the CRT on Windows. Omitting it still
produces a `libctsa.so` that compiles, exports every symbol and packages into the
APK — and then dies at `dlopen` on-device with `cannot locate symbol "log"`. The
Windows host suite passes throughout and cannot catch it. See `_mathLibrariesFor`
in `build.dart`.

So when you touch this hook or add an ABI, do not trust `dart test`. Check the
artifact itself:

```sh
llvm-readelf -d libctsa.so | grep NEEDED        # expect libm.so on Android/Linux
llvm-nm --dynamic --undefined-only libctsa.so   # expect no U log/exp/pow/sin/log10
llvm-nm --dynamic --defined-only libctsa.so     # expect exactly the 5 tseries_* exports
```

These are ABI-independent and prove the link even for an ABI you have no device
for. At runtime, `checkNativeCore()` (public API) is the one-call on-device
equivalent: it throws `TseriesNativeUnavailableException` when the asset cannot
be loaded, as distinct from a merely difficult fit.
