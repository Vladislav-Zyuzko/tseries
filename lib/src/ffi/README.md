# `lib/src/ffi/` — native boundary

This directory holds the **native (FFI) integration** for tseries.

## Golden rule

**Nothing under `lib/src/ffi/` is exported from `lib/tseries.dart`, and no
`dart:ffi` `Pointer` ever appears in a public API.** The FFI machinery is an
implementation detail; consumers (including Sweet Limit) cannot tell native code
exists.

## Layers

```
lib/src/ffi/
  bindings/   Layer 1 — raw bindings to the C shim (native/tseries_ctsa.h),
              generated from ffigen.yaml. Private, never hand-edited, never
              exported. The only place Pointer types live.

  wrapper/    Layer 2 — safe Dart wrapper. Owns native memory (Arena / using,
              freed on every path including throws), null/bounds-checks native
              returns, scans results for NaN/Inf, and translates C status codes
              into typed Dart exceptions (tseries_exceptions.dart). No Pointer
              escapes this layer.
```

The clean, idiomatic result is exposed from `lib/src/domain/arima.dart`, which
is the only layer the public barrel re-exports.

## Memory ownership

Every allocation is scoped to an `Arena` via `using(...)`, so it is freed when
the block exits — including on exceptions. No pointer leaves its owning scope.
The native shim owns and frees its own temporaries (see `native/`).
