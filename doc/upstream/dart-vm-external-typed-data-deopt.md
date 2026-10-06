<!--
Draft of an issue for https://github.com/dart-lang/sdk. NOT SENT — whether and
when to send it is the project owner's decision. Written 2026-10-06 against
Dart 3.10.9 (stable, windows_x64, via Flutter 3.38.10).
-->

# [vm] Deoptimization materializes a scalar-replaced `Pointer.asTypedList` view with a tagged data address

## Summary

When optimized JIT code that holds a non-escaping external typed-data view
(`Pointer<Double>.asTypedList(n)`, allocation-sunk by the optimizer) is
deoptimized, the deoptimizer rebuilds the `_ExternalFloat64Array` through the
generic instance path of `DeferredObject::Fill()` (`runtime/vm/deferred_objects.cc`):
the field at `PointerBase::data_offset()` has no `Field`, so the value is stored
with `SetFieldAtOffset(offset, value)` — i.e. the **tagged** integer. For an
address that fits in a Smi the view's untagged `data_` becomes `address << 1`,
and the next element access through the view reads (or writes) wild memory.
`kPointerCid` has a dedicated case that unboxes the address
(`pointer.SetNativeAddress(Integer::Cast(...).Value())`); external typed data
has none (also on `main` as of 2026-10-06).

AOT is not affected (no deoptimization).

## Evidence

`--trace-deoptimization-verbose`, lazy deopt (CHA invalidation) of the
optimized closure that holds `final diag = diagPtr.asTypedList(4);`:

```
materializing instance of Library:'dart:typed_data' Class: _ExternalFloat64Array@8027147 (…, 2 fields)
    null Field @ offset(8) <- 1735523024032      # == diagPtr.address (0x194153bbca0)
    null Field @ offset(16) <- 4                  # length
```

The access violation that follows is in optimized `_ExternalFloat64Array.[]`:
`mov rdx,[rcx+7]` (data_) → `movsd xmm0,[rdx+rsi*8]` faults with
`rdx = 2 × diagPtr.address` (captured with a vectored exception handler:
fault address `0x3d9d6401900` for a buffer at `0x1eceb200c80`).

## Reproduction

In this repository (`packages/tseries`, before the workaround — any revision
whose `lib/src/ffi/wrapper/ctsa_native.dart` still uses `asTypedList`):

```
dart --no-background-compilation --optimization-counter-threshold=10 \
  run test/support/deopt_materialization_probe.dart
```

20 fits without regressors optimize `sarimaxFit`'s closure; the first fit with
a regressor runs `ExogColumnDefect` code for the first time, which lazily
deoptimizes it while the view is live → 0xC0000005. With default flags the same
happens deterministically after ~1 900 fits (the autoArima length sweep).
With `--optimization-counter-threshold=-1` it never happens.

## Workaround

Do not keep `asTypedList` views in locals of hot code; access native memory
element by element through the `Pointer` (`ptr[i]`), whose materialization is
handled correctly.
