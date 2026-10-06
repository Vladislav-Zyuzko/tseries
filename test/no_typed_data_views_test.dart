/// Regression test for the 0xC0000005 crash found by the autoArima length
/// sweep (tool/auto_arima_length_sweep.dart, case #1948 after #637..#1947).
///
/// Root cause (Dart VM 3.10, JIT only): the optimizer scalar-replaces a
/// `Pointer.asTypedList` view that does not escape. When the optimized code is
/// then deoptimized, the VM rebuilds the view by storing its data address as a
/// tagged integer, so the view points to twice the real address. In
/// `sarimaxFit` the first fit with a regressor lazily deoptimized the closure
/// (class-hierarchy invalidation on the first `ExogColumnDefect` use) while
/// its `diag` view was live; the next `diag[0]` read faulted.
///
/// The wrapper now accesses native buffers only element by element through
/// the `Pointer`. Two guards:
///
/// * a source scan: no `asTypedList` call anywhere in `lib/` (a write view
///   would make the same defect a wild WRITE);
/// * the minimal crashing sequence, in a child VM with the flags that make it
///   deterministic. It crashed before the fix and must exit 0 now.
///
/// The child must not run build hooks: `dart run` would re-bundle the native
/// asset into `.dart_tool/lib/`, which fails while this test process (or any
/// other test file running in parallel) has that library loaded. So the probe
/// is run through `runInChildVm` (test/support/child_vm.dart): kernel from the
/// SDK's `gen_kernel` with the native-assets mapping `dart test` already
/// wrote, executed by the bare VM. Nothing is built or copied.
library;

import 'dart:io';

import 'package:test/test.dart';

import 'support/child_vm.dart';

void main() {
  test('lib/ never creates Pointer.asTypedList views', () {
    final offenders = <String>[];
    for (final entity in Directory('lib').listSync(recursive: true)) {
      if (entity is! File || !entity.path.endsWith('.dart')) continue;
      final lines = entity.readAsLinesSync();
      for (var i = 0; i < lines.length; i++) {
        final code = lines[i].split('//').first;
        if (code.contains('asTypedList')) {
          offenders.add('${entity.path}:${i + 1}: ${lines[i].trim()}');
        }
      }
    }
    expect(offenders, isEmpty);
  });

  test(
    'deoptimizing sarimaxFit on its first regressor fit does not crash',
    () async {
      final result = await runInChildVm(
        'test/support/deopt_materialization_probe.dart',
        vmFlags: const [
          '--no-background-compilation',
          '--optimization-counter-threshold=10',
        ],
      );
      expect(
        result.exitCode,
        0,
        reason:
            'probe process failed\nstdout:\n${result.stdout}\n'
            'stderr:\n${result.stderr}',
      );
      expect(result.stdout as String, contains('PROBE OK'));
    },
    timeout: const Timeout(Duration(minutes: 5)),
  );
}
