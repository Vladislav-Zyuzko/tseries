/// Heavy autoArima gates (tagged `heavy`, skipped by the default run; run
/// with `dart test --tags heavy --run-skipped`):
///
/// * E — determinism: the journal of the six reference series in both modes,
///   plus short seasonal boundary series, is byte-identical in two fresh
///   processes and in a third one that processes the series in reverse
///   order (different heap history).
/// * F — the length rule of spec §5 is stronger than every length guard of
///   the fit path: for every candidate form within the default bounds, at
///   the rule's threshold − 1, threshold and threshold + 1, the fit path
///   never raises its own length error.
@Tags(['heavy'])
library;

import 'dart:io';

import 'package:test/test.dart';

import 'support/child_vm.dart';

Future<String> _journal(
  Directory dir,
  String name, {
  bool reverse = false,
}) async {
  final path = '${dir.path}${Platform.pathSeparator}$name.txt';
  // A child VM that does not re-run the build hooks: `dart run` would fail to
  // re-bundle the native library while this test process has it loaded (see
  // support/child_vm.dart).
  final r = await runInChildVm(
    'tool/auto_arima_journal.dart',
    arguments: ['--journal', path, '--boundary', if (reverse) '--reverse'],
  );
  if (r.exitCode != 0) {
    fail('journal process failed (${r.exitCode}): ${r.stderr}');
  }
  return File(path).readAsStringSync();
}

void main() {
  test(
    'E: journals byte-identical across fresh processes and run order',
    () async {
      final dir = Directory.systemTemp.createTempSync('auto_arima_e_');
      try {
        final a = await _journal(dir, 'a');
        final b = await _journal(dir, 'b');
        final c = await _journal(dir, 'c', reverse: true);
        expect(a.length, greaterThan(10000));
        expect(b, a);
        expect(c, a);
        // ignore: avoid_print
        print(
          'E: ${a.length} bytes, ${a.split('\n').length} lines, identical x3',
        );
      } finally {
        dir.deleteSync(recursive: true);
      }
    },
    timeout: const Timeout(Duration(minutes: 10)),
  );

  test(
    'F: length rule vs the fit path at threshold −1/0/+1',
    () async {
      // Run as a separate process (tool/auto_arima_length_sweep.dart): a
      // native crash must fail this test with the offending case named, not
      // take the test runner down.
      final dir = Directory.systemTemp.createTempSync('auto_arima_f_');
      try {
        final log = '${dir.path}${Platform.pathSeparator}sweep.txt';
        final r = await runInChildVm(
          'tool/auto_arima_length_sweep.dart',
          arguments: [log],
        );
        final lines = File(log).readAsLinesSync();
        final tail = lines.skip(lines.length > 6 ? lines.length - 6 : 0);
        // ignore: avoid_print
        print('F sweep exit ${r.exitCode}; last lines:\n${tail.join('\n')}');
        expect(r.exitCode, 0, reason: 'sweep process died; last case: $tail');
        expect(lines.last, 'DONE');
        expect(
          lines.where((l) => l.contains('DEFECT') || l.contains('NOT APPLIED')),
          isEmpty,
        );
      } finally {
        dir.deleteSync(recursive: true);
      }
    },
    // The 13 608-case sweep alone takes ~19 min on the Windows dev machine
    // (measured 2026-10-06), and heavy test files run in parallel.
    timeout: const Timeout(Duration(minutes: 60)),
  );
}
