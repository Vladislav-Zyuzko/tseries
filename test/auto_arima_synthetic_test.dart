/// Gate C of the autoArima specification (§9.C, revisions 5–6), the local
/// part — run in a child VM by tool/auto_arima_synthetic_study.dart (see
/// there for what is computed). Tagged `heavy`: ~1 h at 100 repetitions.
///
/// AUTO_ARIMA_C_REPS lowers the repetitions for a quick partial run;
/// AUTO_ARIMA_C_EXPORT writes the series and our results for the R run
/// (doc/auto_arima_r_exchange.md).
@Tags(['heavy'])
library;

import 'dart:io';

import 'package:test/test.dart';

import 'support/child_vm.dart';

void main() {
  test('C: synthetic study, local gates', () async {
    final env = Platform.environment;
    final r = await runInChildVm(
      'tool/auto_arima_synthetic_study.dart',
      arguments: [
        '--reps',
        env['AUTO_ARIMA_C_REPS'] ?? '100',
        if (env['AUTO_ARIMA_C_EXPORT'] != null) ...[
          '--export',
          env['AUTO_ARIMA_C_EXPORT']!,
        ],
      ],
    );
    final out = '${r.stdout}';
    // ignore: avoid_print
    print(out);
    expect(r.exitCode, 0, reason: '${r.stderr}');
    expect(out, contains('gate notes: none'));
  }, timeout: const Timeout(Duration(hours: 3)));
}
