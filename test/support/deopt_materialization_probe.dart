/// Child-process probe for test/no_typed_data_views_test.dart: the minimal
/// sequence that crashed the process before the wrapper stopped using
/// `Pointer.asTypedList` views.
///
/// Run under `--no-background-compilation --optimization-counter-threshold=10`
/// so the JIT optimizes `sarimaxFit` during the warm-up fits (no regressors)
/// and then lazily deoptimizes it on the first fit WITH a regressor (the
/// drift column), while the old code held a scalar-replaced `diag` view. The
/// VM rebuilt that view pointing to twice the real address, and the next read
/// faulted (0xC0000005). Prints `PROBE OK` and exits 0 when the sequence
/// completes.
library;

import 'dart:io';
import 'dart:typed_data';

import 'package:tseries/tseries.dart';

import 'arima_simulator.dart';

void main() {
  final warm = Float64List.fromList(
    simulateArima(n: 30, seed: 930, ar: [0.4], seasonalMa: [0.4], period: 4),
  );
  for (var i = 0; i < 20; i++) {
    fitSarimax(
      warm,
      const {},
      order: const ArimaOrder(p: 1, d: 0, q: 1),
      includeMean: false,
      horizon: 1,
    );
  }
  final drifted = Float64List.fromList(
    simulateArima(
      n: 9,
      seed: 279,
      ar: [0.4],
      seasonalMa: [0.4],
      period: 4,
      d: 1,
    ),
  );
  final fit = fitSarimax(
    drifted,
    const {},
    order: const ArimaOrder(p: 0, d: 1, q: 0),
    includeDrift: true,
    horizon: 1,
  );
  stdout.writeln(
    'PROBE OK sigma2=${fit.sigma2} drift=${fit.regressors.single.coefficient}',
  );
}
