import 'package:tseries/tseries.dart';

void main() {
  final glucose = [5.4, 5.8, 6.1, 6.0, 5.7, 5.9, 6.3, 6.1, 5.8, 6.0, 6.2, 6.4];

  // Pure-Dart: smooth the series with a 3-point simple moving average.
  final smoothed = simpleMovingAverage(glucose, window: 3);
  print('series:   $glucose');
  print('SMA(3):   $smoothed');

  // Native (ctsa/FFI): a fixed-order ARIMA forecast with a fallback ladder.
  // Consumers cannot tell this call is backed by native code.
  try {
    final result = arimaForecast(
      glucose,
      orders: const [
        ArimaOrder(p: 1, d: 1, q: 0),
        ArimaOrder(p: 0, d: 1, q: 0),
      ],
      horizon: 3,
    );
    print('order:    ${result.order}');
    print('forecast: ${result.point}');
    print('stderr:   ${result.standardErrors}');
  } on TseriesNumericException catch (e) {
    // The wrapper refuses to hand back non-finite forecasts.
    print('no forecast available: ${e.message}');
  }

  // Automatic order selection: d by KPSS, then a stepwise search over
  // (p, q, constant) minimising AICc. Every candidate is in `trace`.
  final series = [for (var t = 0; t < 60; t++) 10 + 0.1 * t + (7 * t % 5) / 4];
  try {
    final auto = autoArima(series, horizon: 3);
    print(
      'autoArima: ${auto.order} constant=${auto.constant} '
      '(${auto.fitsEvaluated} fits)',
    );
    print('forecast: ${auto.forecast}');
    print('95% band: ${auto.lower} .. ${auto.upper}');
  } on AutoArimaNoModelException catch (e) {
    print('no model selected (${e.reason.name}): ${e.message}');
  }
}
