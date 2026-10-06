/// Shared loaders for the autoArima reference series and the R reference
/// orders (test/fixtures/auto_arima/). Only numbers from external tools are
/// read here.
library;

import 'dart:io';

/// The six reference series and their seasonal periods.
const autoArimaReferenceSeries = <(String, int)>[
  ('air_passengers_log', 12),
  ('wwwusage', 1),
  ('lynx', 1),
  ('lake_huron', 1),
  ('nile', 1),
  ('us_acc_deaths', 12),
];

List<double> loadReferenceSeries(String name) =>
    File('test/fixtures/auto_arima/series_$name.csv')
        .readAsLinesSync()
        .where((l) => l.trim().isNotEmpty)
        .map(double.parse)
        .toList();

/// One order chosen by R's auto.arima (used as an external tool).
typedef ROrder = ({
  int p,
  int d,
  int q,
  int sp,
  int sd,
  int sq,
  int period,
  String constant,
});

/// `name|mode` → order, from reference_r_forecast.txt (format in its header).
Map<String, ROrder> loadROrders() {
  final out = <String, ROrder>{};
  for (final line in File(
    'test/fixtures/auto_arima/reference_r_forecast.txt',
  ).readAsLinesSync()) {
    if (line.isEmpty || line.startsWith('#')) continue;
    final f = line.split('|');
    out['${f[0]}|${f[1]}'] = (
      p: int.parse(f[2]),
      d: int.parse(f[3]),
      q: int.parse(f[4]),
      sp: int.parse(f[5]),
      sd: int.parse(f[6]),
      sq: int.parse(f[7]),
      period: int.parse(f[8]),
      constant: f[9],
    );
  }
  return out;
}

String describeROrder(ROrder o) =>
    '(${o.p},${o.d},${o.q})(${o.sp},${o.sd},${o.sq})[${o.period}] '
    '${o.constant}';
