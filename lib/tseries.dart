/// tseries — a framework-agnostic time-series analysis SDK for Dart.
///
/// This is the public barrel of the package: it re-exports the clean,
/// idiomatic **domain API** and nothing else. Consumers import only this file:
///
/// ```dart
/// import 'package:tseries/tseries.dart';
/// ```
///
/// ## Layering (why this file is thin on purpose)
///
/// The package is built in three layers. Only the top layer is public:
///
/// 1. **Raw bindings** — `ffigen`-generated bindings to the native C core
///    (ctsa), lives under `lib/src/ffi/bindings/`. Private to the package;
///    never exported. `Pointer` types stay here.
/// 2. **Safe wrapper** — Dart code under `lib/src/ffi/wrapper/` that owns
///    native memory, checks nullability/bounds, and translates C error codes
///    into Dart exceptions. Also private. No `Pointer` crosses this boundary.
/// 3. **Domain API** — the pure, idiomatic surface under `lib/src/domain/`.
///    This is the ONLY layer re-exported here.
///
/// The pure-Dart domain primitives and the native-backed ARIMA/SARIMA API are
/// both live. Consumers cannot tell which operations are native — and that is
/// the point.
library;

export 'src/domain/arima.dart';
export 'src/domain/auto_arima.dart';
export 'src/domain/decomposition.dart';
export 'src/domain/holt.dart';
export 'src/domain/information_criterion.dart';
export 'src/domain/kpss.dart' hide kpssStatistic, isConstantSeries, maxAbs;
export 'src/domain/linear.dart';
export 'src/domain/moving_average.dart';
export 'src/domain/native_health.dart';
export 'src/domain/polynomial_roots.dart';
export 'src/domain/sarimax.dart';
