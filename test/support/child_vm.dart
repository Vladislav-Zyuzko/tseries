/// Runs a Dart program of this package in a child VM WITHOUT build hooks.
///
/// `dart run` in a child re-bundles the `ctsa` native asset into
/// `.dart_tool/lib/`, which fails (`PathExistsException`) while the test
/// process — or any test file running in parallel in it — has that library
/// loaded. Instead the program is compiled to kernel with the SDK's
/// `gen_kernel`, given the native-assets mapping `dart test` already wrote
/// (`.dart_tool/native_assets.yaml`, pointing at the library the tests
/// themselves use), and run by the bare VM (`dartvm`). Nothing is built or
/// copied, and [vmFlags] (e.g. JIT tuning flags) can be passed.
library;

import 'dart:io';

String _join(List<String> parts) => parts.join(Platform.pathSeparator);

/// Compiles [script] (a path relative to the package root) once per call and
/// runs it with [arguments]. Fails with a descriptive [StateError] when the
/// SDK layout or the native-assets mapping is not what this relies on.
Future<ProcessResult> runInChildVm(
  String script, {
  List<String> arguments = const [],
  List<String> vmFlags = const [],
}) async {
  final exe = Platform.isWindows ? '.exe' : '';
  final bin = File(Platform.resolvedExecutable).parent.path;
  final sdk = Directory(bin).parent.path;
  final aotRuntime = _join([bin, 'dartaotruntime$exe']);
  final genKernel = _join([bin, 'snapshots', 'gen_kernel_aot.dart.snapshot']);
  final platform = _join([sdk, 'lib', '_internal', 'vm_platform_strong.dill']);
  final vm = _join([bin, 'dartvm$exe']);
  final packages = _join(['.dart_tool', 'package_config.json']);
  final nativeAssets = _join(['.dart_tool', 'native_assets.yaml']);
  for (final f in [
    aotRuntime,
    genKernel,
    platform,
    vm,
    packages,
    nativeAssets,
  ]) {
    if (!File(f).existsSync()) {
      throw StateError('child VM: required file missing: $f');
    }
  }

  final tmp = Directory.systemTemp.createTempSync('tseries_child_vm');
  try {
    final dill = _join([tmp.path, 'program.dill']);
    final compile = await Process.run(aotRuntime, [
      genKernel,
      '--platform',
      platform,
      '--packages',
      packages,
      '--native-assets',
      nativeAssets,
      '-o',
      dill,
      script,
    ]);
    if (compile.exitCode != 0) {
      throw StateError(
        'child VM: gen_kernel failed for $script\n'
        '${compile.stdout}\n${compile.stderr}',
      );
    }
    return await Process.run(vm, [
      ...vmFlags,
      '--packages=$packages',
      dill,
      ...arguments,
    ]);
  } finally {
    tmp.deleteSync(recursive: true);
  }
}
