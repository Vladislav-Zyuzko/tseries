// SPDX-License-Identifier: BSD-3-Clause
//
// Gate E0: counts fused floating-point multiply-add instructions in built
// native binaries. The native core is built with floating-point contraction
// off (hook/build.dart, `_floatingPointFlagsFor`), so the expected count is 0
// on every ABI; anything else means a toolchain or flag change re-enabled FMA
// and cross-platform model selection is no longer what was verified
// (doc/platform_verification.md).
//
// Usage (from packages/tseries):
//
//   dart run tool/fma_count.dart [--objdump <llvm-objdump>] [--verbose]
//       [--hook-objects <build dir>] <binary>...
//
// <binary> is a built libctsa.so (Android/Linux, any ABI), libctsa.dylib
// (macOS/iOS) or an object file. The architecture is read from the binary.
// Windows: pass the hook's build directory with --hook-objects instead of
// ctsa.dll — the DLL links the C runtime statically, and the UCRT math
// functions contain FMA3 variants of their own (see _hookObjects).
//
// Note: `dart run` inside this package runs the build hook first; do not run
// it while `dart test` is running here (both copy .dart_tool/lib/ctsa.dll).
// llvm-objdump is taken from --objdump, else $LLVM_OBJDUMP, else the newest
// Android NDK under $ANDROID_NDK_HOME / $ANDROID_HOME / $ANDROID_SDK_ROOT /
// the default SDK location, else PATH.
//
// Exit code: 0 if every binary has 0 fused instructions, 1 if any has some,
// 2 on a usage or tool error.

import 'dart:io';

/// Fused-multiply-add mnemonics per architecture, as llvm-objdump prints them.
///
/// * aarch64 — scalar `fmadd/fmsub/fnmadd/fnmsub`, vector and by-element
///   `fmla/fmls` (also SVE `fmad/fmsb/fnmad/fnmsb/fnmla/fnmls`).
/// * x86 / x86_64 — FMA3 `vfmadd*/vfmsub*/vfnmadd*/vfnmsub*/vfmaddsub*/
///   vfmsubadd*` and the AMD FMA4 spellings (same prefixes).
/// * 32-bit arm — VFPv4/NEON `vfma/vfms/vfnma/vfnms`. `vmla/vmls` are NOT
///   fused (the product is rounded first) and are not counted.
final Map<String, RegExp> _fusedByArch = {
  'aarch64': RegExp(r'^(fn?m(add|sub)|fml[as]|fn?m(ad|sb)|fnml[as])$'),
  'x86_64': RegExp(r'^vfn?m(add|sub|addsub|subadd)'),
  'x86': RegExp(r'^vfn?m(add|sub|addsub|subadd)'),
  'arm': RegExp(r'^vfn?m[as](\.|$)'),
};

/// Maps llvm-objdump's `file format` / `architecture` strings to a key of
/// [_fusedByArch].
String? _archKey(String header) {
  final h = header.toLowerCase();
  if (h.contains('aarch64') || h.contains('arm64')) return 'aarch64';
  if (h.contains('x86-64') || h.contains('x86_64')) return 'x86_64';
  if (h.contains('i386') ||
      h.contains('coff-i386') ||
      h.contains('elf32-i386')) {
    return 'x86';
  }
  if (h.contains('arm') || h.contains('thumb')) return 'arm';
  return null;
}

final RegExp _insnLine = RegExp(r'^\s*[0-9a-f]+:\s+([a-z][a-z0-9.]*)');
final RegExp _symbolLine = RegExp(r'^[0-9a-f]+ <(.+)>:$');

Future<void> main(List<String> args) async {
  String? objdump;
  var verbose = false;
  final binaries = <String>[];
  for (var i = 0; i < args.length; i++) {
    final a = args[i];
    if (a == '--objdump' && i + 1 < args.length) {
      objdump = args[++i];
    } else if (a == '--hook-objects' && i + 1 < args.length) {
      final objs = _hookObjects(args[++i]);
      if (objs == null) {
        exitCode = 2;
        return;
      }
      binaries.addAll(objs);
    } else if (a == '--verbose' || a == '-v') {
      verbose = true;
    } else if (a == '--help' || a == '-h') {
      stdout.writeln(
        'usage: dart run tool/fma_count.dart [--objdump <llvm-objdump>] '
        '[--verbose] [--hook-objects <build dir>] <binary>...',
      );
      return;
    } else {
      binaries.add(a);
    }
  }
  if (binaries.isEmpty) {
    stderr.writeln('fma_count: no binaries given (see --help)');
    exitCode = 2;
    return;
  }
  objdump ??= _findObjdump();
  if (objdump == null) {
    stderr.writeln(
      'fma_count: llvm-objdump not found; pass --objdump or set LLVM_OBJDUMP',
    );
    exitCode = 2;
    return;
  }

  var anyFused = false;
  for (final path in binaries) {
    if (!File(path).existsSync()) {
      stderr.writeln('fma_count: $path: no such file');
      exitCode = 2;
      return;
    }
    final head = await Process.run(objdump, ['-f', path]);
    if (head.exitCode != 0) {
      stderr.writeln('fma_count: $path: ${head.stderr}');
      exitCode = 2;
      return;
    }
    final arch = _archKey(head.stdout as String);
    if (arch == null) {
      stderr.writeln('fma_count: $path: unknown architecture:\n${head.stdout}');
      exitCode = 2;
      return;
    }
    final dis = await Process.run(objdump, [
      '-d',
      '--no-show-raw-insn',
      path,
    ], stdoutEncoding: systemEncoding);
    if (dis.exitCode != 0) {
      stderr.writeln('fma_count: $path: ${dis.stderr}');
      exitCode = 2;
      return;
    }
    final fused = _fusedByArch[arch]!;
    var total = 0;
    final byMnemonic = <String, int>{};
    final bySymbol = <String, int>{};
    var symbol = '?';
    for (final line in (dis.stdout as String).split('\n')) {
      final s = _symbolLine.firstMatch(line.trimRight());
      if (s != null) {
        symbol = s.group(1)!;
        continue;
      }
      final m = _insnLine.firstMatch(line);
      if (m == null) continue;
      total++;
      final mnemonic = m.group(1)!;
      if (fused.hasMatch(mnemonic)) {
        byMnemonic.update(mnemonic, (n) => n + 1, ifAbsent: () => 1);
        bySymbol.update(symbol, (n) => n + 1, ifAbsent: () => 1);
      }
    }
    final count = byMnemonic.values.fold(0, (a, b) => a + b);
    anyFused |= count > 0;
    final kinds =
        (byMnemonic.entries.toList()
              ..sort((a, b) => b.value.compareTo(a.value)))
            .map((e) => '${e.key} ${e.value}')
            .join(', ');
    stdout.writeln(
      '$path\n  arch $arch, ${File(path).lengthSync()} B, '
      '$total instructions, fused $count${count > 0 ? ' ($kinds)' : ''}',
    );
    if (verbose && bySymbol.isNotEmpty) {
      final rows = bySymbol.entries.toList()
        ..sort((a, b) => b.value.compareTo(a.value));
      for (final r in rows) {
        stdout.writeln('    ${r.value.toString().padLeft(5)}  ${r.key}');
      }
    }
  }
  stdout.writeln(anyFused ? 'E0 FAILED: fused instructions found' : 'E0 OK');
  if (anyFused) exitCode = 1;
}

/// The object files (`<name>.obj` or `<name>.o`) in [dir] for exactly the
/// C sources that hook/build.dart compiles today, so stale objects left in a
/// build directory by older source lists are not counted.
///
/// Used for Windows: ctsa.dll links the C runtime statically, and the UCRT's
/// own math functions carry FMA3 variants (`exp_fma`, `log_fma`, …) that it
/// selects at run time by CPUID — so the DLL as a whole is never 0, while our
/// objects must be.
List<String>? _hookObjects(String dir) {
  final hook = File.fromUri(Platform.script.resolve('../hook/build.dart'));
  if (!hook.existsSync()) {
    stderr.writeln('fma_count: ${hook.path} not found');
    return null;
  }
  final sources = RegExp(
    r"'([^']+)\.c'",
  ).allMatches(hook.readAsStringSync()).map((m) => m.group(1)!.split('/').last);
  final out = <String>[];
  for (final name in sources) {
    final candidates = ['$name.obj', '$name.o']
        .map((f) => File([dir, f].join(Platform.pathSeparator)))
        .where((f) => f.existsSync());
    if (candidates.isEmpty) {
      stderr.writeln('fma_count: no object for $name.c in $dir');
      return null;
    }
    out.add(candidates.first.path);
  }
  return out;
}

String? _findObjdump() {
  final exe = Platform.isWindows ? 'llvm-objdump.exe' : 'llvm-objdump';
  final env = Platform.environment;
  final fromEnv = env['LLVM_OBJDUMP'];
  if (fromEnv != null && File(fromEnv).existsSync()) return fromEnv;

  final ndkRoots = <String>[
    ?env['ANDROID_NDK_HOME'],
    ?env['ANDROID_NDK_ROOT'],
    for (final sdk in [
      ?env['ANDROID_HOME'],
      ?env['ANDROID_SDK_ROOT'],
      if (Platform.isWindows && env['LOCALAPPDATA'] != null)
        '${env['LOCALAPPDATA']}\\Android\\Sdk',
      if (env['HOME'] != null) '${env['HOME']}/Library/Android/sdk',
      if (env['HOME'] != null) '${env['HOME']}/Android/Sdk',
    ])
      ..._ndkVersions('$sdk${Platform.pathSeparator}ndk'),
  ];
  for (final ndk in ndkRoots) {
    final prebuilt = Directory(
      [ndk, 'toolchains', 'llvm', 'prebuilt'].join(Platform.pathSeparator),
    );
    if (!prebuilt.existsSync()) continue;
    for (final host in prebuilt.listSync().whereType<Directory>()) {
      final f = File([host.path, 'bin', exe].join(Platform.pathSeparator));
      if (f.existsSync()) return f.path;
    }
  }

  final which = Process.runSync(Platform.isWindows ? 'where' : 'which', [
    'llvm-objdump',
  ]);
  if (which.exitCode == 0) {
    final first = (which.stdout as String).split(RegExp(r'\r?\n')).first.trim();
    if (first.isNotEmpty) return first;
  }
  return null;
}

/// NDK installations under `<sdk>/ndk`, newest version first.
List<String> _ndkVersions(String ndkDir) {
  final d = Directory(ndkDir);
  if (!d.existsSync()) return const [];
  final versions =
      d.listSync().whereType<Directory>().map((e) => e.path).toList()
        ..sort((a, b) => _compareVersions(b, a));
  return versions;
}

int _compareVersions(String a, String b) {
  List<int> parts(String p) => p
      .split(RegExp(r'[\\/]'))
      .last
      .split('.')
      .map((s) => int.tryParse(s) ?? 0)
      .toList();
  final pa = parts(a), pb = parts(b);
  for (var i = 0; i < pa.length && i < pb.length; i++) {
    if (pa[i] != pb[i]) return pa[i].compareTo(pb[i]);
  }
  return pa.length.compareTo(pb.length);
}
