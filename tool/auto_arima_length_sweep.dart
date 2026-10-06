// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// Gate F sweep of the autoArima specification (§5 rule 1, §9.F): for every
// candidate form (p, q ≤ 5, P, Q ≤ 2, constant) at s ∈ {4, 12}, D ∈ {0, 1},
// d ∈ {0, 1}, evaluates the candidate exactly as autoArima does at the
// length rule's threshold − 1, threshold and threshold + 1, and checks that
// the fit path never raises its own length error (which would be a defect of
// the rule).
//
// Each case is written to the log BEFORE it runs and flushed, so if the
// process dies the last line names the case. After each case all process
// heaps are validated (Windows only: kernel32 HeapValidate), so heap
// corruption is attributed to the first case that causes it.
//
// Usage (from packages/tseries):
//   dart run tool/auto_arima_length_sweep.dart <log> [from] [to] [skip,list]
// Exit code 0 and a final "DONE" line: sweep completed; the log's summary
// line counts fit-path length errors ("defects=").
import 'dart:ffi' as ffi;
import 'dart:io';
import 'dart:typed_data';

import 'package:ffi/ffi.dart';
import 'package:tseries/src/domain/internal/auto_arima_candidate.dart';
import 'package:tseries/tseries.dart';

import '../test/support/arima_simulator.dart';

bool Function() _heapValidator() {
  if (!Platform.isWindows) return () => true;
  final k32 = ffi.DynamicLibrary.open('kernel32.dll');
  final validate = k32
      .lookupFunction<
        ffi.Int32 Function(
          ffi.Pointer<ffi.Void>,
          ffi.Uint32,
          ffi.Pointer<ffi.Void>,
        ),
        int Function(ffi.Pointer<ffi.Void>, int, ffi.Pointer<ffi.Void>)
      >('HeapValidate');
  final heaps = k32
      .lookupFunction<
        ffi.Uint32 Function(ffi.Uint32, ffi.Pointer<ffi.Pointer<ffi.Void>>),
        int Function(int, ffi.Pointer<ffi.Pointer<ffi.Void>>)
      >('GetProcessHeaps');
  final buf = calloc<ffi.Pointer<ffi.Void>>(256);
  return () {
    final count = heaps(256, buf);
    for (var i = 0; i < count && i < 256; i++) {
      if (validate(buf[i], 0, ffi.nullptr) == 0) return false;
    }
    return true;
  };
}

void main(List<String> args) {
  final log = File(args[0]).openSync(mode: FileMode.write);
  final from = args.length > 1 ? int.parse(args[1]) : 0;
  final to = args.length > 2 ? int.parse(args[2]) : 1 << 30;
  final skip = args.length > 3 && args[3].isNotEmpty
      ? args[3].split(',').map(int.parse).toSet()
      : <int>{};
  final heapOk = _heapValidator();
  var idx = -1;
  var evaluated = 0;
  var refused = 0;
  var defects = 0;
  void write(String s) {
    log.writeStringSync(s);
    log.flushSync();
  }

  for (final s in [4, 12]) {
    for (final sD in [0, 1]) {
      for (final d in [0, 1]) {
        for (var p = 0; p <= 5; p++) {
          for (var q = 0; q <= 5; q++) {
            for (var sp = 0; sp <= 2; sp++) {
              for (var sq = 0; sq <= 2; sq++) {
                for (final c in d + sD < 2 ? [false, true] : [false]) {
                  final key = (p: p, q: q, sp: sp, sq: sq, c: c);
                  var t = 1;
                  while (lengthRuleViolation(key, t: t, period: s) != null) {
                    t++;
                  }
                  for (final tt in [t - 1, t, t + 1]) {
                    idx++;
                    if (idx < from || idx > to || skip.contains(idx)) continue;
                    final n = tt + d + s * sD;
                    final y = simulateArima(
                      n: n,
                      seed: 31 * n + p,
                      ar: [0.4],
                      seasonalMa: [0.4],
                      period: s,
                      d: d,
                      seasonalD: sD,
                    );
                    write('#$idx s=$s D=$sD d=$d $key T=$tt n=$n\n');
                    final out = evaluateCandidate(
                      Float64List.fromList(y),
                      key,
                      d: d,
                      seasonalD: sD,
                      period: s,
                      criterion: InformationCriterion.aicc,
                      horizon: 1,
                      rootMargin: 1.01,
                    );
                    evaluated++;
                    final bad =
                        tt >= t &&
                        ((out.detail ?? '').contains('DEFECT') ||
                            out.verdict == CandidateVerdict.rejectedFitError);
                    if (tt < t) {
                      refused++;
                      if (out.verdict != CandidateVerdict.rejectedTooShort ||
                          out.fitCalled) {
                        defects++;
                        write('  RULE NOT APPLIED\n');
                      }
                    }
                    if (bad) defects++;
                    final ok = heapOk();
                    write(
                      '  -> ${out.verdict.name} status=${out.ctsaStatus} '
                      'heap=${ok ? 'ok' : 'CORRUPT'}'
                      '${bad ? ' DEFECT ${out.detail}' : ''}\n',
                    );
                    if (!ok) {
                      write('HEAP CORRUPT after #$idx\n');
                      log.closeSync();
                      exit(3);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  write(
    'SUMMARY evaluated=$evaluated refusedByRule=$refused defects=$defects\n'
    'DONE\n',
  );
  log.closeSync();
}
