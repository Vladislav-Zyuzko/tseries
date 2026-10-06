# `autoArima` — clean-room record

This file is the written record required by §12 of the specification
("Спецификация: автоподбор порядка ARIMA на чистом Dart (`autoArima`)",
Arina, 2026-10-05, revisions 1–3). It states what the implementer of
`autoArima` read, and confirms what they did not open. It was written on
2026-10-05, **before** any `autoArima` code was written, and is updated only
by appending.

## Implementer

A fresh engineering agent instance (role: Vasily / FFI engineer, acting as the
clean-room implementer), 2026-10-05.

## Materials read

1. The specification itself:
   `.claude/agent-memory/arina-ml-timeseries/auto-arima-dart-spec-full.md`
   (revisions 1–3), in full.
2. The current state of this package (HEAD 8427748) — Dart sources under
   `lib/` (`tseries.dart`, `src/domain/*.dart`,
   `src/domain/internal/arima_support.dart`,
   `src/ffi/wrapper/ctsa_native.dart`, `src/ffi/wrapper/tseries_exceptions.dart`),
   `README.md`, the head of `CHANGELOG.md`, `pubspec.yaml`,
   `analysis_options.yaml`, and the existing tests
   `test/fit_path_comparability_test.dart`, `test/short_series_guard_test.dart`
   (style and conventions only), plus grep hits in
   `test/sarimax_golden_test.dart` for the coefficient sign conventions.
3. The fixture numbers in `test/fixtures/auto_arima/`
   (`reference_r_forecast.txt` and the six `series_*.csv`).
4. Published methodology as cited by the specification: Hyndman & Khandakar
   (2008), JSS 27(3), §3.1–3.2; Kwiatkowski, Phillips, Schmidt & Shin (1992)
   (statistic, Table 1, lag rule l4); Hyndman & Athanasopoulos, *Forecasting:
   Principles and Practice* (3rd ed.), §4.3, §9.1, §9.6, §9.8; Hurvich & Tsai
   (1989); Wang, Smith & Hyndman (2006). The implementation follows the
   formulas **as restated in the specification**; nothing beyond them was
   taken from these texts.
5. statsmodels used **as an external tool only** (installed in a throw-away
   virtualenv outside the repository; only its numeric output is stored in
   fixtures): `statsmodels.tsa.stattools.kpss` and
   `statsmodels.tsa.seasonal.seasonal_decompose` as independent references
   for gate A, and `numpy.roots` for the root check. Their API documentation
   (call signatures, parameter meaning) was consulted; their source code was
   not needed and was not read.

## Not opened

I confirm that I did **not** open, read or search:

* the sources of the R packages `forecast`, `tseries`, `urca`, or pmdarima;
* any other implementation of automatic ARIMA order selection in any
  language;
* the ctsa files removed from this repository (`autoutils.c`, `unitroot.c`,
  `seastest.c`, `stl.c`, the auto-ARIMA cluster of `ctsa.c`, `supsmu`), and
  no git history touching them (no `git show`, `git log -p`, `git checkout`
  or `git diff` of c58dc74 or earlier commits for those paths);
* the vendored C sources in `third_party/ctsa/src/` (not needed for this
  task);
* `tool/auto_arima_reference.R` (the script that produced the R reference
  numbers; only its numeric output in the fixture was used);
* the agent-memory files describing the licence audit and the removed code
  (`ctsa-license-audit.md`, `ctsa-latent-defects.md`,
  `tsforge-pubdev-readiness.md` and similar).

### Disclosure

The agent-memory **index** (`MEMORY.md` of the vasily-ffi-engineer memory) is
loaded into every session of this agent automatically. It contains one-line
titles of the files above (for example "auto cluster removed by name-only call
graph", "LGPL+ACM419+R helpers, auto-ARIMA"). Those one-liners carry no
algorithmic content and were not used for the implementation. The
implementation is derived from the specification only.

## Questions to the methodology owner

Where the specification is silent, the decision is marked `[наше]` in the
code and listed in the implementation report for Arina; none of them was
filled in from knowledge of how other implementations behave.

## Addendum (2026-10-06, during implementation)

* mpmath (BSD) was added to the external-tool virtualenv as a 60-digit
  arbiter for the root fixture: on clustered roots near the margin
  `numpy.roots` is itself only ~1e-8 accurate. Only numbers are stored.
* While investigating a process crash found by the length-rule sweep
  (`tool/auto_arima_length_sweep.dart`), the forecasting helper of the
  package's own native shim (`native/tseries_ctsa.c`, `shim_forecast`) was
  read. The shim is part of the current package and explicitly allowed; no
  vendored ctsa source was opened.
* Nothing from the "Not opened" list above was opened at any point.

## Addendum 2 (2026-10-06, iteration on spec revisions 5–6)

* Read: the specification sections changed in revisions 5 and 6 (§1, §2.1,
  §2.2, §2.3, §3.2, §3.4, §5, §9.A–§9.E, §13, §14) in full; the
  black-box outputs of the R probes `test/fixtures/auto_arima/r_probe/out.txt`,
  `out2.txt`, `out4.txt` (numbers only). The probe scripts `q*.R` in the same
  folder were **not** opened.
* Read (package, current state): `test/support/child_vm.dart` and the heavy
  test as changed by 3a57eec (to reuse the child-VM runner);
  `.dart_tool/native_assets.yaml` (to run study programs on a private copy
  of the native library).
* Nothing from the "Not opened" list above was opened.

## Addendum 3 (2026-10-06, spec revision 7)

* Read: specification §15 (revision 7) in full. Nothing else new; nothing
  from the "Not opened" list was opened.

## Addendum 4 (2026-10-06, spec revision 9 — disclosure by the architect)

Written by the architect (Nexus), not by the implementer.

* While preparing revision 9 of the specification, the **methodology owner
  (Arina)** opened the source of the R package `forecast`
  (`R/newarima2.R`) to check the neighbourhood order and the root margin.
  This violated the clean-room rule for the specification side. The
  implementer did not see that source and did not read the revision-9
  section (§17) that describes it. **Correction (same day, after Arina's
  review):** the revision-9 markers in §1 and §5 of the specification, which
  the implementer did read, quote one line of that source (the root-margin
  condition). The quoted behaviour is the one confirmed as a black box
  below; no other content from the source reached the implementer. Those
  markers are to be removed from the specification.
* The architect forwarded to the implementer only the two rules that were
  independently re-derived **as a black box**, with the data stored in
  `test/fixtures/auto_arima/r_probe/root_margin/`:
  1. the order in which the 17 neighbours are visited (seasonal block first,
     then non-seasonal, then the constant toggle) — read off ten stepwise
     traces of `auto.arima(trace = TRUE)` (`traces.txt`);
  2. the root margin 1.01 — every one of the 200 traced models was refitted
     with `Arima(method = "CSS-ML")` and the minimum root modulus of each
     polynomial recorded (`all_roots.csv`): all models R accepted have
     min |root| ≥ 1.0103, all models R rejected for roots have
     min |root| ≤ 1.0082, so the threshold lies in (1.0082, 1.0103).
* Details that Arina learned from the source and that could **not** be
  verified from outputs alone (internal model counters, the exact form of the
  coefficient-variance check, tail trimming of tiny coefficients) were **not**
  forwarded to the implementer and are not part of the specification given
  to them. Where our implementation needs such a rule it keeps its own
  `[наше]` decision (§5: `rejectedDegenerate`).
* The R scripts in `root_margin/` were written by the architect and were not
  given to the implementer; the implementer reads only `traces.txt` and
  `all_roots.csv` (numbers).

## Addendum 4 (2026-10-06, spec revisions 8–9) — provenance flag

* Read: specification §16 (revision 8) in full; from revision 9 only the
  journal line, the §3.2 and §5 markers, and the coordinator's restatement
  of the two changes (neighbour order; `rootMargin` 1.01, seasonal factor
  against 1.01^s). **§17 was not read**: the revision-9 journal line says it
  draws on the source of `forecast` 8.23 (`newarima2.R`), which this process
  forbids the implementer to see.
* Both revision-9 changes are implemented as **black-box facts**: the
  neighbour order is visible in R's `trace = TRUE` output
  (`test/fixtures/auto_arima/r_probe/root_margin/traces.txt`, read), and
  the margin is bracketed by refits of the traced models
  (`all_roots.csv`, read: accepted min|root| ≥ 1.0103, rejected
  1.0017–1.0082; the seasonal factors of R's rejected airline models have
  roots 1.065–1.088 in z^12, below 1.01^12 = 1.127, consistent with the
  seasonal margin 1.01^s). `roots_all.R` was not opened.
* Whether a specification that cites the GPL source still meets §12 is a
  process question for the owner; the implementation does not depend on any
  content beyond the two black-box facts above.
