# AGENT.md — RepeatNet working notes

Context and session history for AI agents working on this repo. Keep this file
up to date as work progresses.

## What this project is

**RepeatNet**: an ab initio centromeric satellite detection algorithm. Single
C program (`repeatnet.c`, ~2100 lines) that takes paired end reads, extracts
12-mer windows (k = `WS`, default 12), finds windows that occur on *both* ends
of the same clone/pair ("double-sided" windows), builds a co-occurrence graph
between windows across pairs, and extracts connected components (candidate
satellite families). See `README.md` for the biology and the full pipeline.

Downstream tools (external): PHRAP (assembly), TRF (consensus), Graphviz
(optional viz).

## Build

```
make            # gcc -O3 -Wall, links -lm  -> ./repeatnet
make debug      # -O1 -g -fsanitize=address,undefined  (used for verification)
make clean
```

`repeatnet --help` lists all options.

## Pipeline (typical run)

1. `repeatnet -i test`  → `test.h11.winlog`, `test.h11.names`, `test.h11.dump`
2. `repeatnet --loadwin test.h11.winlog -m -a 100 -c 100 --compare` → `.matrix`, `.viz`
3. `repeatnet --loadmatrix <.matrix> --loadnames test.h11.names --components`
   → `component-N.{matrix,viz,txt}`
4. `repeatnet --encode <kmer>` / `repeatnet --decode <file>` convert between
   k-mer strings and integer window codes.
5. `repeatnet --loadwin ... --loadnames ... --clones <ids>` → clone names for
   fetching original reads.

Note: the program uses `return 1` as its *success* exit code on most paths
(legacy convention), `return 0` for input errors. Don't treat exit 1 as failure.

## Input formats (current state)

- **FASTA and FASTQ**, auto-detected from the first character (`>` = FASTA,
  `@` = FASTQ). No flag needed.
- **Interleaved** (`-i/--input`): forward+reverse mates in one file, paired by
  matching clone name (FORWARD/REVERSE tag stripped by `end2clone`, then
  hashed). Names of the two mates must share the clone-id part.
- **Two files** (`-f/--forward`, `-r/--reverse`): mates paired by position
  (i-th record of each file), no name matching. Formats detected per file, so
  FASTA-forward + FASTQ-reverse is allowed. Output files are named after the
  forward file.

## Key functions / architecture

- `main` — getopt_long parsing, mode dispatch.
- `detectFormat(FILE*)` — peek first char → `FASTA_FMT`/`FASTQ_FMT`.
- `readSeqRecord(f, fmt, &name, &seq)` — generic single-record reader (FASTA and
  FASTQ). FASTQ: reads header, seq until `+` line, then discards *exactly*
  `strlen(seq)` quality chars — this is what makes `@`/`+` in quality safe.
  Handles multi-line seq/qual.
- `readSingleFasta(in, hashtype)` — **legacy** interleaved FASTA reader, left
  untouched to preserve byte-exact output. Char-by-char parser.
- `readInterleavedFastq(in, hashtype)` — interleaved FASTQ, mirrors
  readSingleFasta's hashing/pairing via readSeqRecord.
- `readPaired(fwd, rev)` — two-file mode.
- `hash()` — 15 selectable hash functions (types 0–14); default 11 (DEK).
- `encode_window2` / `decode_window` — 2-bit k-mer <-> int packing.
- `do_end` — extracts windows for one read end, records double-sided hits.
- `compare()` — builds the co-occurrence graph/matrix.
- `connComp()` — union-find-ish connected components (Gozde's code, `//gozde`).

## Work done in the sessions so far

### Session 1 — refactor for readability/efficiency + memory fixes
Established a byte-for-byte behavioral baseline first (synthetic paired FASTA,
captured md5s of winlog/names/dump/matrix/viz/components), then verified every
change against it. All outputs stayed identical.

Fixes:
- **Heap overflow** in `do_end`: window buffer was `malloc(WS)` but `getWindow`
  writes a `\0` at `window[WS]` → now `WS+1`.
- **`connComp` bug**: the "only v2 clustered" branch appended to
  `lastNodes[seq1Comp]` (stale) instead of `seq2Comp`.
- **Uninitialized `->next`** on newly malloc'd list/interaction nodes in
  `connComp` → caused an ASan SEGV on `--components` (present in original too);
  now every node is NULL-terminated on creation.
- **UB signed left-shift** in `decode_window` → shift on `unsigned`.
- **Uninitialized `names`** in `-t` hashtest path: `malloc`→`calloc`.
- Removed dead code (unreachable compare block in main, dead tails of
  `encode_window`/`revcomp_encoded`, `mds()` stub, unused locals), fixed a
  leaked `malloc`, added missing `fclose`s, fixed a `%ld`/`int` format bug.
- **Efficiency**: `compare()` was O(4^WS²) scanning all ~16.7M window slots
  squared; now collects the few active windows first and compares only those
  (i<j order preserved) → ~3× faster even on tiny data, scales with window
  space. Hoisted `strlen` out of loops in `hash()` and `end2clone`.

### Session 2 — two-file input + getopt overhaul
- Added two-file mode (`-f`/`-r`, positional pairing).
- Replaced the hand-rolled `if/strcmp` argv loop with `getopt_long`; added
  `--help`. **Breaking CLI change**: hash-type flag moved `-h` → `-H/--hash`
  (`-h` was ambiguous with help). Long options added for every stage.
- Verified interleaved FASTA output still byte-identical.

### Session 3 — FASTA/FASTQ autodetection
- Added `detectFormat`, generalized reader to `readSeqRecord`, added
  `readInterleavedFastq`, made `readPaired` detect per-file.
- Verified: FASTQ interleaved, FASTQ two-file, and mixed FASTA/FASTQ all produce
  window histograms identical to the FASTA baseline. Adversarial qualities
  (starting with `@`/`+`, containing `@ + ! I ~`) and multi-line wrapping parse
  correctly. ASan/UBSan clean on all paths.

## How to verify changes (recipe used)

1. `make debug` and run the affected mode; grep stderr for `ERROR|runtime error|SEGV`.
2. Generate synthetic paired reads where forward and the reverse-complement of
   the reverse read share a common motif (so windows are double-sided and
   co-occur across pairs). Seed fixed.
3. The **window→count histogram (`.dump`)** is independent of how pairs are
   indexed, so it must match across interleaved/two-file and FASTA/FASTQ for the
   same underlying reads — the primary cross-check.
4. For pure refactors, compare full output md5s against a saved baseline.

## Known remaining issues / not done

- Pre-existing compiler warnings under `-Wall` (all in legacy code, none in new
  code): `maybe-uninitialized` false positives (`index`, `edgelog`, `nfile`,
  `minhash`, `mincolhash`), and a `format-overflow`/`-Wrestrict` in `trimwins`
  (a `sprintf(fname2, "...", fname2, ...)` self-overlap). Left untouched to
  avoid risky edits to legacy logic. Candidates for a future cleanup.
- `-t` hashtest path is debug-only and FASTA/token oriented; not adapted to FASTQ.
- Remaining `fscanf`/`fread` unused-result warnings (harmless).
- No automated test suite in-repo; verification is manual (see recipe above).
- Nothing has been committed — all changes are in the working tree.

## Conventions when editing

- Match the surrounding (legacy) style in the file; keep new code clean.
- Preserve FASTA byte-exact output — always re-baseline before/after.
- Don't "fix" the `return 1 == success` convention without checking all callers.
