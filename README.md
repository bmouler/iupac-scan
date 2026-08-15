# iupac-scan

[![CI](https://github.com/bmouler/iupac-scan/actions/workflows/ci.yml/badge.svg)](https://github.com/bmouler/iupac-scan/actions/workflows/ci.yml)
![Coverage](https://img.shields.io/badge/coverage-100%25-brightgreen)
![Types](https://img.shields.io/badge/types-mypy%20strict-blue)
![Mutation](https://img.shields.io/badge/mutation-98%25%20killed-brightgreen)
![Python](https://img.shields.io/badge/python-3.11%2B-blue)
[![License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

`iupac-scan` streams FASTA and conventional four-line FASTQ records and finds overlapping degenerate DNA motifs on one or both strands. It has no runtime dependencies and makes no network requests.

## Install

```console
python -m pip install iupac-scan
```

For development:

```console
python -m pip install -e '.[dev]'
```

## Quickstart

```console
$ iupac-scan 'ATN' reads.fasta --both-strands
record_id	start	end	strand	matched_sequence
read1	0	3	+	ATG
read1	4	7	-	CAT
```

Coordinates are zero-based, half-open `[start, end)`, and always refer to the input record's orientation. Output defaults to TSV; `--output-format jsonl` emits one JSON object per match. Use `-` or omit the input path to read stdin, `-o PATH` to write a file, and `--input-format fasta|fastq` to override detection.

Python API:

```python
from iupac_scan import compile_motif, scan_sequence

motif = compile_motif("RYN")
for match in scan_sequence("ACGTN", motif, both_strands=True):
    print(match.start, match.end, match.strand, match.sequence)
```

## Exact ambiguous-symbol semantics

Every IUPAC symbol is a set of canonical DNA bases: `A=A`, `C=C`, `G=G`, `T=T`, `U=T`, `R=AG`, `Y=CT`, `S=CG`, `W=AT`, `K=GT`, `M=AC`, `B=CGT`, `D=AGT`, `H=ACT`, `V=ACG`, and `N=ACGT`. Matching is case-insensitive. A motif position and sequence position match exactly when their sets intersect, so an ambiguous input symbol can match an ambiguous motif symbol. `U` is treated as `T`; this is still a DNA scanner, not an RNA reverse-complement tool.

Matches may overlap. With `--both-strands`, the motif and its reverse complement are scanned independently. A palindromic motif therefore produces both `+` and `-` output at the same coordinates. The matched sequence is reported exactly in the forward input orientation, normalized to uppercase.

## Algorithm

```mermaid
flowchart LR; M[IUPAC motif] --> C[compile: 4 bitmask vectors]; S[FASTA/FASTQ records] --> X[Shift-And state, 1 int]; C --> X; X --> B{both strands?}; B -->|reverse complement| X; X --> O[overlapping matches, zero-based]
```

Compilation converts each motif character to a four-bit `A/C/G/T` mask, then transposes the motif into four Python-integer bit vectors. For each strand, scanning combines those vectors into a 16-entry table indexed by the sequence symbol's four-bit mask, then applies one Shift-And state update and one table lookup per symbol. Motif degeneracy does not expand into concrete motifs or alter asymptotic scan cost; scan memory remains constant apart from emitted matches.

## Reproducible local evidence

`PYTHONPATH=src python benchmarks/benchmark.py` compiles `ATNACAT`, parses 24 wrapped FASTA
records totaling 192,000 bases, scans overlapping matches on both strands, and materializes all
ordered match tuples. Before timing, it requires exact equality with an independent
expansion-free nested-loop oracle.

On an Apple M3 Max with CPython 3.11.12 on 2026-08-15, 13 samples after three warmups measured
frozen baseline `d9d6db50b342` at **89.047 ms** median and this implementation at
**30.541 ms**, a **2.916x speedup**. Both runs found 237 matches and produced SHA-256
`3700c5424eb93982b201ad55d39c98b191480afe6ec71c5fb5ff25b73c2401a1`. Fixture
generation and interpreter startup are excluded; motif compilation, parsing, both scans, and
output materialization are included. These are local timings; rerun with `PYTHONPATH` pointed
at the desired source worktree.

## Verification

### Mutation testing

From the repository root, reproduce the mutation run with:

```console
source .venv/bin/activate
mutmut run
mutmut results
```

The run generated **498 mutants**, of which **491 were killed (98.59%)**. The seven remaining survivors were individually reviewed and are behavior-equivalent, not missed mutants. There were **zero suspicious mutants and zero timeouts**.

| Behavior-equivalent rationale | Count |
| --- | ---: |
| Bitmask bits discarded by subsequent shifts | 2 |
| Identifier split `maxsplit` variant with identical output | 1 |
| Initialization sentinel variants with identical observable behavior | 2 |
| Duplicate characters in `rstrip` character sets | 2 |

## Limitations

- FASTQ parsing intentionally accepts the widespread four-line form only; wrapped sequence or quality lines are rejected.
- FASTA sequence lines may wrap, but whitespace inside sequence lines is not removed and is reported as an invalid DNA symbol.
- Record IDs beginning with `=`, `+`, `-`, or `@` are rejected so default TSV output is safe from spreadsheet formula interpretation.
- Quality scores are length-validated but otherwise ignored.
- Inputs are decoded as UTF-8 and symbols outside the documented alphabet are errors.
- Very long motifs remain correct, but Python big-integer work grows with motif length.
