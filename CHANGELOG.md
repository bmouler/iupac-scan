# Changelog

## [Unreleased]

## [1.0.1] - 2026-08-15

- Replaced per-base compatibility accumulation with a compiled 16-entry mask table while preserving IUPAC, overlap, and strand-order semantics.
- Added a deterministic end-to-end FASTA parse and both-strand scan benchmark with an independent nested-loop oracle.


## [1.0.0] - 2026-08-12

First stable release.

- Bit-parallel scanning of degenerate IUPAC DNA motifs.
- Added a deterministic property-based suite covering the scanner against a naive oracle, reverse-complement symmetry, and all-window matching.
- Reached a 98.59% mutation score: 491 of 498 mutants killed, with seven behavior-equivalent survivors documented.
- Adopted strict mypy checking and shipped inline typing metadata.
- Broadened CI to Linux and macOS across Python 3.11–3.13.
