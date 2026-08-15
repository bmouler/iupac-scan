"""Deterministic end-to-end benchmark for the documented FASTA scan path."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import statistics
import time
from collections.abc import Callable

from iupac_scan import compile_motif, scan_sequence
from iupac_scan.core import IUPAC_MASKS
from iupac_scan.formats import Record, parse_records

MOTIF = "ATNACAT"
SEED = 20260815


def _fixture(
    records: int, sequence_length: int
) -> tuple[tuple[str, ...], tuple[Record, ...]]:
    """Build fixed wrapped FASTA input outside the timed region."""
    random_source = random.Random(SEED)
    lines: list[str] = []
    expected_records: list[Record] = []
    plus_overlaps = "ATGACATGACATGACAT"
    minus_overlaps = "ATGTCATGTCATGTCAT"
    for record_index in range(records):
        sequence = list(random_source.choices("ACGT", k=sequence_length))
        for position in range(499 + record_index % 17, sequence_length, 997):
            sequence[position] = "NRY"[(position + record_index) % 3]
        plus_start = sequence_length // 4 + record_index % 23
        minus_start = 3 * sequence_length // 5 + record_index % 29
        sequence[plus_start : plus_start + len(plus_overlaps)] = plus_overlaps
        sequence[minus_start : minus_start + len(minus_overlaps)] = minus_overlaps
        text = "".join(sequence)
        identifier = f"record-{record_index:03d}"
        expected_records.append(Record(identifier, text))
        lines.append(f">{identifier} deterministic benchmark record\n")
        lines.extend(
            f"{text[offset : offset + 80]}\n" for offset in range(0, len(text), 80)
        )
    return tuple(lines), tuple(expected_records)


def _complement(mask: int) -> int:
    return ((mask & 1) << 3) | ((mask & 2) << 1) | ((mask & 4) >> 1) | ((mask & 8) >> 3)


def _reference_matches(records: tuple[Record, ...]) -> tuple[tuple[object, ...], ...]:
    """Compute an expansion-free nested-loop oracle for exact match objects."""
    motif_masks = tuple(IUPAC_MASKS[symbol] for symbol in MOTIF)
    strand_masks = (
        motif_masks,
        tuple(_complement(mask) for mask in reversed(motif_masks)),
    )
    matches: list[tuple[object, ...]] = []
    for record in records:
        normalized = record.sequence.upper()
        sequence_masks = tuple(IUPAC_MASKS[symbol] for symbol in normalized)
        for strand, masks in zip(("+", "-"), strand_masks, strict=True):
            width = len(masks)
            for start in range(len(sequence_masks) - width + 1):
                if all(
                    sequence_masks[start + offset] & mask
                    for offset, mask in enumerate(masks)
                ):
                    matches.append(
                        (
                            record.identifier,
                            start,
                            start + width,
                            strand,
                            normalized[start : start + width],
                        )
                    )
    return tuple(matches)


def _scan(lines: tuple[str, ...]) -> tuple[tuple[object, ...], ...]:
    """Run the complete public compile, FASTA parse, scan, and materialize path."""
    motif = compile_motif(MOTIF)
    return tuple(
        (record.identifier, match.start, match.end, match.strand, match.sequence)
        for record in parse_records(lines)
        for match in scan_sequence(record.sequence, motif, both_strands=True)
    )


def _time(action: Callable[[], object], warmups: int, samples: int) -> list[float]:
    for _ in range(warmups):
        action()
    timings = []
    for _ in range(samples):
        started = time.perf_counter()
        action()
        timings.append(time.perf_counter() - started)
    return timings


def run(
    records: int, sequence_length: int, warmups: int, samples: int
) -> dict[str, object]:
    """Validate and measure one seeded end-to-end workload."""
    lines, expected_records = _fixture(records, sequence_length)
    parsed_records = tuple(parse_records(lines))
    if parsed_records != expected_records:
        raise RuntimeError("FASTA parsing differs from the exact fixture records")
    expected = _reference_matches(expected_records)
    observed = _scan(lines)
    if observed != expected:
        raise RuntimeError("public scan differs from the exact nested-loop oracle")
    timings = _time(lambda: _scan(lines), warmups, samples)
    encoded = json.dumps(expected, separators=(",", ":"), ensure_ascii=True).encode()
    return {
        "seed": SEED,
        "records": records,
        "sequence_length_per_record": sequence_length,
        "total_bases": records * sequence_length,
        "fasta_lines": len(lines),
        "motif": MOTIF,
        "motif_length": len(MOTIF),
        "both_strands": True,
        "matches": len(expected),
        "overlapping_matches": any(
            first[0] == second[0]
            and first[3] == second[3]
            and int(second[1]) < int(first[2])
            for first, second in zip(expected, expected[1:], strict=False)
        ),
        "warmups": warmups,
        "samples": samples,
        "median_seconds": statistics.median(timings),
        "min_seconds": min(timings),
        "max_seconds": max(timings),
        "timings_seconds": timings,
        "equivalence": (
            "exact parsed records and ordered match tuples equal nested-loop oracle"
        ),
        "match_sha256": hashlib.sha256(encoded).hexdigest(),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", type=int, default=24)
    parser.add_argument("--sequence-length", type=int, default=8_000)
    parser.add_argument("--warmups", type=int, default=3)
    parser.add_argument("--samples", type=int, default=13)
    args = parser.parse_args()
    if (
        args.records < 1
        or args.sequence_length < 100
        or args.warmups < 0
        or args.samples < 11
    ):
        parser.error(
            "records must be positive, sequence length at least 100, warmups "
            "nonnegative, and samples at least 11"
        )
    print(
        json.dumps(
            run(args.records, args.sequence_length, args.warmups, args.samples),
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
