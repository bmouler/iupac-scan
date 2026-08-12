"""Deterministic Shift-And versus nested-loop benchmark."""

from __future__ import annotations

import argparse
import json
import random
import statistics
import sys
import time
from collections.abc import Callable
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from iupac_scan import compile_motif, scan_sequence  # noqa: E402
from iupac_scan.core import IUPAC_MASKS  # noqa: E402

MOTIF = "BDHVNRYKMSWBDHVNRYKMSWBDHVNRYKMS"
SEED = 20260812


def naive_scan(sequence: str, motif: str) -> list[int]:
    """Return compatible starts using a correct expansion-free nested loop."""
    sequence_masks = [IUPAC_MASKS[symbol] for symbol in sequence]
    motif_masks = [IUPAC_MASKS[symbol] for symbol in motif]
    return [
        start
        for start in range(len(sequence_masks) - len(motif_masks) + 1)
        if all(
            sequence_masks[start + offset] & motif_mask
            for offset, motif_mask in enumerate(motif_masks)
        )
    ]


def _time(action: Callable[[], object], repeats: int) -> list[float]:
    samples = []
    for _ in range(repeats):
        started = time.perf_counter()
        action()
        samples.append(time.perf_counter() - started)
    return samples


def run(length: int, repeats: int) -> dict[str, int | float | str]:
    """Run and validate one seeded benchmark workload."""
    random_source = random.Random(SEED)
    sequence = "".join(random_source.choices("ACGT", k=length))
    insertion = length // 2
    concrete_motif = "".join(
        ("A" if IUPAC_MASKS[symbol] & 1 else "C" if IUPAC_MASKS[symbol] & 2 else "G")
        for symbol in MOTIF
    )
    sequence = (
        sequence[:insertion] + concrete_motif + sequence[insertion + len(MOTIF) :]
    )
    compiled = compile_motif(MOTIF)
    expected = naive_scan(sequence, MOTIF)
    observed = [match.start for match in scan_sequence(sequence, compiled)]
    if observed != expected:
        raise RuntimeError("benchmark implementations disagree")
    shift_samples = _time(lambda: list(scan_sequence(sequence, compiled)), repeats)
    naive_samples = _time(lambda: naive_scan(sequence, MOTIF), repeats)
    shift_median = statistics.median(shift_samples)
    naive_median = statistics.median(naive_samples)
    return {
        "seed": SEED,
        "sequence_length": length,
        "motif": MOTIF,
        "motif_length": len(MOTIF),
        "matches": len(expected),
        "repeats": repeats,
        "shift_and_median_seconds": shift_median,
        "naive_median_seconds": naive_median,
        "speedup": naive_median / shift_median,
        "result": "identical coordinates",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, default=200_000)
    parser.add_argument("--repeats", type=int, default=7)
    args = parser.parse_args()
    if args.length < len(MOTIF) or args.repeats < 1:
        parser.error("length must cover the motif and repeats must be positive")
    print(json.dumps(run(args.length, args.repeats), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
