"""Command-line interface for iupac-scan."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from contextlib import AbstractContextManager, nullcontext
from pathlib import Path
from typing import Final, TextIO

from .core import compile_motif, scan_sequence
from .formats import parse_records


class _Arguments(argparse.Namespace):
    motif: str
    input: str
    both_strands: bool
    input_format: str
    output_format: str
    output: str


_FIELDS: Final = ("record_id", "start", "end", "strand", "matched_sequence")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="iupac-scan",
        description="Scan FASTA or FASTQ records for an IUPAC DNA motif.",
    )
    parser.add_argument("motif", help="non-empty motif using the IUPAC DNA alphabet")
    parser.add_argument(
        "input", nargs="?", default="-", help="input FASTA/FASTQ path (default: stdin)"
    )
    parser.add_argument(
        "--both-strands", action="store_true", help="also scan the reverse complement"
    )
    parser.add_argument(
        "--input-format", choices=("auto", "fasta", "fastq"), default="auto"
    )
    parser.add_argument("--output-format", choices=("tsv", "jsonl"), default="tsv")
    parser.add_argument(
        "-o", "--output", default="-", help="output path (default: stdout)"
    )
    return parser


def _open_input(path: str) -> AbstractContextManager[TextIO]:
    if path == "-":
        return nullcontext(sys.stdin)
    return Path(path).open(encoding="utf-8")


def _open_output(path: str) -> AbstractContextManager[TextIO]:
    if path == "-":
        return nullcontext(sys.stdout)
    return Path(path).open("w", encoding="utf-8", newline="")


def _write_match(
    output: TextIO,
    output_format: str,
    record_id: str,
    start: int,
    end: int,
    strand: str,
    sequence: str,
) -> None:
    values = (record_id, start, end, strand, sequence)
    if output_format == "tsv":
        output.write("\t".join(map(str, values)) + "\n")
    else:
        output.write(
            json.dumps(dict(zip(_FIELDS, values, strict=False)), separators=(",", ":"))
            + "\n"
        )


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line scanner and return a process exit status."""
    args = _Arguments()
    _parser().parse_args(argv, namespace=args)
    try:
        motif = compile_motif(args.motif)
        with _open_input(args.input) as input_file, _open_output(args.output) as output:
            if args.output_format == "tsv":
                output.write("\t".join(_FIELDS) + "\n")
            for record in parse_records(input_file, args.input_format):
                for match in scan_sequence(
                    record.sequence, motif, both_strands=args.both_strands
                ):
                    _write_match(
                        output,
                        args.output_format,
                        record.identifier,
                        match.start,
                        match.end,
                        match.strand,
                        match.sequence,
                    )
    except BrokenPipeError:
        return 0
    except (OSError, UnicodeError, ValueError) as error:
        print(f"iupac-scan: error: {error}", file=sys.stderr)
        return 2
    return 0
