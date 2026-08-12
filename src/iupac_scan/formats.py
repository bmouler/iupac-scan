"""Streaming FASTA and FASTQ parsing."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from itertools import chain


@dataclass(frozen=True, slots=True)
class Record:
    """A sequence record."""

    identifier: str
    sequence: str


_SPREADSHEET_FORMULA_PREFIXES = frozenset("=+-@")


def _identifier(header: str, line_number: int) -> str:
    identifier = header.strip().split(maxsplit=1)[0] if header.strip() else ""
    if not identifier:
        raise ValueError(f"missing record ID on line {line_number}")
    if identifier[0] in _SPREADSHEET_FORMULA_PREFIXES:
        raise ValueError(f"unsafe formula prefix in record ID on line {line_number}")
    return identifier


def parse_fasta(lines: Iterable[str]) -> Iterator[Record]:
    """Yield FASTA records without loading the input file."""
    identifier: str | None = None
    sequence_parts: list[str] = []
    for line_number, raw_line in enumerate(lines, 1):
        line = raw_line.rstrip("\r\n")
        if line.startswith(">"):
            if identifier is not None:
                if not sequence_parts:
                    raise ValueError(f"record {identifier!r} has no sequence")
                yield Record(identifier, "".join(sequence_parts))
            identifier = _identifier(line[1:], line_number)
            sequence_parts = []
        elif not line:
            continue
        elif identifier is None:
            raise ValueError(f"expected FASTA header on line {line_number}")
        else:
            sequence_parts.append(line.strip())
    if identifier is None:
        raise ValueError("FASTA input contains no records")
    if not sequence_parts:
        raise ValueError(f"record {identifier!r} has no sequence")
    yield Record(identifier, "".join(sequence_parts))


def parse_fastq(lines: Iterable[str]) -> Iterator[Record]:
    """Yield conventional four-line FASTQ records."""
    iterator = iter(enumerate(lines, 1))
    found = False
    while True:
        try:
            header_number, header_raw = next(iterator)
        except StopIteration:
            break
        header = header_raw.rstrip("\r\n")
        if not header and not found:
            continue
        found = True
        if not header.startswith("@"):
            raise ValueError(f"expected FASTQ header on line {header_number}")
        identifier = _identifier(header[1:], header_number)
        fields: list[tuple[int, str]] = []
        for label in ("sequence", "separator", "quality"):
            try:
                number, raw = next(iterator)
            except StopIteration as error:
                raise ValueError(
                    f"truncated FASTQ record {identifier!r}: missing {label} line"
                ) from error
            fields.append((number, raw.rstrip("\r\n")))
        sequence = fields[0][1]
        separator_number, separator = fields[1]
        quality = fields[2][1]
        if not sequence:
            raise ValueError(f"record {identifier!r} has no sequence")
        if not separator.startswith("+"):
            raise ValueError(f"expected FASTQ '+' separator on line {separator_number}")
        if len(quality) != len(sequence):
            raise ValueError(
                f"record {identifier!r} has sequence length {len(sequence)} "
                f"but quality length {len(quality)}"
            )
        yield Record(identifier, sequence)
    if not found:
        raise ValueError("FASTQ input contains no records")


def parse_records(lines: Iterable[str], file_format: str = "auto") -> Iterator[Record]:
    """Detect or apply an input format and stream records."""
    if file_format == "fasta":
        yield from parse_fasta(lines)
        return
    if file_format == "fastq":
        yield from parse_fastq(lines)
        return
    if file_format != "auto":
        raise ValueError(f"unsupported input format {file_format!r}")

    iterator = iter(lines)
    prefix: list[str] = []
    for line in iterator:
        prefix.append(line)
        stripped = line.lstrip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            yield from parse_fasta(chain(prefix, iterator))
            return
        if stripped.startswith("@"):
            yield from parse_fastq(chain(prefix, iterator))
            return
        raise ValueError("cannot detect input format: expected '>' or '@'")
    raise ValueError("input contains no records")
