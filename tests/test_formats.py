from __future__ import annotations

import pytest

from iupac_scan.formats import Record, parse_fasta, parse_fastq, parse_records


def test_fasta_streams_multiline_records_and_ignores_blank_lines() -> None:
    lines = iter([">one description\n", "AC\n", "\n", "GT\n", ">two\n", "NN\n"])
    records = parse_fasta(lines)
    assert next(records) == Record("one", "ACGT")
    assert next(records) == Record("two", "NN")
    with pytest.raises(StopIteration):
        next(records)


def test_fastq_streams_records() -> None:
    lines = ["\n", "@one description\n", "ACGT\n", "+\n", "!!!!\n"]
    assert list(parse_fastq(lines)) == [Record("one", "ACGT")]


@pytest.mark.parametrize(
    ("parser", "lines", "message"),
    [
        (parse_fasta, [], "FASTA input contains no records"),
        (parse_fasta, ["AC\n"], "expected FASTA header on line 1"),
        (parse_fasta, [">\n", "AC\n"], "missing record ID on line 1"),
        (parse_fasta, [">a\n", ">b\n", "AC\n"], "record 'a' has no sequence"),
        (parse_fasta, [">a\n"], "record 'a' has no sequence"),
        (parse_fasta, [">=formula\n", "AC\n"], "unsafe formula prefix"),
        (parse_fastq, [], "FASTQ input contains no records"),
        (parse_fastq, ["bad\n"], "expected FASTQ header on line 1"),
        (parse_fastq, ["@\n", "A\n", "+\n", "!\n"], "missing record ID"),
        (parse_fastq, ["@a\n"], "missing sequence line"),
        (parse_fastq, ["@a\n", "A\n"], "missing separator line"),
        (parse_fastq, ["@a\n", "A\n", "+\n"], "missing quality line"),
        (parse_fastq, ["@a\n", "\n", "+\n", "\n"], "has no sequence"),
        (parse_fastq, ["@a\n", "A\n", "x\n", "!\n"], "expected FASTQ '+'"),
        (
            parse_fastq,
            ["@a\n", "AA\n", "+\n", "!\n"],
            "sequence length 2 but quality length 1",
        ),
    ],
)
def test_malformed_records_raise_useful_errors(parser, lines, message: str) -> None:
    with pytest.raises(ValueError, match=message.replace("+", "\\+")):
        list(parser(lines))


def test_auto_detection_and_explicit_formats() -> None:
    assert list(parse_records(["\n", ">a\n", "A\n"])) == [Record("a", "A")]
    assert list(parse_records(["@a\n", "A\n", "+\n", "!\n"])) == [Record("a", "A")]
    assert list(parse_records([">a\n", "A\n"], "fasta")) == [Record("a", "A")]
    assert list(parse_records(["@a\n", "A\n", "+\n", "!\n"], "fastq")) == [
        Record("a", "A")
    ]


@pytest.mark.parametrize(
    ("lines", "file_format", "message"),
    [
        ([], "auto", "input contains no records"),
        (["\n", "garbage\n"], "auto", "cannot detect input format"),
        ([], "sam", "unsupported input format 'sam'"),
    ],
)
def test_record_dispatch_errors(lines, file_format: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        list(parse_records(lines, file_format))
