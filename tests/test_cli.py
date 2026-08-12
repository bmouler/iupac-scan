from __future__ import annotations

import io
import json
import subprocess
import sys
from pathlib import Path

from iupac_scan import cli

EXPECTED_HELP = (
    "usage: iupac-scan [-h] [--both-strands] "
    "[--input-format {auto,fasta,fastq}]\n"
    """                  [--output-format {tsv,jsonl}] [-o OUTPUT]
                  motif [input]

Scan FASTA or FASTQ records for an IUPAC DNA motif.

positional arguments:
  motif                 non-empty motif using the IUPAC DNA alphabet
  input                 input FASTA/FASTQ path (default: stdin)

options:
  -h, --help            show this help message and exit
  --both-strands        also scan the reverse complement
  --input-format {auto,fasta,fastq}
  --output-format {tsv,jsonl}
  -o OUTPUT, --output OUTPUT
                        output path (default: stdout)
"""
)


def test_parser_help_documents_the_complete_cli_contract() -> None:
    assert cli._parser().format_help() == EXPECTED_HELP


def test_jsonl_is_compact() -> None:
    output = io.StringIO()
    cli._write_match(output, "jsonl", "r", 1, 2, "+", "A")
    assert output.getvalue() == (
        '{"record_id":"r","start":1,"end":2,"strand":"+","matched_sequence":"A"}\n'
    )


def test_explicit_input_format_is_honored(monkeypatch, capsys) -> None:
    monkeypatch.setattr(sys, "stdin", io.StringIO("@r\nA\n+\n!\n"))
    assert cli.main(["A", "--input-format", "fasta"]) == 2
    assert "expected FASTA header on line 1" in capsys.readouterr().err


def test_file_io_is_utf8_under_an_ascii_locale(tmp_path: Path) -> None:
    source = tmp_path / "reads.fa"
    output = tmp_path / "matches.tsv"
    source.write_bytes(">réad\nA\n".encode())
    command = Path(sys.executable).with_name("iupac-scan")
    completed = subprocess.run(
        [str(command), "A", str(source), "--output", str(output)],
        capture_output=True,
        text=True,
        check=False,
        env={
            "PATH": str(command.parent),
            "PYTHONUTF8": "0",
            "LC_ALL": "C",
            "LANG": "C",
        },
    )
    assert completed.returncode == 0, completed.stderr
    assert output.read_text(encoding="utf-8").splitlines()[1] == "réad\t0\t1\t+\tA"


def test_open_helpers_pin_utf8_and_newline_policy(monkeypatch, tmp_path: Path) -> None:
    observed: list[tuple[str, tuple[object, ...], dict[str, object]]] = []
    original_open = Path.open

    def recording_open(path: Path, *args: object, **kwargs: object):
        observed.append((str(path), args, kwargs))
        return original_open(path, *args, **kwargs)

    input_path = tmp_path / "in.fa"
    output_path = tmp_path / "out.tsv"
    input_path.write_text(">a\nA\n", encoding="utf-8")
    monkeypatch.setattr(Path, "open", recording_open)
    with cli._open_input(str(input_path)):
        pass
    with cli._open_output(str(output_path)):
        pass
    assert observed == [
        (str(input_path), (), {"encoding": "utf-8"}),
        (str(output_path), ("w",), {"encoding": "utf-8", "newline": ""}),
    ]


def test_main_reads_stdin_and_writes_tsv(monkeypatch, capsys) -> None:
    monkeypatch.setattr(sys, "stdin", io.StringIO(">r1\nAAAA\n"))
    assert cli.main(["AA"]) == 0
    captured = capsys.readouterr()
    assert captured.err == ""
    assert captured.out.splitlines() == [
        "record_id\tstart\tend\tstrand\tmatched_sequence",
        "r1\t0\t2\t+\tAA",
        "r1\t1\t3\t+\tAA",
        "r1\t2\t4\t+\tAA",
    ]


def test_main_reads_fastq_and_writes_jsonl_file(tmp_path: Path) -> None:
    source = tmp_path / "reads.fastq"
    output = tmp_path / "matches.jsonl"
    source.write_text("@read/1 note\nATGACAT\n+\nIIIIIII\n", encoding="utf-8")
    assert (
        cli.main(
            [
                "ATG",
                str(source),
                "--both-strands",
                "--output-format",
                "jsonl",
                "--output",
                str(output),
            ]
        )
        == 0
    )
    rows = [
        json.loads(line) for line in output.read_text(encoding="utf-8").splitlines()
    ]
    assert rows == [
        {
            "record_id": "read/1",
            "start": 0,
            "end": 3,
            "strand": "+",
            "matched_sequence": "ATG",
        },
        {
            "record_id": "read/1",
            "start": 4,
            "end": 7,
            "strand": "-",
            "matched_sequence": "CAT",
        },
    ]


def test_main_reports_data_and_io_errors(tmp_path: Path, capsys) -> None:
    bad = tmp_path / "bad.fa"
    bad.write_text(">x\nAX\n", encoding="utf-8")
    assert cli.main(["A", str(bad)]) == 2
    assert "invalid sequence symbol 'X' at offset 1" in capsys.readouterr().err
    assert cli.main(["A", str(tmp_path / "missing.fa")]) == 2
    assert "error:" in capsys.readouterr().err


def test_main_reports_invalid_motif_before_opening_input(capsys) -> None:
    assert cli.main(["X", "missing.fa"]) == 2
    assert "invalid motif symbol 'X'" in capsys.readouterr().err


def test_main_handles_broken_pipe(monkeypatch) -> None:
    class BrokenOutput(io.StringIO):
        def write(self, value: str) -> int:
            raise BrokenPipeError

    monkeypatch.setattr(sys, "stdout", BrokenOutput())
    monkeypatch.setattr(sys, "stdin", io.StringIO(">x\nA\n"))
    assert cli.main(["A"]) == 0


def test_open_helpers_use_standard_streams(monkeypatch) -> None:
    input_stream = io.StringIO("")
    output_stream = io.StringIO()
    monkeypatch.setattr(sys, "stdin", input_stream)
    monkeypatch.setattr(sys, "stdout", output_stream)
    with cli._open_input("-") as observed_input:
        assert observed_input is input_stream
    with cli._open_output("-") as observed_output:
        assert observed_output is output_stream


def test_module_end_to_end_from_clean_directory(tmp_path: Path) -> None:
    source = tmp_path / "input.fa"
    source.write_text(">alpha description\nACGTT\n", encoding="utf-8")
    command = Path(sys.executable).with_name("iupac-scan")
    completed = subprocess.run(
        [
            str(command),
            "RY",
            str(source),
            "--output-format",
            "jsonl",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert [json.loads(line) for line in completed.stdout.splitlines()] == [
        {
            "record_id": "alpha",
            "start": 0,
            "end": 2,
            "strand": "+",
            "matched_sequence": "AC",
        },
        {
            "record_id": "alpha",
            "start": 2,
            "end": 4,
            "strand": "+",
            "matched_sequence": "GT",
        },
    ]
