from __future__ import annotations

import io
import json
import subprocess
import sys
from pathlib import Path

from iupac_scan import cli


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
