from __future__ import annotations

import runpy
import sys

import pytest


def test_module_entrypoint_exits_with_main_status(monkeypatch, tmp_path) -> None:
    source = tmp_path / "input.fa"
    source.write_text(">x\nA\n", encoding="utf-8")
    monkeypatch.setattr(sys, "argv", ["iupac-scan", "A", str(source)])
    with pytest.raises(SystemExit, match="0"):
        runpy.run_module("iupac_scan", run_name="__main__")
