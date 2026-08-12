"""Bit-parallel IUPAC DNA motif scanning."""

from .core import CompiledMotif, Match, compile_motif, scan_sequence

__all__ = ["CompiledMotif", "Match", "compile_motif", "scan_sequence"]
