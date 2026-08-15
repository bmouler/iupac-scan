"""Compile and scan degenerate DNA motifs with Shift-And."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from typing import Final, Literal, TypeAlias

IUPAC_MASKS: Final[Mapping[str, int]] = {
    "A": 0b0001,
    "C": 0b0010,
    "G": 0b0100,
    "T": 0b1000,
    "U": 0b1000,
    "R": 0b0101,
    "Y": 0b1010,
    "S": 0b0110,
    "W": 0b1001,
    "K": 0b1100,
    "M": 0b0011,
    "B": 0b1110,
    "D": 0b1101,
    "H": 0b1011,
    "V": 0b0111,
    "N": 0b1111,
}

_MASK_NAMES: Final = " ACMGRSVTWYHKDBN"

Strand: TypeAlias = Literal["+", "-"]


@dataclass(frozen=True, slots=True)
class CompiledMotif:
    """A motif represented by compatibility bit-vectors for Shift-And."""

    motif: str
    position_masks: tuple[int, ...]
    shift_masks: tuple[int, int, int, int]
    accept_bit: int

    def reverse_complement(self) -> CompiledMotif:
        """Return the compiled reverse-complement motif."""
        masks = tuple(_complement(mask) for mask in reversed(self.position_masks))
        return _compile_masks(_mask_name(mask) for mask in masks)


@dataclass(frozen=True, slots=True)
class Match:
    """One motif occurrence in sequence coordinates."""

    start: int
    end: int
    strand: Strand
    sequence: str


def _invalid_symbol_error(kind: str, text: str, index: int) -> ValueError:
    symbol = text[index]
    return ValueError(f"invalid {kind} symbol {symbol!r} at offset {index}")


def _validate_iupac(text: str, kind: str) -> str:
    normalized = text.upper()
    if not normalized:
        raise ValueError(f"{kind} must not be empty")
    for index, symbol in enumerate(normalized):
        if symbol not in IUPAC_MASKS:
            raise _invalid_symbol_error(kind, text, index)
    return normalized


def _complement(mask: int) -> int:
    return (
        ((mask & 0b0001) << 3)
        | ((mask & 0b0010) << 1)
        | ((mask & 0b0100) >> 1)
        | ((mask & 0b1000) >> 3)
    )


def _mask_name(mask: int) -> str:
    return _MASK_NAMES[mask]


def _compile_masks(symbols: Iterable[str]) -> CompiledMotif:
    motif = "".join(symbols)
    position_masks = tuple(IUPAC_MASKS[symbol] for symbol in motif)
    compatibility: list[int] = []
    for base_bit in (1, 2, 4, 8):
        vector = 0
        for position, motif_mask in enumerate(position_masks):
            if motif_mask & base_bit:
                vector |= 1 << position
        compatibility.append(vector)
    return CompiledMotif(
        motif,
        position_masks,
        (
            compatibility[0],
            compatibility[1],
            compatibility[2],
            compatibility[3],
        ),
        1 << (len(motif) - 1),
    )


def compile_motif(motif: str) -> CompiledMotif:
    """Validate and compile a non-empty IUPAC motif."""
    return _compile_masks(iter(_validate_iupac(motif, "motif")))


def _scan_one(
    sequence: str, compiled: CompiledMotif, strand: Strand
) -> Iterator[Match]:
    state = 0
    length = len(compiled.motif)
    a, c, g, t = compiled.shift_masks
    compatible_by_mask = (
        0,
        a,
        c,
        a | c,
        g,
        a | g,
        c | g,
        a | c | g,
        t,
        a | t,
        c | t,
        a | c | t,
        g | t,
        a | g | t,
        c | g | t,
        a | c | g | t,
    )
    masks = IUPAC_MASKS
    accept_bit = compiled.accept_bit
    for end, symbol in enumerate(sequence):
        state = ((state << 1) | 1) & compatible_by_mask[masks[symbol]]
        if state & accept_bit:
            start = end - length + 1
            yield Match(start, end + 1, strand, sequence[start : end + 1])


def scan_sequence(
    sequence: str,
    motif: str | CompiledMotif,
    *,
    both_strands: bool = False,
) -> Iterator[Match]:
    """Yield overlapping matches, optionally on both DNA strands.

    Two IUPAC symbols are compatible when their represented base sets intersect.
    Coordinates always refer to the supplied sequence. Palindromic motifs produce
    separate ``+`` and ``-`` matches when both strands are requested.
    """
    normalized = _validate_iupac(sequence, "sequence")
    compiled = compile_motif(motif) if isinstance(motif, str) else motif
    yield from _scan_one(normalized, compiled, "+")
    if both_strands:
        yield from _scan_one(normalized, compiled.reverse_complement(), "-")
