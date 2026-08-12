from __future__ import annotations

import pytest

from iupac_scan import CompiledMotif, Match, compile_motif, scan_sequence


@pytest.mark.parametrize(
    ("symbol", "bases"),
    [
        ("A", "A"),
        ("C", "C"),
        ("G", "G"),
        ("T", "T"),
        ("U", "T"),
        ("R", "AG"),
        ("Y", "CT"),
        ("S", "CG"),
        ("W", "AT"),
        ("K", "GT"),
        ("M", "AC"),
        ("B", "CGT"),
        ("D", "AGT"),
        ("H", "ACT"),
        ("V", "ACG"),
        ("N", "ACGT"),
    ],
)
def test_full_iupac_alphabet_has_set_intersection_semantics(
    symbol: str, bases: str
) -> None:
    observed = {base for base in "ACGT" if list(scan_sequence(base, symbol))}
    assert observed == set(bases)


def test_compile_motif_normalizes_and_exposes_masks() -> None:
    compiled = compile_motif("aRn")
    assert isinstance(compiled, CompiledMotif)
    assert compiled.motif == "ARN"
    assert compiled.position_masks == (1, 5, 15)
    assert compiled.shift_masks == (7, 4, 6, 4)
    assert compiled.accept_bit == 4


def test_overlapping_matches_and_ambiguous_sequence_symbols() -> None:
    assert list(scan_sequence("AANA", "ANA")) == [
        Match(0, 3, "+", "AAN"),
        Match(1, 4, "+", "ANA"),
    ]


def test_both_strands_reverse_complement_and_palindrome_duplicates() -> None:
    assert list(scan_sequence("ATGACAT", "ATG", both_strands=True)) == [
        Match(0, 3, "+", "ATG"),
        Match(4, 7, "-", "CAT"),
    ]
    palindrome = list(scan_sequence("AT", "AT", both_strands=True))
    assert [match.strand for match in palindrome] == ["+", "-"]


def test_reverse_complement_handles_each_single_base_mask() -> None:
    assert compile_motif("ACGT").reverse_complement().motif == "ACGT"


def test_string_motif_is_compiled_and_lowercase_sequence_is_normalized() -> None:
    assert list(scan_sequence("acgt", "cG")) == [Match(1, 3, "+", "CG")]


@pytest.mark.parametrize(
    ("call", "message"),
    [
        (lambda: compile_motif(""), "motif must not be empty"),
        (lambda: compile_motif("AX"), "invalid motif symbol 'X' at offset 1"),
        (lambda: list(scan_sequence("", "A")), "sequence must not be empty"),
        (
            lambda: list(scan_sequence("AZ", "A")),
            "invalid sequence symbol 'Z' at offset 1",
        ),
    ],
)
def test_invalid_input_has_precise_error(call: object, message: str) -> None:
    with pytest.raises(ValueError, match=message.replace("'", "\\'")):
        call()  # type: ignore[operator]
