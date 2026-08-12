from collections import Counter

from hypothesis import assume, given
from hypothesis import strategies as st

from iupac_scan import compile_motif, scan_sequence

_IUPAC_BASES = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "U": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}
_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
_SEQUENCES = st.text(alphabet="ACGTNacgt", max_size=60)
_MOTIFS = st.text(
    alphabet="ACGTURYSWKMBDHVN",
    min_size=1,
    max_size=8,
)


@st.composite
def _sequence_and_valid_width(
    draw: st.DrawFn,
) -> tuple[str, int]:
    sequence = draw(st.text(alphabet="ACGT", min_size=1, max_size=60))
    width = draw(st.integers(min_value=1, max_value=min(8, len(sequence))))
    return sequence, width


def _naive_matches(
    sequence: str, motif: str, *, both_strands: bool
) -> set[tuple[int, int, str]]:
    sequence_sets = [_IUPAC_BASES[symbol] for symbol in sequence.upper()]
    motif_sets = [_IUPAC_BASES[symbol] for symbol in motif]
    oriented_motifs = [(motif_sets, "+")]
    if both_strands:
        reverse_complement = [
            frozenset(_COMPLEMENT[base] for base in bases)
            for bases in reversed(motif_sets)
        ]
        oriented_motifs.append((reverse_complement, "-"))

    matches = set()
    width = len(motif_sets)
    for oriented, strand in oriented_motifs:
        for start in range(len(sequence_sets) - width + 1):
            if all(
                sequence_sets[start + offset] & motif_bases
                for offset, motif_bases in enumerate(oriented)
            ):
                matches.add((start, start + width, strand))
    return matches


def _reverse_complement(sequence: str) -> str:
    return "".join(_COMPLEMENT[symbol] for symbol in reversed(sequence.upper()))


@given(sequence=_SEQUENCES, motif=_MOTIFS, both_strands=st.booleans())
def test_scan_matches_naive_iupac_oracle(
    sequence: str, motif: str, both_strands: bool
) -> None:
    assume(sequence)
    observed = {
        (match.start, match.end, match.strand)
        for match in scan_sequence(
            sequence, compile_motif(motif), both_strands=both_strands
        )
    }
    assert observed == _naive_matches(sequence, motif, both_strands=both_strands)


@given(sequence=_SEQUENCES, motif=_MOTIFS)
def test_both_strands_is_symmetric_under_sequence_reverse_complement(
    sequence: str, motif: str
) -> None:
    assume(sequence)
    length = len(sequence)
    original = {
        (match.start, match.end, match.strand)
        for match in scan_sequence(sequence, compile_motif(motif), both_strands=True)
    }
    reversed_observed = {
        (match.start, match.end, match.strand)
        for match in scan_sequence(
            _reverse_complement(sequence),
            compile_motif(motif),
            both_strands=True,
        )
    }
    mirrored = {
        (length - end, length - start, "-" if strand == "+" else "+")
        for start, end, strand in original
    }
    assert reversed_observed == mirrored


@given(case=_sequence_and_valid_width())
def test_wildcard_motif_matches_every_window_per_strand(
    case: tuple[str, int],
) -> None:
    sequence, width = case
    matches = scan_sequence(sequence, compile_motif("N" * width), both_strands=True)
    counts = Counter(match.strand for match in matches)
    expected = max(0, len(sequence) - width + 1)
    assert counts["+"] == expected
    assert counts["-"] == expected
    assert counts.total() == 2 * expected
