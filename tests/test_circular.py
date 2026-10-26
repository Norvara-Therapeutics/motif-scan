"""Circular topology: a motif may span the join."""

import pytest

from motif_scan.alphabet import reverse_complement
from motif_scan.patterns import Motif, PatternError, compile_motif
from motif_scan.search import count_uncertain, find


def motif(pattern, syntax="iupac", name="m", alphabet="nucleotide"):
    return compile_motif(
        Motif(name=name, pattern=pattern, syntax=syntax, alphabet=alphabet)
    )


def test_a_motif_spanning_the_join_is_found():
    # TTCGAA split across the join: the sequence ends GAA and begins TTC.
    sequence = "GAA" + "CCCCCC" + "TTC"
    (hit,) = find(sequence, motif("TTCGAA"), circular=True)
    assert (hit.start, hit.end) == (10, 3)
    assert hit.spans_origin
    assert hit.matched == "TTCGAA"


def test_the_same_motif_is_not_found_on_a_linear_sequence():
    sequence = "GAA" + "CCCCCC" + "TTC"
    assert find(sequence, motif("TTCGAA"), circular=False) == []


def test_a_hit_is_not_reported_twice_from_the_appended_start():
    sequence = "ACGTACGTACGT"
    linear = find(sequence, motif("ACGT"), circular=False)
    circular = find(sequence, motif("ACGT"), circular=True)
    # The circle adds the wrapped occurrences and keeps the rest once each.
    assert len(circular) >= len(linear)
    assert len(circular) == len({(hit.start, hit.strand) for hit in circular})


def test_the_reverse_strand_also_wraps():
    motif_text = "GCTCGCTC"
    sequence = (
        reverse_complement(motif_text)[4:] + "T" * 12 + reverse_complement(motif_text)[:4]
    )
    (hit,) = [
        hit
        for hit in find(sequence, motif(motif_text), circular=True)
        if hit.spans_origin
    ]
    assert hit.strand == "-"
    assert hit.matched == motif_text


def test_window_count_on_a_circle_is_the_whole_length():
    # Every position is a possible start on a circle.
    sequence = "ACGT" * 4
    hits = find(sequence, motif("ACGTACGT"), circular=True, strand="forward")
    assert len(hits) == len(sequence) // 4


def test_a_motif_longer_than_the_sequence_still_has_no_hits():
    assert find("ACGT", motif("ACGTACGTACGT"), circular=True) == []


def test_uncertain_windows_are_counted_across_the_join():
    sequence = "AAN" + "C" * 6 + "TTA"
    assert count_uncertain(sequence, motif("TTAAAA"), circular=True) >= 1


def test_an_anchored_pattern_cannot_be_searched_on_a_circle():
    with pytest.raises(PatternError, match="no termini"):
        find(
            "MAKDEL",
            motif("[KR]-D-E-L>", "prosite", alphabet="protein"),
            circular=True,
        )


def test_an_unbounded_pattern_cannot_be_searched_on_a_circle():
    with pytest.raises(PatternError, match="unbounded length"):
        find("ACGTACGT", motif("A.+T", "regex"), circular=True)


def test_a_bounded_prosite_pattern_can_be():
    hits = find("TTTACGTTTT", motif("A-x(1,2)-T", "prosite"), circular=True)
    assert hits


def test_the_join_is_only_crossed_once_for_a_short_sequence():
    # The appended window is capped at the length of the sequence, so a motif
    # longer than the molecule cannot wrap around it repeatedly.
    assert find("AT", motif("ATATAT"), circular=True) == []


def test_overlap_suppression_takes_in_the_join():
    # AAA tiles around a circle of four As: every start overlaps the first hit,
    # including the ones that reach it the long way round.
    # (TTT does not occur, so the reverse strand contributes nothing here.)
    assert [
        (hit.start, hit.end)
        for hit in find("AAAA", motif("AAA"), circular=True, overlapping=False)
    ] == [(1, 3)]
    # Without suppression every start is a hit, including those that wrap.
    assert len(find("AAAA", motif("AAA"), circular=True)) == 4


def test_suppressed_hits_on_a_circle_share_no_residue():
    from motif_scan.search import covered_by

    sequence = "ACGTACGTACGTACG"
    hits = find(sequence, motif("GTAC"), circular=True, overlapping=False)
    for strand in {hit.strand for hit in hits}:
        occupied: set[int] = set()
        for hit in (hit for hit in hits if hit.strand == strand):
            residues = covered_by(hit, len(sequence))
            assert not residues & occupied
            occupied |= residues


def test_covered_residues_of_a_wrapping_hit():
    from motif_scan.search import Hit, covered_by

    assert covered_by(Hit("r", "m", 4, 2, "+", "AAA"), 4) == {3, 0, 1}
    assert covered_by(Hit("r", "m", 1, 3, "+", "AAA"), 4) == {0, 1, 2}


def test_overlap_suppression_still_keeps_both_strands():
    # A hit on each strand covering the same residues is two findings.
    sequence = "GCTCGCTC"
    hits = find(sequence, motif("GCTCGCTC"), overlapping=False)
    forward = find(sequence, motif("GCTCGCTC"), strand="forward", overlapping=False)
    assert len(hits) >= len(forward)


def test_a_bounded_pattern_cannot_lap_the_molecule():
    # A(2,8) on a circle of four adenines: greedy matching would run round the
    # circle and past its own start, reporting eight residues in a span of four.
    hits = find("AAAA", motif("A(2,8)", "prosite"), circular=True, strand="forward")
    assert hits
    for hit in hits:
        assert hit.length <= 4
        assert len(hit.matched) == hit.length
        assert hit.matched == "AAAA"


def test_coordinates_and_matched_text_always_agree_on_a_circle():
    sequence = "ACGTACGT"
    for pattern in ("A-x(1,6)-T", "A(1,9)", "[ACGT](2,12)"):
        compiled = motif(pattern, "prosite")
        for hit in find(sequence, compiled, circular=True):
            covered = (
                hit.end - hit.start + 1
                if not hit.spans_origin
                else len(sequence) - hit.start + 1 + hit.end
            )
            assert covered == len(hit.matched), (pattern, hit)


def test_a_bounded_pattern_shorter_than_the_molecule_is_unaffected():
    hits = find("ACGTACGTACGT", motif("A-x(1,2)-T", "prosite"), circular=True)
    assert hits
    assert all(3 <= hit.length <= 4 for hit in hits)
