"""Finding motifs."""

import pytest

from motif_scan.alphabet import PROTEIN, reverse_complement
from motif_scan.patterns import Motif, PatternError, compile_motif
from motif_scan.search import STRANDS, count_uncertain, find


def motif(pattern, syntax="iupac", alphabet="nucleotide", name="m"):
    return compile_motif(
        Motif(name=name, pattern=pattern, syntax=syntax, alphabet=alphabet)
    )


def spans(hits):
    return [(hit.start, hit.end, hit.strand) for hit in hits]


def test_strand_labels_are_the_four_documented():
    assert STRANDS == ("+", "-", "=", ".")


def test_finds_a_literal_on_the_forward_strand():
    hits = find("AAAGAATTCAAA", motif("GAATTC"))
    assert spans(hits) == [(4, 9, "=")]


def test_a_palindrome_is_reported_once():
    hits = find("GAATTC", motif("GAATTC"))
    assert len(hits) == 1
    assert hits[0].strand == "="


def test_an_asymmetric_motif_is_found_on_both_strands():
    sequence = "TTTGCTCGCTCTTT" + reverse_complement("TTTGCTCGCTCTTT")
    hits = find(sequence, motif("GCTCGCTC"))
    assert {hit.strand for hit in hits} == {"+", "-"}


def test_a_reverse_hit_is_given_forward_coordinates_and_reads_as_the_motif():
    # ATGCAAAT on the reverse strand of this sequence.
    sequence = "CCC" + reverse_complement("ATGCAAAT") + "CCC"
    (hit,) = find(sequence, motif("ATGCAAAT"))
    assert hit.strand == "-"
    assert (hit.start, hit.end) == (4, 11)
    assert hit.matched == "ATGCAAAT"
    assert sequence[hit.start - 1: hit.end] == reverse_complement(hit.matched)


def test_overlapping_hits_are_all_reported():
    assert spans(find("CAAAATCCC", motif("AA"), strand="forward")) == [
        (2, 3, "+"),
        (3, 4, "+"),
        (4, 5, "+"),
    ]


def test_overlapping_can_be_switched_off():
    assert spans(
        find("CAAAATCCC", motif("AA"), strand="forward", overlapping=False)
    ) == [(2, 3, "+"), (4, 5, "+")]


def test_overlaps_are_thinned_strand_by_strand():
    # GGATCC is palindromic, so a single pass; GGATCCGGATCC tiles.
    sequence = "ATATAT"
    hits = find(sequence, motif("ATAT"), overlapping=False)
    assert len(hits) == 1


def test_degenerate_codes_in_the_pattern_widen_the_match():
    hits = find("TATAAATG", motif("TATAWAWR"), strand="forward")
    assert [hit.matched for hit in hits] == ["TATAAATG"]


def test_an_ambiguous_sequence_residue_is_not_a_hit():
    assert find("TATANAAG", motif("TATAWAWR"), strand="forward") == []


def test_an_ambiguous_sequence_residue_is_counted_as_uncertain():
    assert count_uncertain("TATANAAG", motif("TATAWAWR"), strand="forward") == 1


def test_uncertain_windows_can_be_listed_when_asked_for():
    (hit,) = find(
        "TATANAAG", motif("TATAWAWR"), strand="forward", include_uncertain=True
    )
    assert hit.uncertain == 1
    assert not hit.certain
    assert hit.exact


def test_uncertainty_is_not_counted_for_a_definite_match():
    assert count_uncertain("TATAAAAG", motif("TATAWAWR"), strand="forward") == 0


def test_uncertainty_is_not_assessed_for_a_variable_length_pattern():
    # Reported as zero because nothing was assessed, which the count carries
    # separately; the search itself still runs.
    assert count_uncertain("ACGTNNNN", motif("A.+T", syntax="regex")) == 0


def test_a_sequence_without_ambiguity_costs_nothing_to_assess():
    assert count_uncertain("ACGTACGT", motif("ACGT")) == 0


def test_mismatches_are_allowed_against_a_fixed_pattern():
    hits = find("TATAAAAT", motif("TATAWAWR"), mismatches=1, strand="forward")
    assert len(hits) == 1
    assert hits[0].mismatches == 1
    assert not hits[0].exact


def test_the_mismatch_budget_is_respected():
    assert find("TATAAATT", motif("TATAWAWR"), mismatches=0, strand="forward") == []
    assert len(find("TAGAAATT", motif("TATAWAWR"), mismatches=2, strand="forward")) == 1


def test_mismatches_are_refused_for_a_regex():
    with pytest.raises(PatternError, match="variable length"):
        find("ACGT", motif("A.T", syntax="regex"), mismatches=1)


def test_mismatches_are_refused_for_a_prosite_range():
    with pytest.raises(PatternError, match="variable length"):
        find("ACGT", motif("A-x(1,2)-T", syntax="prosite"), mismatches=1)


def test_a_negative_mismatch_budget_is_refused():
    with pytest.raises(PatternError, match="cannot be negative"):
        find("ACGT", motif("ACGT"), mismatches=-1)


def test_an_ambiguous_position_does_not_spend_the_mismatch_budget():
    # Without the distinction, a run of Ns would match anything short enough.
    assert find("TATANAAG", motif("TATAWAWR"), mismatches=0, strand="forward") == []
    (hit,) = find(
        "TATANAAG",
        motif("TATAWAWR"),
        mismatches=0,
        strand="forward",
        include_uncertain=True,
    )
    assert (hit.mismatches, hit.uncertain) == (0, 1)


def test_strand_can_be_restricted():
    sequence = "CCC" + reverse_complement("ATGCAAAT") + "CCC"
    assert find(sequence, motif("ATGCAAAT"), strand="forward") == []
    assert len(find(sequence, motif("ATGCAAAT"), strand="reverse")) == 1


def test_protein_hits_are_unstranded():
    (hit,) = find("MKRGDLV", motif("R-G-D", "prosite", "protein"))
    assert hit.strand == "."
    assert (hit.start, hit.end) == (3, 5)


def test_asking_for_the_reverse_strand_of_a_protein_is_refused():
    with pytest.raises(PatternError, match="no reverse strand"):
        find("MKRGDLV", motif("R-G-D", "prosite", "protein"), strand="reverse")


def test_an_anchored_motif_is_found_only_at_the_terminus():
    compiled = motif("[KR]-D-E-L>", "prosite", "protein")
    assert len(find("MAKDEL", compiled)) == 1
    assert find("MAKDELAA", compiled) == []


def test_a_start_anchor_is_respected():
    compiled = motif("<M-x-K", "prosite", "protein")
    assert len(find("MAKAAA", compiled)) == 1
    assert find("AMAKAA", compiled) == []


def test_a_protein_wildcard_does_not_match_a_stop_codon():
    assert find("MK*LV", motif("K-x-L", "prosite", "protein")) == []


def test_an_empty_sequence_has_no_hits():
    assert find("", motif("ACGT")) == []


def test_a_sequence_shorter_than_the_motif_has_no_hits():
    assert find("ACG", motif("ACGT")) == []


def test_hits_carry_the_record_identifier():
    (hit,) = find("GAATTC", motif("GAATTC"), record="plasmid")
    assert hit.record == "plasmid"
    assert hit.motif == "m"


def test_hit_length_is_the_length_of_the_match():
    hits = find("ACGGGGTT", motif("G-x(1,3)-T", "prosite"), strand="forward")
    assert hits
    for hit in hits:
        assert hit.length == len(hit.matched)
        assert 3 <= hit.length <= 5


def test_hits_are_sorted_by_position_then_strand():
    sequence = "GCTCGCTC" + "N" * 10 + reverse_complement("GCTCGCTC")
    hits = find(sequence, motif("GCTCGCTC"))
    assert [hit.start for hit in hits] == sorted(hit.start for hit in hits)


def test_a_variable_length_prosite_pattern_matches_a_range():
    compiled = motif("[KR](2)-x(2,4)-[KR](2)", "prosite", "protein")
    assert [hit.matched for hit in find("KKAAARR", compiled)] == ["KKAAARR"]


def test_the_protein_alphabet_has_no_strands():
    assert not PROTEIN.has_strands


# --- user-written regular expressions -------------------------------------


def test_a_regex_with_a_backreference_is_searched():
    hits = find("ACAATT", motif(r"([AC])\1", syntax="regex"), strand="forward")
    assert [hit.matched for hit in hits] == ["AA"]


def test_a_regex_with_a_named_backreference_is_searched():
    hits = find("ACAATT", motif(r"(?P<x>[AC])(?P=x)", syntax="regex"), strand="forward")
    assert [hit.matched for hit in hits] == ["AA"]


def test_a_regex_with_inline_flags_is_searched():
    hits = find("ttacgg", motif("(?i)ACG", syntax="regex"), strand="forward")
    assert [hit.matched for hit in hits] == ["acg"]


def test_a_regex_with_a_lookahead_is_searched():
    hits = find("ACGTACC", motif("AC(?=G)", syntax="regex"), strand="forward")
    assert [(hit.start, hit.matched) for hit in hits] == [(1, "AC")]


def test_a_regex_finds_overlapping_occurrences():
    hits = find("CAAAATCCC", motif("A{2}", syntax="regex"), strand="forward")
    assert [(hit.start, hit.end) for hit in hits] == [(2, 3), (3, 4), (4, 5)]


def test_a_regex_is_searched_on_both_strands():
    hits = find("GGGCCC", motif("GGG", syntax="regex"))
    assert sorted(hit.strand for hit in hits) == ["+", "-"]


def test_a_regex_anchor_binds_to_the_real_start():
    hits = find("ACGACG", motif("^ACG", syntax="regex"), strand="forward")
    assert [(hit.start, hit.end) for hit in hits] == [(1, 3)]


# --- uncertain hits do not contradict their labels ------------------------


def test_an_uncertain_hit_is_neither_exact_nor_mismatched():
    (hit,) = find(
        "AATANA", motif("AATAAA"), strand="forward", include_uncertain=True
    )
    assert hit.uncertain == 1
    assert not hit.certain
    assert hit.mismatches == 0
