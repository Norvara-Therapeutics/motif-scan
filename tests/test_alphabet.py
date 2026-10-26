"""Residue alphabets, codes and alphabet detection."""

import pytest

from motif_scan import alphabet as ab


def test_every_nucleotide_code_stands_for_at_least_one_base():
    for code, bases in ab.NUCLEOTIDE_CODES.items():
        assert bases, code
        assert bases <= frozenset("ACGT")


def test_uracil_is_thymine_not_an_ambiguity():
    assert ab.NUCLEOTIDE.expand("U") == frozenset("T")
    assert not ab.has_ambiguity("AUGC", ab.NUCLEOTIDE)


def test_n_stands_for_every_base():
    assert ab.NUCLEOTIDE.expand("N") == frozenset("ACGT")


def test_x_stands_for_the_twenty_standard_amino_acids():
    assert len(ab.PROTEIN.expand("X")) == 20
    assert "U" not in ab.PROTEIN.expand("X")
    assert "*" not in ab.PROTEIN.expand("X")


def test_protein_ambiguity_codes():
    assert ab.PROTEIN.expand("B") == frozenset("DN")
    assert ab.PROTEIN.expand("Z") == frozenset("EQ")
    assert ab.PROTEIN.expand("J") == frozenset("IL")


def test_stop_stands_only_for_itself():
    assert ab.PROTEIN.expand("*") == frozenset("*")


def test_definite_codes_stand_for_one_residue():
    for code in ab.NUCLEOTIDE.definite:
        assert len(ab.NUCLEOTIDE.expand(code)) == 1
    assert "N" not in ab.NUCLEOTIDE.definite
    assert "X" not in ab.PROTEIN.definite


def test_wildcard_of_each_alphabet():
    assert ab.NUCLEOTIDE.anything == frozenset("ACGT")
    assert len(ab.PROTEIN.anything) == 20


def test_expand_rejects_a_character_that_is_not_a_code():
    with pytest.raises(KeyError, match="not a nucleotide code"):
        ab.NUCLEOTIDE.expand("Q")


def test_expand_is_case_insensitive():
    assert ab.NUCLEOTIDE.expand("r") == ab.NUCLEOTIDE.expand("R")


def test_has_ambiguity_finds_lowercase_codes():
    assert ab.has_ambiguity("acgtn", ab.NUCLEOTIDE)
    assert not ab.has_ambiguity("acgt", ab.NUCLEOTIDE)


def test_has_ambiguity_in_protein():
    assert ab.has_ambiguity("MKXV", ab.PROTEIN)
    assert not ab.has_ambiguity("MKAV", ab.PROTEIN)


def test_complement_preserves_case_and_codes():
    assert ab.complement("acgtRY") == "tgcaYR"


def test_reverse_complement_round_trips():
    sequence = "ACGTRYSWKMBDHVN"
    assert ab.reverse_complement(ab.reverse_complement(sequence)) == sequence


def test_palindrome_recognises_a_degenerate_one():
    assert ab.is_palindrome("GAATTC")
    assert ab.is_palindrome("GGWCC")
    assert not ab.is_palindrome("GAGCGAGC")


def test_empty_pattern_is_not_a_palindrome():
    assert not ab.is_palindrome("")


def test_matching_is_containment_not_overlap():
    # An ambiguous motif position accepts a definite base.
    assert ab.matches("N", "A", ab.NUCLEOTIDE)
    # An ambiguous sequence residue does not satisfy a definite position.
    assert not ab.matches("A", "N", ab.NUCLEOTIDE)
    # It satisfies a position that covers everything it might be.
    assert ab.matches("N", "N", ab.NUCLEOTIDE)
    assert ab.matches("R", "A", ab.NUCLEOTIDE)
    assert not ab.matches("R", "Y", ab.NUCLEOTIDE)


def test_could_match_is_overlap():
    assert ab.could_match("A", "N", ab.NUCLEOTIDE)
    assert ab.could_match("R", "N", ab.NUCLEOTIDE)
    assert not ab.could_match("A", "C", ab.NUCLEOTIDE)
    assert not ab.could_match("R", "Y", ab.NUCLEOTIDE)


def test_detect_reads_dna_as_nucleotide():
    assert ab.detect("ACGTACGTACGTACGTNNNN") == "nucleotide"


def test_detect_reads_a_peptide_as_protein():
    assert ab.detect("MKQLEDKVEELLSKNYHLENEVARLKKLV") == "protein"


def test_detect_settles_on_an_amino_acid_with_no_nucleotide_reading():
    assert ab.detect("ACGTACGTE") == "protein"


def test_detect_has_nothing_to_judge_in_an_empty_sequence():
    assert ab.detect("") is None


def test_detect_is_case_insensitive():
    assert ab.detect("acgtacgtacgt") == "nucleotide"


def test_detection_share_is_the_published_rule():
    assert ab.NUCLEOTIDE_SHARE == 0.85


def test_protein_only_residues_have_no_nucleotide_reading():
    assert ab.PROTEIN_ONLY & frozenset(ab.NUCLEOTIDE_CODES) == frozenset()
    assert "E" in ab.PROTEIN_ONLY
    assert "X" in ab.PROTEIN_ONLY


def test_every_nucleotide_code_is_also_an_amino_acid_letter():
    # The reason detection has to be by composition rather than by letter.
    assert frozenset(ab.NUCLEOTIDE_CODES) <= frozenset(ab.PROTEIN_CODES)


def test_get_by_name():
    assert ab.get("protein") is ab.PROTEIN
    with pytest.raises(KeyError, match="not an alphabet"):
        ab.get("rna")


def test_nucleotide_has_strands_and_protein_does_not():
    assert ab.NUCLEOTIDE.has_strands
    assert not ab.PROTEIN.has_strands
