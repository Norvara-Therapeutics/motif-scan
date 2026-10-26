"""Counting, expecting and distributing hits."""

from collections import Counter

import pytest

from motif_scan.patterns import Motif, compile_motif
from motif_scan.search import Hit, find
from motif_scan.summary import (
    DEFAULT_BINS,
    MotifCount,
    distribute,
    expected_count,
    frequencies,
    tally,
)


def motif(pattern, syntax="iupac", alphabet="nucleotide", name="m"):
    return compile_motif(
        Motif(name=name, pattern=pattern, syntax=syntax, alphabet=alphabet)
    )


EVEN = Counter({"A": 250, "C": 250, "G": 250, "T": 250})


def count(hits=0, exact=0, expected=None, residues=1000, **extra):
    """A MotifCount with only the fields a test cares about."""
    fields = dict(
        motif="m",
        hits=hits,
        forward=hits,
        reverse=0,
        palindromic=0,
        unstranded=0,
        exact=exact,
        with_mismatches=0,
        listed_uncertain=0,
        uncertain=0,
        assessed=True,
        expected=expected,
        residues=residues,
    )
    fields.update(extra)
    return MotifCount(**fields)


def test_frequencies_are_over_the_residues_that_are_known():
    frequency = frequencies(Counter({"A": 30, "C": 10, "N": 60}), motif("A"))
    assert frequency == {"A": 0.75, "C": 0.25}


def test_frequencies_count_uracil_as_thymine():
    frequency = frequencies(Counter({"U": 10, "A": 10}), motif("A"))
    assert frequency == {"T": 0.5, "A": 0.5}


def test_frequencies_are_none_without_a_definite_residue():
    assert frequencies(Counter({"N": 10}), motif("A")) is None


def test_expected_count_of_a_literal_under_even_composition():
    # One of 4^6 windows, over 1000 - 6 + 1 windows, on each of two strands.
    expected = expected_count(motif("AATAAA"), EVEN, 1000)
    assert expected == pytest.approx(2 * 995 / 4096, rel=1e-9)


def test_expected_count_halves_for_one_strand():
    both = expected_count(motif("AATAAA"), EVEN, 1000)
    one = expected_count(motif("AATAAA"), EVEN, 1000, strand="forward")
    assert one == pytest.approx(both / 2)


def test_a_palindrome_is_expected_once_not_twice():
    expected = expected_count(motif("GAATTC"), EVEN, 1000)
    assert expected == pytest.approx(995 / 4096, rel=1e-9)


def test_degenerate_positions_raise_the_expectation():
    # CANNTG: four fixed positions at 1/4 and two free ones.
    expected = expected_count(motif("CANNTG"), EVEN, 1000, strand="forward")
    assert expected == pytest.approx(995 / 256, rel=1e-9)


def test_composition_moves_the_expectation():
    gc_rich = Counter({"A": 100, "T": 100, "C": 400, "G": 400})
    at_rich = Counter({"A": 400, "T": 400, "C": 100, "G": 100})
    gc_motif = motif("GGGCCC")
    assert expected_count(gc_motif, gc_rich, 1000) > expected_count(
        gc_motif, at_rich, 1000
    )


def test_mismatches_raise_the_expectation_exactly():
    # With one mismatch allowed over six independent positions at p = 1/4:
    # P = p^6 + 6 p^5 (1 - p).
    one_strand = expected_count(
        motif("AATAAA"), EVEN, 1000, mismatches=1, strand="forward"
    )
    probability = (0.25 ** 6) + 6 * (0.25 ** 5) * 0.75
    assert one_strand == pytest.approx(995 * probability, rel=1e-9)


def test_expectation_matches_an_exhaustive_count_on_a_small_alphabet():
    # Enumerate every 8-mer and count the ones the pattern accepts, which must
    # agree with the dynamic programme the expectation uses.
    from itertools import product

    compiled = motif("TATAWAWR")
    accepted = sum(
        1
        for word in product("ACGT", repeat=8)
        if all(residue in allowed for residue, allowed in zip(word, compiled.positions))
    )
    expected = expected_count(compiled, EVEN, 8, strand="forward")
    assert expected == pytest.approx(accepted / 4 ** 8, rel=1e-9)


def test_expectation_with_mismatches_matches_an_exhaustive_count():
    from itertools import product

    compiled = motif("GAATTC")
    accepted = sum(
        1
        for word in product("ACGT", repeat=6)
        if sum(
            residue not in allowed
            for residue, allowed in zip(word, compiled.positions)
        )
        <= 1
    )
    expected = expected_count(compiled, EVEN, 6, mismatches=1, strand="forward")
    assert expected == pytest.approx(accepted / 4 ** 6, rel=1e-9)


def test_a_circular_sequence_has_one_window_per_residue():
    linear = expected_count(motif("AATAAA"), EVEN, 1000, strand="forward")
    circular = expected_count(
        motif("AATAAA"), EVEN, 1000, strand="forward", circular=True
    )
    assert circular == pytest.approx(linear * 1000 / 995)


def test_an_anchored_pattern_has_one_window():
    compiled = motif("[KR]-D-E-L>", "prosite", "protein")
    composition = Counter({residue: 10 for residue in "ACDEFGHIKLMNPQRSTVWY"})
    assert expected_count(compiled, composition, 500) == pytest.approx(
        (2 / 20) * (1 / 20) ** 3
    )


def test_no_expectation_for_a_variable_length_pattern():
    assert expected_count(motif("A.+T", "regex"), EVEN, 1000) is None


def test_no_expectation_without_a_known_residue():
    assert expected_count(motif("AATAAA"), Counter({"N": 1000}), 1000) is None


def test_a_sequence_shorter_than_the_motif_expects_nothing():
    assert expected_count(motif("AATAAA"), EVEN, 3) == 0.0


def test_the_count_fields_do_not_contradict_their_labels():
    # An uncertain window with no mismatch is neither exact nor mismatched.
    compiled = motif("AATAAA")
    hits = find("AATANA", compiled, strand="forward", include_uncertain=True)
    tallied = tally(compiled, hits, 6, Counter("AATANA"))
    assert tallied.hits == 1
    assert tallied.exact == 0
    assert tallied.with_mismatches == 0
    assert tallied.listed_uncertain == 1


def test_tally_counts_strands_and_exactness():
    compiled = motif("GCTCGCTC")
    sequence = "GCTCGCTC" + "A" * 20 + "GAGCGAGC"
    hits = find(sequence, compiled)
    count = tally(compiled, hits, len(sequence), Counter(sequence))
    assert count.hits == 2
    assert count.forward == 1
    assert count.reverse == 1
    assert count.exact == 2
    assert count.with_mismatches == 0
    assert count.assessed


def test_tally_reports_density_per_thousand_residues():
    assert count(hits=5, exact=5).per_kb == pytest.approx(5.0)


def test_density_of_an_empty_record_is_zero():
    assert count(residues=0).per_kb == 0.0


def test_enrichment_is_observed_over_expected():
    assert count(hits=8, exact=8, expected=2.0).enrichment == 4.0


def test_enrichment_is_undefined_without_an_expectation():
    assert count(hits=8, exact=8).enrichment is None
    assert count(hits=8, exact=8, expected=0.0).enrichment is None


def test_uncertainty_is_not_assessed_for_a_variable_length_pattern():
    compiled = motif("A.+T", "regex")
    count = tally(compiled, [], 100, EVEN)
    assert not count.assessed
    assert count.uncertain == 0


def test_distribution_places_hits_by_their_start():
    hits = [
        Hit("r", "m", 1, 2, "+", "AC"),
        Hit("r", "m", 99, 100, "+", "AC"),
    ]
    distribution = distribute(hits, ["m"], 100, bins=2)
    (name, counts) = distribution.counts[0]
    assert name == "m"
    assert counts == (1, 1)


def test_distribution_keeps_a_row_for_a_motif_without_hits():
    distribution = distribute([], ["a", "b"], 100, bins=4)
    assert [name for name, _ in distribution.counts] == ["a", "b"]
    assert distribution.busiest == 0


def test_distribution_is_never_finer_than_one_residue_per_bin():
    distribution = distribute([], ["m"], 5, bins=DEFAULT_BINS)
    assert distribution.bins == 5


def test_distribution_of_an_empty_record():
    distribution = distribute([], ["m"], 0, bins=10)
    assert distribution.bins == 0
    assert distribution.width == 0.0


def test_distribution_spans_cover_the_record():
    distribution = distribute([], ["m"], 100, bins=4)
    assert distribution.span(0) == (1, 25)
    assert distribution.span(3)[1] == 100


def test_distribution_ignores_a_motif_it_was_not_given():
    distribution = distribute([Hit("r", "other", 1, 2, "+", "AC")], ["m"], 100, 4)
    assert distribution.busiest == 0


# --- the expectation is over the windows that could be called -------------


def test_assessable_windows_skips_those_holding_an_ambiguity_code():
    from motif_scan.alphabet import NUCLEOTIDE
    from motif_scan.summary import assessable_windows

    assert assessable_windows("ANNNNNNNNN", 1, NUCLEOTIDE) == 1
    assert assessable_windows("ACGTACGT", 4, NUCLEOTIDE) == 5
    # Of ACNT CNTA NTAC TACG ACGT, the last two carry no N.
    assert assessable_windows("ACNTACGT", 4, NUCLEOTIDE) == 2
    assert assessable_windows("NNNN", 2, NUCLEOTIDE) == 0


def test_assessable_windows_wrap_on_a_circle():
    from motif_scan.alphabet import NUCLEOTIDE
    from motif_scan.summary import assessable_windows

    assert assessable_windows("ACGT", 2, NUCLEOTIDE, circular=True) == 4
    assert assessable_windows("ACGT", 2, NUCLEOTIDE) == 3
    assert assessable_windows("ACGN", 2, NUCLEOTIDE, circular=True) == 2


def test_assessable_windows_of_a_protein():
    from motif_scan.alphabet import PROTEIN
    from motif_scan.summary import assessable_windows

    assert assessable_windows("MKXLV", 2, PROTEIN) == 2


def test_the_expectation_can_be_taken_over_given_windows():
    compiled = motif("A")
    every = expected_count(compiled, EVEN, 10, strand="forward")
    assessable = expected_count(compiled, EVEN, 10, strand="forward", windows=1)
    assert every == pytest.approx(10 * 0.25)
    assert assessable == pytest.approx(0.25)


def test_ambiguity_does_not_inflate_the_expectation():
    # One adenine then nine Ns: one hit, and one window that could be called.
    from motif_scan.alphabet import NUCLEOTIDE
    from motif_scan.summary import assessable_windows

    sequence = "ANNNNNNNNN"
    compiled = motif("A")
    hits = find(sequence, compiled, strand="forward")
    windows = assessable_windows(sequence, 1, NUCLEOTIDE)
    tallied = tally(
        compiled,
        hits,
        len(sequence),
        Counter(sequence),
        strand="forward",
        windows=windows,
    )
    assert tallied.hits == 1
    assert tallied.expected == pytest.approx(1.0)
    assert tallied.enrichment == pytest.approx(1.0)
