"""Reading and compiling patterns."""

import pytest

from motif_scan.alphabet import NUCLEOTIDE, PROTEIN
from motif_scan.patterns import (
    GROUP,
    MAX_REPEAT,
    SYNTAXES,
    Motif,
    PatternError,
    compile_motif,
    complement_positions,
    named,
    parse_iupac,
    parse_prosite,
)


def nucleotide(pattern, syntax="iupac", name="m"):
    return compile_motif(Motif(name=name, pattern=pattern, syntax=syntax))


def protein(pattern, syntax="prosite", name="m"):
    return compile_motif(
        Motif(name=name, pattern=pattern, syntax=syntax, alphabet="protein")
    )


def test_syntaxes_are_the_three_documented():
    assert SYNTAXES == ("iupac", "prosite", "regex")


def test_a_motif_needs_a_name_and_a_pattern():
    with pytest.raises(PatternError, match="must have a name"):
        Motif(name="", pattern="ACGT")
    with pytest.raises(PatternError, match="pattern is empty"):
        Motif(name="m", pattern="")


def test_an_unknown_syntax_is_refused():
    with pytest.raises(PatternError, match="not a pattern syntax"):
        Motif(name="m", pattern="ACGT", syntax="glob")


def test_iupac_gives_one_set_per_position():
    positions = parse_iupac("TATAWAWR", NUCLEOTIDE)
    assert len(positions) == 8
    assert positions[4] == frozenset("AT")
    assert positions[7] == frozenset("AG")


def test_iupac_is_case_insensitive():
    assert parse_iupac("acgt", NUCLEOTIDE) == parse_iupac("ACGT", NUCLEOTIDE)


def test_iupac_refuses_a_character_that_is_not_a_code():
    with pytest.raises(PatternError, match="not nucleotide codes"):
        parse_iupac("ACGTQ", NUCLEOTIDE)


def test_iupac_points_at_prosite_for_bracket_notation():
    with pytest.raises(PatternError, match="needs --prosite"):
        parse_iupac("[ST]AC", NUCLEOTIDE)


def test_iupac_fixed_length_and_mismatch_support():
    compiled = nucleotide("CANNTG")
    assert compiled.fixed_length == 6
    assert compiled.supports_mismatches


def test_palindrome_is_decided_over_the_degenerate_sets():
    assert nucleotide("CANNTG").palindromic
    assert nucleotide("GGWCC").palindromic
    assert not nucleotide("GCTCGCTC").palindromic


def test_a_protein_motif_is_never_palindromic():
    assert not protein("R-G-D").palindromic


def test_complement_positions_mirrors_the_pattern():
    positions = parse_iupac("GAATTC", NUCLEOTIDE)
    assert complement_positions(positions, NUCLEOTIDE) == positions
    assert complement_positions(parse_iupac("AC", NUCLEOTIDE), NUCLEOTIDE) == (
        frozenset("G"),
        frozenset("T"),
    )


def test_complement_positions_is_empty_for_one_stranded_alphabet():
    assert complement_positions(parse_iupac("RGD", PROTEIN), PROTEIN) == ()


def test_prosite_residue_sets_and_exclusions():
    positions, _, _, _, longest = parse_prosite("N-{P}-[ST]-{P}", PROTEIN)
    assert len(positions) == 4
    assert positions[0] == frozenset("N")
    assert "P" not in positions[1]
    assert positions[2] == frozenset("ST")
    assert longest == 4


def test_prosite_wildcard_and_exact_repeat():
    positions, _, _, _, _ = parse_prosite("[ST]-x(2)-[DE]", PROTEIN)
    assert len(positions) == 4
    assert positions[1] == PROTEIN.anything
    assert positions[3] == frozenset("DE")


def test_prosite_repeat_range_has_no_fixed_length():
    positions, body, _, _, longest = parse_prosite("[KR](2)-x(10,12)-[KR](3,5)", PROTEIN)
    assert positions == ()
    assert "{10,12}" in body
    assert longest == 19


def test_prosite_anchors():
    _, _, start, end, _ = parse_prosite("<M-x(2)-K>", PROTEIN)
    assert start and end
    compiled = protein("[KR]-D-E-L>")
    assert compiled.anchor_end and not compiled.anchor_start
    assert compiled.anchored


def test_prosite_trailing_period_is_accepted():
    assert parse_prosite("R-G-D.", PROTEIN)[0] == parse_prosite("R-G-D", PROTEIN)[0]


def test_prosite_refuses_an_empty_token():
    with pytest.raises(PatternError, match="empty token"):
        parse_prosite("R--D", PROTEIN)


def test_prosite_refuses_an_unreadable_token():
    with pytest.raises(PatternError, match="not a PROSITE token"):
        parse_prosite("R-(2)-D", PROTEIN)


def test_prosite_refuses_an_unknown_residue():
    with pytest.raises(PatternError, match="not protein codes"):
        parse_prosite("[B1]-G", PROTEIN)


def test_prosite_refuses_an_exclusion_of_everything():
    with pytest.raises(PatternError, match="excludes every residue"):
        parse_prosite("{ACDEFGHIKLMNPQRSTVWY}-G", PROTEIN)


def test_prosite_refuses_a_backwards_range():
    with pytest.raises(PatternError, match="runs backwards"):
        parse_prosite("x(4,2)", PROTEIN)


def test_prosite_refuses_a_zero_repeat():
    with pytest.raises(PatternError, match="repeat of zero"):
        parse_prosite("x(0)", PROTEIN)


def test_prosite_refuses_an_enormous_repeat():
    with pytest.raises(PatternError, match="exceeds the limit"):
        parse_prosite(f"x({MAX_REPEAT + 1})", PROTEIN)


def test_prosite_works_for_nucleotide_too():
    positions, _, _, _, _ = parse_prosite("[AG]-x(3)-T-A-T-A", NUCLEOTIDE)
    assert len(positions) == 8
    assert positions[0] == frozenset("AG")
    assert positions[1] == frozenset("ACGT")


def test_prosite_empty_pattern_is_refused():
    with pytest.raises(PatternError, match="pattern is empty"):
        parse_prosite("<>", PROTEIN)


def test_regex_compiles_and_has_no_fixed_length():
    compiled = nucleotide("GA+TTC", syntax="regex")
    assert compiled.fixed_length is None
    assert not compiled.supports_mismatches
    assert compiled.maximum is None


def test_regex_that_cannot_be_compiled_is_refused():
    with pytest.raises(PatternError, match="not a valid regular expression"):
        nucleotide("GA(TTC", syntax="regex")


def test_a_pattern_matching_an_empty_string_is_refused():
    with pytest.raises(PatternError, match="empty string"):
        nucleotide("A*", syntax="regex")


def test_a_user_regex_is_compiled_untouched():
    # Nothing is wrapped around it, so its own groups keep their numbers and
    # its own inline flags stay at the start of the expression.
    assert not nucleotide("AC", syntax="regex").wrapped
    assert nucleotide("([AC])\\1", syntax="regex").regex.match("AA")
    assert nucleotide("(?i)ac", syntax="regex").regex.match("ac")
    assert nucleotide("(?P<x>A)(?P=x)", syntax="regex").regex.match("AA")
    assert nucleotide(f"(?P<{GROUP}>AC)", syntax="regex").regex.match("AC")


def test_a_pattern_this_package_builds_carries_its_own_group():
    compiled = nucleotide("ACGT")
    assert compiled.wrapped
    assert compiled.matched(compiled.regex.search("ACGT")) == "ACGT"


def test_the_compiled_regex_finds_overlapping_matches():
    compiled = nucleotide("AA")
    assert [match.start() for match in compiled.regex.finditer("AAAA")] == [0, 1, 2]


def test_compiling_against_the_wrong_alphabet_is_refused():
    with pytest.raises(PatternError, match="cannot be searched in"):
        compile_motif(
            Motif(name="m", pattern="R-G-D", syntax="prosite", alphabet="protein"),
            NUCLEOTIDE,
        )


def test_errors_are_prefixed_with_the_motif_name():
    with pytest.raises(PatternError, match="^kozak:"):
        compile_motif(Motif(name="kozak", pattern="GCCQCC"))


def test_the_character_class_admits_ambiguity_codes_it_covers():
    # A position permitting every base must admit a sequence N, which is not
    # one of A, C, G or T.
    assert nucleotide("N").regex.search("N")
    assert not nucleotide("A").regex.search("N")


def test_named_reads_a_bare_pattern_as_its_own_name():
    motif = named("TATAWAWR", "iupac", "nucleotide", 1)
    assert motif.name == "TATAWAWR"
    assert motif.pattern == "TATAWAWR"


def test_named_reads_name_equals_pattern():
    motif = named("kozak=GCCRCCATGG", "iupac", "nucleotide", 1)
    assert motif.name == "kozak"
    assert motif.pattern == "GCCRCCATGG"


def test_named_names_a_long_pattern_after_its_position():
    motif = named("A" * 40, "regex", "nucleotide", 3)
    assert motif.name == "regex-3"


def test_named_refuses_an_empty_name_or_pattern():
    with pytest.raises(PatternError, match="name before"):
        named("=ACGT", "iupac", "nucleotide", 1)
    with pytest.raises(PatternError, match="pattern after"):
        named("x=", "iupac", "nucleotide", 1)
    with pytest.raises(PatternError, match="pattern is empty"):
        named("  ", "iupac", "nucleotide", 1)
