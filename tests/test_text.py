"""Formatting helpers."""

from motif_scan import text


def test_thousands_groups():
    assert text.thousands(4679) == "4,679"


def test_plural_agrees_the_noun():
    assert text.plural(1, "hit") == "1 hit"
    assert text.plural(2, "hit") == "2 hits"
    assert text.plural(0, "hit") == "0 hits"


def test_plural_takes_an_irregular_suffix():
    assert text.plural(2, "analys", "es") == "2 analyses"


def test_number_formats_or_reports_missing():
    assert text.number(1.5) == "1.50"
    assert text.number(None) == "n/a"
    assert text.number(None, missing="-") == "-"


def test_ratio_keeps_a_small_value_informative():
    assert text.ratio(1.234) == "1.23"
    assert text.ratio(0.004) == "4e-03"
    assert text.ratio(None) == "n/a"
    assert text.ratio(0.0) == "0.00"


def test_residues_uses_the_unit_given():
    assert text.residues(735, "aa") == "735 aa"


def test_span_writes_a_range():
    assert text.span(1, 34) == "1..34"


def test_truncate_marks_what_it_shortened():
    assert text.truncate("abcdef", 4) == "abc\u2026"
    assert text.truncate("abc", 4) == "abc"
    assert text.truncate("abc", 0) == "abc"


def test_wrap_sequence_keeps_every_residue():
    lines = text.wrap_sequence("ACGT" * 40, 60)
    assert "".join(lines) == "ACGT" * 40
    assert len(lines[0]) == 60
