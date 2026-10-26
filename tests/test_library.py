"""The built-in catalogue and motif files."""

import io

import pytest

from motif_scan import library
from motif_scan.patterns import PatternError, compile_motif


def test_every_catalogue_motif_compiles():
    for motif in library.all_motifs():
        compile_motif(motif)


def test_every_catalogue_motif_cites_a_source():
    for motif in library.all_motifs():
        assert motif.source, motif.name
        assert motif.description, motif.name


def test_every_catalogue_motif_has_a_category():
    assert all(motif.category for motif in library.all_motifs())
    assert set(library.categories()) == {
        motif.category for motif in library.all_motifs()
    }


def test_catalogue_names_are_unique():
    names = [motif.name for motif in library.all_motifs()]
    assert len(names) == len(set(names))


def test_both_alphabets_are_represented():
    alphabets = {motif.alphabet for motif in library.all_motifs()}
    assert alphabets == {"nucleotide", "protein"}


def test_every_protein_motif_is_written_in_prosite_syntax():
    for motif in library.all_motifs():
        if motif.alphabet == "protein":
            assert motif.syntax == "prosite", motif.name


def test_get_is_case_insensitive():
    assert library.get("tata").name == "TATA"
    assert library.get("  Kozak ").name == "Kozak"


def test_get_suggests_a_near_name():
    with pytest.raises(KeyError, match="did you mean"):
        library.get("kozak-")


def test_get_refuses_an_unknown_name():
    with pytest.raises(KeyError, match="not a motif in the catalogue"):
        library.get("zzz")


def test_resolve_keeps_the_order_asked_for():
    chosen = library.resolve(["polyA", "TATA"])
    assert [motif.name for motif in chosen] == ["polyA", "TATA"]


def test_resolve_applies_a_repeated_name_once():
    assert len(library.resolve(["TATA", "tata"])) == 1


def test_select_filters_by_category_and_alphabet():
    promoters = library.select(category="core-promoter")
    assert promoters
    assert all(motif.category == "core-promoter" for motif in promoters)
    proteins = library.select(alphabet="protein")
    assert all(motif.alphabet == "protein" for motif in proteins)


def test_select_searches_name_description_and_pattern():
    assert library.select(search="kozak")
    assert library.select(search="glycosylation")
    assert library.select(search="TATAWAWR")
    assert library.select(search="zzz") == []


def test_the_aav_elements_are_present():
    names = {motif.name for motif in library.all_motifs()}
    assert {"AAV-RBE", "AAV-trs", "PLA2"} <= names


def test_reading_a_motif_file():
    text = (
        "# a comment\n"
        "name\tpattern\tsyntax\talphabet\n"
        "mine\tGGWCC\tiupac\tnucleotide\n"
        "yours\tR-G-D\tprosite\tprotein\n"
    )
    motifs = library.read_motifs(io.StringIO(text))
    assert [motif.name for motif in motifs] == ["mine", "yours"]
    assert motifs[1].alphabet == "protein"
    assert motifs[0].category == "file"


def test_a_motif_file_needs_only_a_name_and_a_pattern():
    motifs = library.read_motifs(io.StringIO("name\tpattern\nmine\tGGWCC\n"))
    assert motifs[0].syntax == "iupac"
    assert motifs[0].alphabet == "nucleotide"


def test_a_motif_file_without_the_required_columns_is_refused():
    # Columns that exist but are not enough.
    with pytest.raises(PatternError, match="must name at least"):
        library.read_motifs(io.StringIO("syntax\tdescription\niupac\tx\n"))
    # Columns that do not exist at all are named instead.
    with pytest.raises(PatternError, match="unknown column"):
        library.read_motifs(io.StringIO("motif\tseq\nmine\tGGWCC\n"))


def test_a_motif_file_with_an_unknown_column_is_refused():
    with pytest.raises(PatternError, match="unknown column"):
        library.read_motifs(io.StringIO("name\tpattern\tcolour\nm\tGG\tred\n"))


def test_a_motif_file_row_without_a_pattern_is_refused():
    with pytest.raises(PatternError, match="both a name and a pattern"):
        library.read_motifs(io.StringIO("name\tpattern\nmine\t\n"))


def test_a_motif_file_with_an_unknown_syntax_is_refused():
    with pytest.raises(PatternError, match="not a pattern syntax"):
        library.read_motifs(io.StringIO("name\tpattern\tsyntax\nm\tGG\tglob\n"))


def test_a_motif_file_with_an_unknown_alphabet_is_refused():
    with pytest.raises(PatternError, match="not an alphabet"):
        library.read_motifs(io.StringIO("name\tpattern\talphabet\nm\tGG\trna\n"))


def test_an_empty_motif_file_is_refused():
    with pytest.raises(PatternError, match="holds no motifs"):
        library.read_motifs(io.StringIO("\n# only a comment\n"))


def test_a_header_only_motif_file_is_refused():
    with pytest.raises(PatternError, match="header only"):
        library.read_motifs(io.StringIO("name\tpattern\n"))


def test_load_reads_a_file(tmp_path):
    path = tmp_path / "m.tsv"
    path.write_text("name\tpattern\nmine\tGGWCC\n")
    assert library.load(path)[0].pattern == "GGWCC"


def test_unique_refuses_a_repeated_name():
    motifs = library.resolve(["TATA"]) + [
        library.Motif(name="TATA", pattern="TTTT")
    ]
    with pytest.raises(PatternError, match="more than one motif is named"):
        library.unique(motifs)


def test_unique_keeps_an_identical_repeat_once():
    motifs = library.resolve(["TATA"]) * 2
    assert len(library.unique(motifs)) == 1


def test_the_bundled_motif_file_reads_and_compiles(examples):
    motifs = library.load(examples / "payload_motifs.tsv")
    assert motifs
    for motif in motifs:
        compile_motif(motif)


def test_the_bundled_loxp_pattern_matches_the_canonical_site(examples):
    from motif_scan.search import find

    motifs = {motif.name: motif for motif in library.load(examples / "payload_motifs.tsv")}
    canonical = "ATAACTTCGTATAATGTATGCTATACGAAGTTAT"
    (hit,) = find(canonical, compile_motif(motifs["loxP"]))
    assert (hit.start, hit.end) == (1, len(canonical))


def test_the_readme_states_the_catalogue_size():
    # The figure in the README is a claim about this file, and drifts silently
    # when a motif is added.
    from pathlib import Path

    readme = Path(__file__).resolve().parent.parent / "README.md"
    if not readme.is_file():
        pytest.skip("the README is not present in this build")
    text = readme.read_text(encoding="utf-8")
    counts = (
        len(library.all_motifs()),
        len(library.select(alphabet="nucleotide")),
        len(library.select(alphabet="protein")),
    )
    assert (
        f"{counts[0]} regulatory and protein elements" in text
    ), f"the README should state {counts[0]} elements"
    assert f"{counts[1]} nucleotide, {counts[2]} protein" in text


def test_a_repeated_column_name_is_refused():
    # csv.DictReader keeps the last of two columns of one name, which would
    # rename every motif in the file without saying so.
    with pytest.raises(PatternError, match="names a column more than once"):
        library.read_motifs(io.StringIO("name\tpattern\tname\nm\tAA\tx\n"))


def test_a_row_with_too_many_fields_is_refused():
    with pytest.raises(PatternError, match="3 field\\(s\\) against 2 column"):
        library.read_motifs(io.StringIO("name\tpattern\nm\tAA\textra\n"))


def test_a_row_with_too_few_fields_is_refused():
    with pytest.raises(PatternError, match="1 field\\(s\\) against 2 column"):
        library.read_motifs(io.StringIO("name\tpattern\nm\n"))


def test_the_line_number_of_a_malformed_row_is_reported():
    with pytest.raises(PatternError, match="line 3"):
        library.read_motifs(io.StringIO("name\tpattern\nfine\tAA\nbroken\tAA\tx\n"))


def test_a_rectangular_file_with_empty_optional_cells_is_accepted():
    motifs = library.read_motifs(
        io.StringIO("name\tpattern\tcategory\nm\tAA\t\n")
    )
    assert motifs[0].category == ""
