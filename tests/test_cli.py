"""The command line."""

import json

import pytest

from motif_scan import __version__
from motif_scan.cli import build_parser, ineffective_options, main

DNA = ">seq a test record\nGCCACCATGGAATTCAAATATAAAAGGGCTCGCTCAAATAAACCC\n"
PEPTIDE = ">pep a peptide\nMKQLEDKVEELLSKNYHLENEVARLKKLVGERGDHDKKYKKKRKV\n"


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def run(argv, capsys):
    status = main(argv)
    captured = capsys.readouterr()
    return status, captured.out, captured.err


def test_version(capsys):
    with pytest.raises(SystemExit) as exit:
        main(["--version"])
    assert exit.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_a_command_is_required():
    with pytest.raises(SystemExit):
        main([])


def test_scan_reports_hits(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, out, _ = run(["scan", "-i", path, "-m", "polyA,Kozak"], capsys)
    assert status == 0
    assert "polyA" in out
    assert "AATAAA" in out


def test_scan_needs_a_motif(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, _, err = run(["scan", "-i", path], capsys)
    assert status == 1
    assert "no motif was asked for" in err


def test_scan_rejects_an_unknown_motif(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, _, err = run(["scan", "-i", path, "-m", "zzz"], capsys)
    assert status == 1
    assert "not a motif in the catalogue" in err


def test_scan_accepts_a_pattern_on_the_command_line(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, out, _ = run(["scan", "-i", path, "-p", "kozak=GCCRCCATGG"], capsys)
    assert status == 0
    assert "kozak" in out


def test_an_unnamed_pattern_is_named_after_itself(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    _, out, _ = run(["scan", "-i", path, "-p", "GGGCTCGCTC"], capsys)
    assert "GGGCTCGCTC" in out


def test_scan_accepts_a_prosite_pattern(tmp_path, capsys):
    path = write(tmp_path, "p.fasta", PEPTIDE)
    status, out, _ = run(["scan", "-i", path, "--prosite", "rgd=R-G-D"], capsys)
    assert status == 0
    assert "rgd" in out


def test_scan_accepts_a_regular_expression(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, out, _ = run(["scan", "-i", path, "--regex", "runs=A{3,}"], capsys)
    assert status == 0
    assert "runs" in out


def test_a_regex_matching_an_empty_string_is_refused(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, _, err = run(["scan", "-i", path, "--regex", "A*"], capsys)
    assert status == 1
    assert "empty string" in err


def test_mismatches_against_a_regex_are_refused(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, _, err = run(
        ["scan", "-i", path, "--regex", "A{3}", "--mismatches", "1"], capsys
    )
    assert status == 1
    assert "variable length" in err


def test_mismatches_find_a_variant_consensus(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    _, exact, _ = run(["scan", "-i", path, "-m", "TATA", "-f", "tsv"], capsys)
    _, tolerant, _ = run(
        ["scan", "-i", path, "-m", "TATA", "--mismatches", "1", "-f", "tsv"], capsys
    )
    assert len(tolerant.splitlines()) > len(exact.splitlines())


def test_all_applies_the_catalogue_of_the_right_alphabet(tmp_path, capsys):
    path = write(tmp_path, "p.fasta", PEPTIDE)
    status, out, _ = run(["scan", "-i", path, "--all", "--limit", "0"], capsys)
    assert status == 0
    assert "N-glycosylation" in out
    assert "TATA" not in out


def test_category_filters_the_catalogue(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, out, _ = run(
        ["scan", "-i", path, "--all", "--category", "core-promoter"], capsys
    )
    assert status == 0
    assert "TATA" in out
    assert "polyA" not in out


def test_category_without_all_is_refused(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, _, err = run(
        ["scan", "-i", path, "-m", "TATA", "--category", "core-promoter"], capsys
    )
    assert status == 1
    assert "--category applies to --all" in err


def test_an_unknown_category_is_refused(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, _, err = run(["scan", "-i", path, "--all", "--category", "zzz"], capsys)
    assert status == 1
    assert "not a category" in err


def test_a_protein_motif_against_nucleotide_is_refused(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, _, err = run(["scan", "-i", path, "-m", "RGD"], capsys)
    assert status == 1
    assert "these motifs are not" in err


def test_the_alphabet_is_reported_not_assumed(tmp_path, capsys):
    path = write(tmp_path, "p.fasta", PEPTIDE)
    _, out, _ = run(["scan", "-i", path, "--prosite", "R-G-D"], capsys)
    assert "protein" in out


def test_the_alphabet_can_be_forced(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", ">x\nGGGGGGGGGGGG\n")
    status, out, _ = run(
        ["scan", "-i", path, "--type", "protein", "--prosite", "G(3)"], capsys
    )
    assert status == 0
    assert "protein" in out


def test_mixed_alphabets_are_refused(tmp_path, capsys):
    path = write(tmp_path, "m.fasta", DNA + PEPTIDE)
    status, _, err = run(["scan", "-i", path, "-m", "polyA"], capsys)
    assert status == 1
    assert "not all the same kind of sequence" in err


def test_circular_is_refused_for_protein(tmp_path, capsys):
    path = write(tmp_path, "p.fasta", PEPTIDE)
    status, _, err = run(
        ["scan", "-i", path, "--prosite", "R-G-D", "--circular"], capsys
    )
    assert status == 1
    assert "closed molecule of nucleic acid" in err


def test_strand_can_be_restricted(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    _, both, _ = run(["scan", "-i", path, "-p", "GCTCGCTC", "-f", "tsv"], capsys)
    _, forward, _ = run(
        ["scan", "-i", path, "-p", "GCTCGCTC", "--strand", "forward", "-f", "tsv"],
        capsys,
    )
    assert len(both.splitlines()) >= len(forward.splitlines())


def test_tsv_output_has_the_documented_header(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    _, out, _ = run(["scan", "-i", path, "-m", "polyA", "-f", "tsv"], capsys)
    assert out.splitlines()[0].split("\t")[:5] == [
        "record",
        "motif",
        "start",
        "end",
        "strand",
    ]


def test_json_output_is_valid_and_complete(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    _, out, _ = run(["scan", "-i", path, "-m", "polyA,TATA", "-f", "json"], capsys)
    data = json.loads(out)
    assert data["motif_scan"] == __version__
    assert data["alphabet"] == "nucleotide"
    assert [motif["name"] for motif in data["motifs"]] == ["polyA", "TATA"]
    record = data["records"][0]
    assert record["identifier"] == "seq"
    assert record["counts"][0]["motif"] == "polyA"
    assert record["hits"][0]["start"] >= 1


def test_format_none_writes_nothing(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, out, err = run(["scan", "-i", path, "-m", "polyA", "-f", "none"], capsys)
    assert status == 0
    assert out == ""
    assert "writes nothing" in err


def test_output_writes_every_hit(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    hits = tmp_path / "hits.tsv"
    run(
        ["scan", "-i", path, "-m", "polyA", "-f", "none", "-o", str(hits)], capsys
    )
    lines = hits.read_text().splitlines()
    assert lines[0].startswith("record\tmotif")
    assert len(lines) > 1


def test_counts_file_carries_the_expectation(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    counts = tmp_path / "counts.tsv"
    run(
        ["scan", "-i", path, "-m", "polyA", "-f", "none", "--counts", str(counts)],
        capsys,
    )
    header, row = counts.read_text().splitlines()[:2]
    fields = dict(zip(header.split("\t"), row.split("\t")))
    assert fields["motif"] == "polyA"
    assert float(fields["expected"]) > 0
    assert fields["uncertainty_assessed"] == "yes"


def test_report_is_written(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    report = tmp_path / "r.html"
    run(
        ["scan", "-i", path, "-m", "polyA", "-f", "none", "--report", str(report)],
        capsys,
    )
    html = report.read_text()
    assert html.startswith("<!DOCTYPE html>")
    assert "<script" not in html


def test_output_may_not_overwrite_an_input(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, _, err = run(["scan", "-i", path, "-m", "polyA", "-o", path], capsys)
    assert status == 1
    assert "would overwrite an input file" in err
    assert path.endswith("s.fasta")
    assert "GCCACC" in open(path).read()


def test_two_outputs_may_not_share_a_file(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    same = str(tmp_path / "both.tsv")
    status, _, err = run(
        ["scan", "-i", path, "-m", "polyA", "-o", same, "--counts", same], capsys
    )
    assert status == 1
    assert "same file" in err


def test_several_inputs_are_labelled_by_file(tmp_path, capsys):
    first = write(tmp_path, "a.fasta", DNA)
    second = write(tmp_path, "b.fasta", DNA)
    _, out, _ = run(
        ["scan", "-i", first, "-i", second, "-m", "polyA", "-f", "tsv"], capsys
    )
    records = {line.split("\t")[0] for line in out.splitlines()[1:]}
    assert len(records) == 2
    assert all(":" in record for record in records)


def test_a_missing_file_is_reported(tmp_path, capsys):
    status, _, err = run(
        ["scan", "-i", str(tmp_path / "nope.fasta"), "-m", "polyA"], capsys
    )
    assert status == 1
    assert "no such file" in err


def test_motifs_file_is_read(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    motifs = write(tmp_path, "m.tsv", "name\tpattern\nmine\tAATAAA\n")
    status, out, _ = run(
        ["scan", "-i", path, "--motifs-file", motifs], capsys
    )
    assert status == 0
    assert "mine" in out


def test_two_motifs_of_the_same_name_are_refused(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, _, err = run(
        ["scan", "-i", path, "-m", "polyA", "-p", "polyA=AATAAA"], capsys
    )
    assert status == 1
    assert "more than one motif is named" in err


def test_sequence_view_prints_the_residues(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    _, out, _ = run(["scan", "-i", path, "-m", "polyA", "--sequence"], capsys)
    assert "GCCACCATGG" in out


def test_sequence_view_respects_the_highlight_limit(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    _, out, _ = run(
        ["scan", "-i", path, "-m", "polyA", "--sequence", "--highlight-limit", "5"],
        capsys,
    )
    assert "not shown" in out


def test_the_map_can_be_switched_off(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    _, with_map, _ = run(["scan", "-i", path, "-m", "polyA"], capsys)
    _, without, _ = run(["scan", "-i", path, "-m", "polyA", "--no-map"], capsys)
    assert len(without.splitlines()) < len(with_map.splitlines())


def test_limit_truncates_the_listing_and_says_so(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    _, out, _ = run(
        ["scan", "-i", path, "-p", "AA", "--limit", "1"], capsys
    )
    assert "not shown" in out


def test_sort_by_motif(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    _, out, _ = run(
        ["scan", "-i", path, "-m", "polyA,TATA", "--sort", "motif", "-f", "tsv"],
        capsys,
    )
    names = [line.split("\t")[1] for line in out.splitlines()[1:]]
    assert names == sorted(names, key=["polyA", "TATA"].index)


def test_uncertain_positions_are_counted_and_reported(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", ">x\nAATANAAAAACCC\n")
    _, out, _ = run(["scan", "-i", path, "-m", "polyA"], capsys)
    assert "uncertain" in out


def test_uncertain_windows_can_be_listed(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", ">x\nAATANACCC\n")
    _, plain, _ = run(["scan", "-i", path, "-m", "polyA", "-f", "tsv"], capsys)
    _, listed, _ = run(
        ["scan", "-i", path, "-m", "polyA", "--include-uncertain", "-f", "tsv"],
        capsys,
    )
    assert len(listed.splitlines()) > len(plain.splitlines())


def test_negative_numbers_are_refused():
    with pytest.raises(SystemExit):
        main(["scan", "-i", "x", "-m", "polyA", "--mismatches", "-1"])


def test_a_zero_width_map_is_refused():
    with pytest.raises(SystemExit):
        main(["scan", "-i", "x", "-m", "polyA", "--map-width", "0"])


def test_motifs_lists_the_catalogue(capsys):
    status, out, _ = run(["motifs"], capsys)
    assert status == 0
    assert "TATA" in out
    assert "N-glycosylation" in out


def test_motifs_filters_by_alphabet(capsys):
    _, out, _ = run(["motifs", "--type", "protein"], capsys)
    assert "N-glycosylation" in out
    assert "TATA" not in out


def test_motifs_filters_by_category_and_search(capsys):
    _, out, _ = run(["motifs", "--category", "vector"], capsys)
    assert "AAV-RBE" in out
    _, found, _ = run(["motifs", "--search", "kozak"], capsys)
    assert "Kozak" in found


def test_motifs_refuses_filters_that_match_nothing(capsys):
    status, _, err = run(["motifs", "--search", "zzz"], capsys)
    assert status == 1
    assert "no catalogue motif matches" in err


def test_motifs_json_round_trips(capsys):
    _, out, _ = run(["motifs", "--type", "protein", "-f", "json"], capsys)
    data = json.loads(out)
    assert any(entry["name"] == "RGD" for entry in data)
    assert all(entry["source"] for entry in data)


def test_motifs_tsv_can_be_read_back_as_a_motif_file(tmp_path, capsys):
    from motif_scan import library

    out = tmp_path / "cat.tsv"
    run(["motifs", "-f", "none", "-o", str(out)], capsys)
    # The 'length' column is informational, so it is dropped on the way back in.
    lines = out.read_text().splitlines()
    header = lines[0].split("\t")
    keep = [index for index, name in enumerate(header) if name != "length"]
    trimmed = "\n".join(
        "\t".join(line.split("\t")[index] for index in keep) for line in lines
    )
    motifs = library.read_motifs(trimmed.splitlines(keepends=True).__iter__())
    assert len(motifs) == len(library.all_motifs())


def test_codes_shows_both_alphabets_and_the_notations(capsys):
    status, out, _ = run(["codes"], capsys)
    assert status == 0
    assert "nucleotide" in out and "protein" in out
    assert "Pattern notations" in out
    assert "prosite" in out


def test_codes_filters_by_alphabet(capsys):
    _, out, _ = run(["codes", "--type", "nucleotide", "-f", "tsv"], capsys)
    assert "protein" not in out


def test_codes_json(capsys):
    _, out, _ = run(["codes", "-f", "json"], capsys)
    data = json.loads(out)
    assert {"alphabet", "code", "stands_for", "residues"} == set(data[0])


def test_ineffective_options_are_reported(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    _, _, err = run(
        ["scan", "-i", path, "-m", "polyA", "-f", "tsv", "--limit", "3"], capsys
    )
    assert "--limit affects only" in err


def test_ineffective_options_are_silent_when_everything_applies():
    parser = build_parser()
    args = parser.parse_args(["scan", "-i", "x", "-m", "polyA", "--limit", "3"])
    assert ineffective_options(args) == []


def test_bins_only_affect_the_report():
    parser = build_parser()
    args = parser.parse_args(["scan", "-i", "x", "-m", "polyA", "--bins", "5"])
    assert any("--bins" in note for note in ineffective_options(args))


# --- output path guard ----------------------------------------------------


def test_output_may_not_overwrite_the_motif_file(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    motifs = write(tmp_path, "m.tsv", "name\tpattern\nmine\tAATAAA\n")
    original = open(motifs).read()
    status, _, err = run(
        ["scan", "-i", path, "--motifs-file", motifs, "-f", "none", "-o", motifs],
        capsys,
    )
    assert status == 1
    assert "would overwrite an input file" in err
    assert open(motifs).read() == original


def test_the_report_may_not_overwrite_the_motif_file(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    motifs = write(tmp_path, "m.tsv", "name\tpattern\nmine\tAATAAA\n")
    status, _, err = run(
        [
            "scan", "-i", path, "--motifs-file", motifs, "-f", "none",
            "--report", motifs,
        ],
        capsys,
    )
    assert status == 1
    assert "would overwrite an input file" in err


def test_two_outputs_may_not_share_a_file_written_differently(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    plain = str(tmp_path / "out.tsv")
    dotted = str(tmp_path / "." / "out.tsv")
    status, _, err = run(
        ["scan", "-i", path, "-m", "polyA", "-f", "none", "-o", plain,
         "--counts", dotted],
        capsys,
    )
    assert status == 1
    assert "written to the same file" in err
    assert not (tmp_path / "out.tsv").exists()


def test_an_output_hard_linked_to_an_input_is_refused(tmp_path, capsys):
    import os

    path = write(tmp_path, "s.fasta", DNA)
    link = str(tmp_path / "linked.tsv")
    os.link(path, link)
    status, _, err = run(
        ["scan", "-i", path, "-m", "polyA", "-f", "none", "-o", link], capsys
    )
    assert status == 1
    assert "would overwrite an input file" in err
    assert "GCCACC" in open(path).read()


def test_outputs_to_different_files_are_allowed(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, _, _ = run(
        [
            "scan", "-i", path, "-m", "polyA", "-f", "none",
            "-o", str(tmp_path / "hits.tsv"),
            "--counts", str(tmp_path / "counts.tsv"),
            "--report", str(tmp_path / "r.html"),
        ],
        capsys,
    )
    assert status == 0
    assert (tmp_path / "hits.tsv").exists()
    assert (tmp_path / "counts.tsv").exists()
    assert (tmp_path / "r.html").exists()


def test_standard_input_and_output_are_not_treated_as_a_collision(tmp_path):
    from motif_scan.cli import _guard_outputs, build_parser

    args = build_parser().parse_args(
        ["scan", "-i", "-", "-m", "polyA", "-o", "-", "--counts", "-"]
    )
    _guard_outputs(args)


# --- duplicate record identifiers -----------------------------------------


def test_duplicate_identifiers_get_their_own_sequence_in_the_report(
    tmp_path, capsys
):
    import re

    path = write(tmp_path, "d.fasta", ">dup\nAAAAAAAAAA\n>dup\nCCCCCCCCCC\n")
    report = tmp_path / "d.html"
    status, _, _ = run(
        ["scan", "-i", path, "-p", "AA", "-f", "none", "--report", str(report)],
        capsys,
    )
    assert status == 0
    html = report.read_text()
    drawn = [
        re.sub(r"<[^>]+>", "", block.split("</div>")[0]).strip()
        for block in html.split('<div class="seqmap">')[1:]
    ]
    # Only the first record has AA hits, and the residues shown beside them
    # must be its own.
    assert drawn
    assert all("C" not in text for text in drawn)
    assert "AAAAAAAAAA" in drawn[0]


def test_duplicate_identifiers_keep_their_own_residues_in_the_sequence_view(
    tmp_path, capsys
):
    path = write(tmp_path, "d.fasta", ">dup\nAAAAAAAAAA\n>dup\nCCCCCCCCCC\n")
    _, out, _ = run(
        ["scan", "-i", path, "-p", "AA", "--sequence", "--no-map"], capsys
    )
    assert "AAAAAAAAAA" in out


# --- overlap suppression ---------------------------------------------------


def test_no_overlapping_withholds_the_expectation(tmp_path, capsys):
    path = write(tmp_path, "p.fasta", ">p\nAAAAA\n")
    counts = tmp_path / "k.tsv"
    run(
        ["scan", "-i", path, "--type", "protein", "-p", "AA", "-f", "none",
         "--no-overlapping", "--counts", str(counts)],
        capsys,
    )
    header, row = counts.read_text().splitlines()[:2]
    fields = dict(zip(header.split("\t"), row.split("\t")))
    assert fields["hits"] == "2"
    assert fields["expected"] == ""
    assert fields["observed_over_expected"] == ""


def test_overlapping_hits_still_get_an_expectation(tmp_path, capsys):
    path = write(tmp_path, "p.fasta", ">p\nAAAAA\n")
    counts = tmp_path / "k.tsv"
    run(
        ["scan", "-i", path, "--type", "protein", "-p", "AA", "-f", "none",
         "--counts", str(counts)],
        capsys,
    )
    header, row = counts.read_text().splitlines()[:2]
    fields = dict(zip(header.split("\t"), row.split("\t")))
    assert fields["hits"] == "4"
    assert float(fields["expected"]) > 0


def test_the_terminal_says_why_the_expectation_is_missing(tmp_path, capsys):
    path = write(tmp_path, "p.fasta", ">p\nAAAAA\n")
    _, out, _ = run(
        ["scan", "-i", path, "--type", "protein", "-p", "AA", "--no-overlapping"],
        capsys,
    )
    assert "no expected count is reported" in out


def test_the_report_says_why_the_expectation_is_missing(tmp_path, capsys):
    path = write(tmp_path, "p.fasta", ">p\nAAAAA\n")
    report = tmp_path / "r.html"
    run(
        ["scan", "-i", path, "--type", "protein", "-p", "AA", "-f", "none",
         "--no-overlapping", "--report", str(report)],
        capsys,
    )
    assert "No expectation is reported at all" in report.read_text()


# --- the expectation basis -------------------------------------------------


def test_ambiguity_does_not_inflate_the_expectation(tmp_path, capsys):
    path = write(tmp_path, "a.fasta", ">a\nANNNNNNNNN\n")
    counts = tmp_path / "k.tsv"
    run(
        ["scan", "-i", path, "-p", "A", "--strand", "forward", "-f", "none",
         "--counts", str(counts)],
        capsys,
    )
    header, row = counts.read_text().splitlines()[:2]
    fields = dict(zip(header.split("\t"), row.split("\t")))
    assert fields["hits"] == "1"
    assert fields["uncertain"] == "9"
    assert float(fields["expected"]) == pytest.approx(1.0)
    assert float(fields["observed_over_expected"]) == pytest.approx(1.0)


def test_including_uncertain_windows_counts_them_on_both_sides(tmp_path, capsys):
    path = write(tmp_path, "a.fasta", ">a\nANNNNNNNNN\n")
    counts = tmp_path / "k.tsv"
    run(
        ["scan", "-i", path, "-p", "A", "--strand", "forward", "-f", "none",
         "--include-uncertain", "--counts", str(counts)],
        capsys,
    )
    header, row = counts.read_text().splitlines()[:2]
    fields = dict(zip(header.split("\t"), row.split("\t")))
    assert fields["hits"] == "10"
    assert float(fields["expected"]) == pytest.approx(10.0)


def test_a_clean_record_counts_every_window(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    counts = tmp_path / "k.tsv"
    run(
        ["scan", "-i", path, "-m", "polyA", "-f", "none", "--counts", str(counts)],
        capsys,
    )
    row = dict(
        zip(*(line.split("\t") for line in counts.read_text().splitlines()[:2]))
    )
    assert float(row["expected"]) > 0


# --- uncertain hits --------------------------------------------------------


def test_an_uncertain_hit_is_not_reported_as_exact(tmp_path, capsys):
    path = write(tmp_path, "u.fasta", ">u\nAATANACCC\n")
    counts = tmp_path / "k.tsv"
    run(
        ["scan", "-i", path, "-m", "polyA", "--include-uncertain", "-f", "none",
         "--counts", str(counts)],
        capsys,
    )
    fields = dict(
        zip(*(line.split("\t") for line in counts.read_text().splitlines()[:2]))
    )
    assert fields["hits"] == "1"
    assert fields["exact"] == "0"
    assert fields["with_mismatches"] == "0"
    assert fields["listed_uncertain"] == "1"


def test_the_report_does_not_claim_uncertain_windows_are_excluded(
    tmp_path, capsys
):
    path = write(tmp_path, "u.fasta", ">u\nAATANACCC\n")
    report = tmp_path / "u.html"
    run(
        ["scan", "-i", path, "-m", "polyA", "--include-uncertain", "-f", "none",
         "--report", str(report)],
        capsys,
    )
    html = report.read_text()
    assert "listed among the hits below" in html
    assert "not counted as hits" not in html
    assert "uncertain" in html


def test_the_report_says_uncertain_windows_are_excluded_by_default(
    tmp_path, capsys
):
    path = write(tmp_path, "u.fasta", ">u\nAATANACCC\n")
    report = tmp_path / "u.html"
    run(
        ["scan", "-i", path, "-m", "polyA", "-f", "none", "--report", str(report)],
        capsys,
    )
    assert "not counted as hits" in report.read_text()


def test_the_terminal_says_what_became_of_uncertain_windows(tmp_path, capsys):
    path = write(tmp_path, "u.fasta", ">u\nAATANACCC\n")
    _, plain, _ = run(["scan", "-i", path, "-m", "polyA", "--no-map"], capsys)
    assert "not counted as hits" in plain
    _, listed, _ = run(
        ["scan", "-i", path, "-m", "polyA", "--include-uncertain", "--no-map"],
        capsys,
    )
    assert "listed among the hits below" in listed


# --- regular expressions on the command line -------------------------------


def test_a_regex_with_an_equals_sign_is_not_split_into_a_name(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    status, out, _ = run(
        ["scan", "-i", path, "--regex", "(?P<x>A)(?P=x)", "-f", "tsv"], capsys
    )
    assert status == 0
    assert "AA" in out


def test_a_named_regex_still_takes_its_name(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    _, out, _ = run(
        ["scan", "-i", path, "--regex", "runs=A{3}", "-f", "tsv"], capsys
    )
    assert "runs" in out


def test_regex_features_python_accepts_are_accepted(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    for pattern in (r"([AC])\1", "(?i)acc", "AC(?=C)"):
        status, _, err = run(
            ["scan", "-i", path, "--regex", pattern, "-f", "none"], capsys
        )
        assert status == 0, (pattern, err)


# --- motif files -----------------------------------------------------------


def test_a_repeated_column_in_a_motif_file_is_refused(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    motifs = write(tmp_path, "m.tsv", "name\tpattern\tname\nm\tAA\tx\n")
    status, _, err = run(
        ["scan", "-i", path, "--motifs-file", motifs, "-f", "none"], capsys
    )
    assert status == 1
    assert "names a column more than once" in err


def test_a_ragged_motif_file_is_refused(tmp_path, capsys):
    path = write(tmp_path, "s.fasta", DNA)
    motifs = write(tmp_path, "m.tsv", "name\tpattern\nm\tAA\textra\n")
    status, _, err = run(
        ["scan", "-i", path, "--motifs-file", motifs, "-f", "none"], capsys
    )
    assert status == 1
    assert "against 2 column" in err
