"""The self-contained HTML report."""

from collections import Counter

from motif_scan import __version__
from motif_scan.patterns import compile_motif
from motif_scan.report import render_html
from motif_scan.search import find
from motif_scan.summary import RecordScan, Scan, distribute, tally

SEQUENCE = "GCCACCATGGAATTCAAAAAAAAAATATAAAAGGGCTCGCTCAAATAAACCC"


def build(circular=False, motifs=("Kozak", "TATA", "polyA"), sequence=SEQUENCE,
          highlight_limit=20000, limit=500, identifier="test", description="a record"):
    from motif_scan import library

    compiled = [compile_motif(library.get(name)) for name in motifs]
    hits = []
    counts = []
    for motif in compiled:
        found = find(
            sequence, motif, mismatches=1, circular=circular, record=identifier
        )
        hits.extend(found)
        counts.append(
            tally(
                motif,
                found,
                len(sequence),
                Counter(sequence),
                mismatches=1,
                circular=circular,
            )
        )
    hits.sort(key=lambda hit: (hit.start, hit.strand, hit.motif))
    record = RecordScan(
        identifier=identifier,
        description=description,
        source="",
        residues=len(sequence),
        alphabet="nucleotide",
        circular=circular,
        ambiguous=0,
        hits=tuple(hits),
        counts=tuple(counts),
        distribution=distribute(hits, [m.name for m in compiled], len(sequence), 8),
    )
    scan = Scan(
        name="test.fasta",
        records=(record,),
        motifs=tuple(compiled),
        parameters=(("Input", "test.fasta"),),
        mismatches=1,
    )
    return render_html(
        scan,
        [sequence],
        __version__,
        width=30,
        highlight_limit=highlight_limit,
        limit=limit,
    )


def test_report_is_a_complete_document():
    html = build()
    assert html.startswith("<!DOCTYPE html>")
    assert html.rstrip().endswith("</html>")


def test_report_loads_nothing_from_the_network():
    html = build()
    assert "http://" not in html.replace("http://www.w3.org", "")
    assert "https://" not in html
    assert "src=" not in html


def test_report_uses_no_javascript():
    html = build()
    assert "<script" not in html
    assert "javascript:" not in html
    assert "onclick" not in html


def test_report_embeds_its_own_style():
    assert "<style>" in build()


def test_report_is_deterministic():
    assert build() == build()


def test_report_is_titled_after_its_subject():
    assert "<title>Motif scan of test.fasta</title>" in build()


def test_report_names_every_motif_applied():
    html = build()
    for name in ("Kozak", "TATA", "polyA"):
        assert name in html


def test_report_records_the_version():
    assert f"motif-scan {__version__}" in build()


def test_report_states_the_expectation_caveat():
    html = build()
    assert "not a test" in html
    assert "independent residues" in html


def test_report_states_the_consensus_caveat():
    html = build()
    assert "blunt instrument" in html
    assert "matrix model" in html


def test_report_cites_the_source_of_each_catalogue_motif():
    html = build()
    assert "Kozak 1987" in html
    assert "Bucher 1990" in html


def test_report_highlights_the_sequence():
    html = build()
    assert "<mark" in html
    assert 'class="seqmap"' in html


def test_report_does_not_highlight_a_record_past_the_limit():
    html = build(highlight_limit=10)
    assert "past the" in html
    assert 'class="seqmap"' not in html


def test_hits_link_into_the_sequence_map_only_when_it_is_drawn():
    assert 'href="#r0-' in build()
    assert 'href="#r0-' not in build(highlight_limit=10)


def test_report_carries_filter_controls_for_every_motif():
    html = build()
    for index in range(3):
        assert f'id="r0-m{index}"' in html
        assert f"#r0-m{index}:not(:checked) ~ table tr.m-{index}" in html


def test_report_carries_a_strand_filter_for_the_strands_present():
    html = build()
    assert 'id="r0-s-plus"' in html
    assert "#r0-s-plus:not(:checked) ~ table tr.s-plus" in html


def test_report_carries_a_mismatch_filter():
    html = build()
    assert 'id="r0-inexact"' in html
    assert "with mismatches" in html


def test_report_says_filtering_does_not_change_the_counts():
    assert "changes what is listed, not what" in build()


def test_rows_are_classed_by_motif_and_strand():
    html = build()
    assert 'class="m-0 s-' in html or 'class="m-1 s-' in html


def test_report_draws_the_distribution_chart():
    html = build()
    assert 'class="chart"' in html
    assert "equal bins along the record" in html


def test_report_draws_the_scaled_map():
    assert 'class="map"' in build()


def test_circular_report_mentions_the_join():
    assert "join of a circular" in build(circular=True)


def test_report_without_hits_says_so():
    html = build(sequence="CCCCCCCCCCCCCCCCCCCC")
    assert "No motif was found in this record" in html


def test_report_truncates_a_long_hit_list_and_says_so():
    html = build(limit=1)
    assert "Showing the first 1 hit" in html


def test_report_escapes_its_input():
    html = build(identifier="<script>x</script>", description="<b>bold</b>")
    assert "<script>x</script>" not in html
    assert "&lt;script&gt;" in html
    assert "<b>bold</b>" not in html


def test_report_reports_uncertain_positions():
    from motif_scan import library
    from motif_scan.search import count_uncertain

    compiled = compile_motif(library.get("polyA"))
    sequence = "AATANACCC"
    uncertain = count_uncertain(sequence, compiled)
    count = tally(
        compiled, [], len(sequence), Counter(sequence), uncertain=uncertain
    )
    record = RecordScan(
        "r", "", "", len(sequence), "nucleotide", False, 1,
        (), (count,), distribute([], ["polyA"], len(sequence), 4),
    )
    scan = Scan("x", (record,), (compiled,))
    html = render_html(scan, [sequence], __version__)
    assert "could not be called either way" in html
    assert "not counted as hits" in html


def test_a_multi_record_scan_gets_a_totals_table():
    from motif_scan import library

    compiled = (compile_motif(library.get("polyA")),)
    records = tuple(
        RecordScan(
            f"r{index}", "", "", 10, "nucleotide", False, 0,
            (), (tally(compiled[0], [], 10, Counter("AATAAACCCC")),),
            distribute([], ["polyA"], 10, 4),
        )
        for index in range(2)
    )
    html = render_html(Scan("x", records, compiled), [], __version__)
    assert html.count("<h2>Motifs</h2>") == 1
    assert "r0" in html and "r1" in html


def test_each_filter_input_is_immediately_followed_by_its_label():
    # The control is drawn by its label through an adjacent-sibling rule, so
    # anything between the two leaves a checkbox that cannot be seen.
    import re

    html = build()
    inputs = re.findall(r'<input type="checkbox" class="filter" id="([^"]+)" checked>'
                        r'<label for="([^"]+)">', html)
    assert inputs
    assert all(box == target for box, target in inputs)
    assert html.count('class="filter"') == len(inputs)


def test_the_controls_sit_before_the_table_they_filter():
    html = build()
    section = html.split('<section class="hits">')[1].split("</section>")[0]
    assert section.index('class="filter"') < section.index("<table>")
