"""The bundled examples, against what the literature says is in them.

These are the tests that would catch a scanner that runs but is wrong. Each
expectation comes from a published description of the sequence rather than from
a previous run of this tool, so a change in coordinates, strand convention or
degenerate matching breaks them.
"""

import pytest

from motif_scan import library
from motif_scan.patterns import compile_motif
from motif_scan.reader import read_records
from motif_scan.search import find


def scan(path, motif, **kwargs):
    (record,) = read_records(path)
    return record, find(
        record.sequence, compile_motif(library.get(motif)), record=record.identifier,
        **kwargs,
    )


def positions(hits, strand=None):
    return sorted(
        hit.start for hit in hits if strand is None or hit.strand == strand
    )


# --- human beta-globin, NG_059281.1, transcription starts at 5001 ---------

TSS = 5001


def test_hbb_ccaat_box_sits_at_minus_76(examples):
    _, hits = scan(examples / "hbb.fasta", "CCAAT")
    assert TSS - 76 in positions(hits)


def test_hbb_has_the_duplicated_caccc_box_at_minus_90_and_minus_105(examples):
    _, hits = scan(examples / "hbb.fasta", "CACCC")
    found = positions(hits)
    assert TSS - 90 in found
    assert TSS - 105 in found


def test_hbb_tata_box_is_a_variant_an_exact_search_misses(examples):
    # The beta-globin TATA box is CATAAAA, one base from the TATAWAWR
    # consensus. This is the case mismatch tolerance exists for.
    _, exact = scan(examples / "hbb.fasta", "TATA")
    promoter = [start for start in positions(exact) if TSS - 40 < start < TSS]
    assert promoter == []

    _, tolerant = scan(examples / "hbb.fasta", "TATA", mismatches=1)
    promoter = [start for start in positions(tolerant) if TSS - 40 < start < TSS]
    assert promoter == [TSS - 32]


def test_hbb_tata_variant_reads_as_published(examples):
    _, hits = scan(examples / "hbb.fasta", "TATA", mismatches=1)
    (hit,) = [hit for hit in hits if hit.start == TSS - 32]
    assert hit.matched == "CATAAAAG"
    assert hit.mismatches == 1


def test_hbb_polyadenylation_signal_lies_inside_the_gene(examples):
    record, hits = scan(examples / "hbb.fasta", "polyA")
    # The gene occupies 5001..6608; its poly(A) signal is near the 3' end.
    assert 6582 in positions(hits, strand="+")
    assert len(record) == 10106


# --- adeno-associated virus 2, NC_001401.2 --------------------------------


def test_aav2_rep_binding_element_is_in_both_inverted_terminal_repeats(examples):
    record, hits = scan(examples / "aav2.fasta", "AAV-RBE")
    found = positions(hits)
    # Two GCTC repeats resolve to two overlapping hits per arm, and the ITR is
    # itself a palindrome, so the element appears on both strands of both ends.
    assert found[:4] == [22, 26, 93, 97]
    assert found[-4:] == [4576, 4580, 4647, 4651]
    assert len(record) == 4679


def test_aav2_rbe_hits_face_both_ways(examples):
    _, hits = scan(examples / "aav2.fasta", "AAV-RBE")
    assert positions(hits, strand="+") == [22, 26, 4576, 4580]
    assert positions(hits, strand="-") == [93, 97, 4647, 4651]


def test_aav2_terminal_resolution_site_is_at_each_end(examples):
    _, hits = scan(examples / "aav2.fasta", "AAV-trs")
    assert positions(hits) == [120, 4554]
    assert positions(hits, strand="-") == [120]


def test_the_trs_is_only_required_on_the_nicked_strand(examples):
    # Brister and Muzyczka 1999: the minimal site is 3'-CCGGT/TG-5', needed
    # only on the strand that is nicked, which is why it is found once per end
    # rather than on both strands of each.
    _, hits = scan(examples / "aav2.fasta", "AAV-trs")
    assert len(hits) == 2
    assert {hit.matched for hit in hits} == {"GTTGGCC"}


def test_aav2_promoter_tata_boxes_sit_at_p5_p19_and_p40(examples):
    # The three promoters are named for their positions in map units, which on
    # a 4,679 base genome puts them near bases 230, 890 and 1,870.
    record, hits = scan(examples / "aav2.fasta", "TATA", mismatches=1)
    found = positions(hits, strand="+")
    assert 255 in found
    assert 843 in found
    assert 1823 in found
    for start, expected in ((255, 5), (843, 19), (1823, 40)):
        units = round(100 * start / len(record))
        assert abs(units - expected) <= 3


def test_aav2_polyadenylation_signal_is_shared_by_the_transcripts(examples):
    _, hits = scan(examples / "aav2.fasta", "polyA")
    # The single poly(A) site used by all AAV2 transcripts lies near map unit
    # 96; the other hit is in the body of the genome.
    assert 4424 in positions(hits, strand="+")


def test_aav2_has_no_kozak_consensus_at_full_strictness(examples):
    _, hits = scan(examples / "aav2.fasta", "Kozak")
    assert hits == []


# --- AAV2 VP1 capsid protein, UniProt P03135 ------------------------------


def test_vp1_carries_the_phospholipase_motif_in_its_unique_region(examples):
    record, hits = scan(examples / "aav2_vp1.fasta", "PLA2")
    (hit,) = hits
    assert hit.matched == "HDKAY"
    assert hit.start == 75
    # VP1u is the first ~137 residues, which VP2 and VP3 do not carry.
    assert hit.end < 137
    assert len(record) == 735


def test_vp1_hits_are_unstranded(examples):
    _, hits = scan(examples / "aav2_vp1.fasta", "PLA2")
    assert all(hit.strand == "." for hit in hits)


def test_vp1_n_glycosylation_sequons(examples):
    _, hits = scan(examples / "aav2_vp1.fasta", "N-glycosylation")
    assert positions(hits) == [223, 335, 382, 408, 496, 705]


def test_asking_for_a_reverse_strand_of_the_capsid_is_refused(examples):
    from motif_scan.patterns import PatternError

    with pytest.raises(PatternError, match="no reverse strand"):
        scan(examples / "aav2_vp1.fasta", "PLA2", strand="reverse")


# --- the catalogue against the examples -----------------------------------


def test_every_nucleotide_motif_can_be_applied_to_the_genome(examples):
    (record,) = read_records(examples / "aav2.fasta")
    for motif in library.select(alphabet="nucleotide"):
        find(record.sequence, compile_motif(motif), record=record.identifier)


def test_every_protein_motif_can_be_applied_to_the_capsid(examples):
    (record,) = read_records(examples / "aav2_vp1.fasta")
    for motif in library.select(alphabet="protein"):
        find(record.sequence, compile_motif(motif), record=record.identifier)


def test_the_examples_are_read_in_the_alphabet_they_are_written_in(examples):
    from motif_scan.alphabet import detect

    for name, expected in (
        ("aav2.fasta", "nucleotide"),
        ("hbb.fasta", "nucleotide"),
        ("aav2_vp1.fasta", "protein"),
    ):
        (record,) = read_records(examples / name)
        assert detect(record.sequence) == expected, name
