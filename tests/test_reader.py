"""FASTA input."""

import bz2
import gzip
import io
import lzma

import pytest

from motif_scan.alphabet import NUCLEOTIDE, PROTEIN
from motif_scan.reader import (
    Record,
    SequenceError,
    check_alphabet,
    parse_fasta,
    read_records,
    residues_outside,
)

FASTA = ">one first record\nACGTACGT\nACGT\n>two\nTTTT\n"


def records(text):
    return list(parse_fasta(io.StringIO(text)))


def test_reads_identifier_description_and_sequence():
    first, second = records(FASTA)
    assert first.identifier == "one"
    assert first.description == "first record"
    assert first.sequence == "ACGTACGTACGT"
    assert second.identifier == "two"
    assert second.description == ""


def test_sequence_is_uppercased_and_whitespace_removed():
    (record,) = records(">x\nacg t\n a cgt\n")
    assert record.sequence == "ACGTACGT"


def test_length_and_header():
    (record,) = records(">x a b\nACGT\n")
    assert len(record) == 4
    assert record.header == "x a b"


def test_label_is_qualified_by_file_when_given():
    assert Record("x", "", "ACGT", "f.fa").label == "f.fa:x"
    assert Record("x", "", "ACGT").label == "x"


def test_blank_lines_are_ignored():
    (record,) = records(">x\n\nACGT\n\n")
    assert record.sequence == "ACGT"


def test_byte_order_mark_is_stripped():
    (record,) = records("\ufeff>x\nACGT\n")
    assert record.identifier == "x"


def test_header_without_identifier_is_refused():
    with pytest.raises(SequenceError, match="no identifier"):
        records(">\nACGT\n")


def test_sequence_before_a_header_is_refused():
    with pytest.raises(SequenceError, match="before any"):
        records("ACGT\n>x\nACGT\n")


def test_input_without_records_is_refused():
    with pytest.raises(SequenceError, match="no FASTA records"):
        records("\n\n")


def test_prose_is_refused_rather_than_scanned():
    with pytest.raises(SequenceError, match="not sequence data"):
        records(">x\nthe quick brown fox, which is not sequence\n")


def test_the_offending_characters_are_named():
    with pytest.raises(SequenceError, match=r"'1'"):
        records(">x\nACGT123\n")


def test_many_offending_characters_are_abbreviated():
    with pytest.raises(SequenceError, match=r"\.\.\."):
        records(">x\n1234567\n")


def test_both_alphabets_are_accepted_at_read_time():
    # Which alphabet applies is settled after reading, so a peptide and a
    # nucleotide sequence both parse here.
    (peptide,) = records(">x\nMKQLEDKVEELL\n")
    assert peptide.sequence == "MKQLEDKVEELL"


def test_check_alphabet_accepts_matching_residues():
    check_alphabet(Record("x", "", "ACGTN"), NUCLEOTIDE)
    check_alphabet(Record("x", "", "MKQLE"), PROTEIN)


def test_check_alphabet_names_residues_from_the_other_alphabet():
    with pytest.raises(SequenceError, match="try --type protein"):
        check_alphabet(Record("x", "", "MKQLEDK"), NUCLEOTIDE)


def test_check_alphabet_suggests_n_for_an_x_in_nucleotide():
    with pytest.raises(SequenceError, match="written N, not X"):
        check_alphabet(Record("x", "", "ACGTXACGT"), NUCLEOTIDE)


def test_residues_outside_lists_them_sorted():
    assert residues_outside("ACGTEF", NUCLEOTIDE) == ["E", "F"]


def test_composition_counts_each_code():
    counts = Record("x", "", "AACGN").composition()
    assert counts["A"] == 2
    assert counts["N"] == 1


def test_ambiguous_counts_codes_standing_for_several_residues():
    assert Record("x", "", "ACGTNRY").ambiguous(NUCLEOTIDE) == 3
    assert Record("x", "", "MKQLX").ambiguous(PROTEIN) == 1


def test_reads_from_a_path(tmp_path):
    path = tmp_path / "x.fasta"
    path.write_text(FASTA)
    assert len(read_records(path)) == 2


def test_limit_stops_early(tmp_path):
    path = tmp_path / "x.fasta"
    path.write_text(FASTA)
    assert len(read_records(path, limit=1)) == 1


def test_label_is_recorded(tmp_path):
    path = tmp_path / "x.fasta"
    path.write_text(FASTA)
    assert read_records(path, label="x.fasta")[0].source == "x.fasta"


@pytest.mark.parametrize(
    "compress", [gzip.compress, bz2.compress, lzma.compress]
)
def test_compression_is_recognised_from_the_magic_bytes(tmp_path, compress):
    path = tmp_path / "x.bin"  # deliberately not named .gz, .bz2 or .xz
    path.write_bytes(compress(FASTA.encode()))
    assert len(read_records(path)) == 2


def test_a_missing_file_is_reported(tmp_path):
    with pytest.raises(FileNotFoundError, match="no such file"):
        read_records(tmp_path / "absent.fasta")


def test_a_directory_is_reported(tmp_path):
    with pytest.raises(IsADirectoryError, match="not a file"):
        read_records(tmp_path)


def test_a_truncated_gzip_is_reported_as_an_input_error(tmp_path):
    path = tmp_path / "x.fasta.gz"
    path.write_bytes(gzip.compress(FASTA.encode())[:-4])
    with pytest.raises(SequenceError, match="truncated|cannot be read"):
        read_records(path)


def test_utf16_is_reported_as_not_being_fasta(tmp_path):
    path = tmp_path / "x.fasta"
    path.write_bytes(FASTA.encode("utf-16"))
    with pytest.raises(SequenceError, match="not UTF-8"):
        read_records(path)
