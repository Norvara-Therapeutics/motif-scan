"""Text maps."""

from motif_scan import mapview
from motif_scan.search import Hit


def hit(start, end, motif="m", strand="+", matched=None):
    return Hit("r", motif, start, end, strand, matched or "A" * (end - start + 1))


def test_legend_gives_each_motif_a_symbol():
    marks = mapview.legend(["a", "b"])
    assert marks.of("a") == "a"
    assert marks.of("b") == "b"
    assert marks.of("absent") == mapview.SPARE


def test_legend_lines_name_every_motif():
    lines = mapview.legend(["TATA"]).lines()
    assert lines == ["  a  TATA"]


def test_legend_falls_back_when_motifs_outnumber_symbols():
    marks = mapview.legend([str(index) for index in range(len(mapview.SYMBOLS) + 2)])
    assert marks.symbols[-1][1] == mapview.SPARE


def test_column_of_scales_positions():
    assert mapview.column_of(1, 100, 10) == 0
    assert mapview.column_of(100, 100, 10) == 9
    assert mapview.column_of(51, 100, 10) == 5


def test_column_of_an_empty_record():
    assert mapview.column_of(1, 0, 10) == 0


def test_track_marks_a_single_hit_and_a_crowd():
    line = mapview.track([hit(1, 2), hit(2, 3)], 100, 10)
    assert line[0] == mapview.CROWDED
    line = mapview.track([hit(1, 2)], 100, 10)
    assert line[0] == mapview.MARK


def test_track_of_nothing_is_blank():
    assert mapview.track([], 100, 10).strip() == ""


def test_ruler_marks_the_topology():
    (line, labels) = mapview.ruler(100, 20, circular=True)
    assert line.startswith("o") and line.endswith("o")
    (line, labels) = mapview.ruler(100, 20, circular=False)
    assert line.startswith("[") and line.endswith("[")
    assert labels.startswith("1") and labels.endswith("100")


def test_ruler_never_exceeds_the_width_asked_for():
    for width in range(1, 12):
        (line, _) = mapview.ruler(100, width, False)
        assert len(line) <= width


def test_tracks_keep_a_row_for_a_motif_without_hits():
    lines = mapview.tracks([hit(1, 4, "a")], ["a", "b"], 100, width=20)
    assert any(line.strip().startswith("a") for line in lines)
    assert any(line.strip().startswith("b") for line in lines)


def test_tracks_end_with_a_ruler():
    lines = mapview.tracks([], ["a"], 100, width=20)
    assert lines[-1].strip().endswith("100")


def test_tracks_of_an_empty_record_are_empty():
    assert mapview.tracks([], ["a"], 0, width=20) == []


def test_blocks_print_the_residues_with_the_gutter():
    lines = mapview.blocks("ACGT" * 20, [hit(1, 4)], width=40)
    assert lines[0].strip().startswith("1 ACGT")
    assert "41 " in lines[2]


def test_blocks_underline_a_hit_with_its_symbol():
    lines = mapview.blocks("ACGTACGT", [hit(1, 4)], width=8)
    assert lines[1].strip().startswith("aaa")


def test_blocks_mark_direction():
    forward = mapview.blocks("ACGTACGT", [hit(1, 4, strand="+")], width=8)
    assert forward[1].rstrip().endswith(">")
    reverse = mapview.blocks("ACGTACGT", [hit(1, 4, strand="-")], width=8)
    assert reverse[1].strip().startswith("<")


def test_blocks_give_each_motif_its_own_line():
    marks = mapview.legend(["a", "b"])
    lines = mapview.blocks(
        "ACGTACGT", [hit(1, 4, "a"), hit(3, 6, "b")], width=8, names=marks
    )
    assert len(lines) == 3


def test_blocks_split_a_hit_that_wraps():
    # end < start means the hit spans the join of a circle.
    wrapped = Hit("r", "m", 7, 2, "+", "ACGT")
    lines = mapview.blocks("ACGTACGT", [wrapped], width=8)
    annotation = lines[1]
    assert annotation.rstrip().endswith("a") or annotation.rstrip().endswith(">")
    assert "a" in annotation


def test_blocks_of_an_empty_sequence():
    assert mapview.blocks("", [hit(1, 2)]) == []
