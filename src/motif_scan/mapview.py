"""Text maps of where the motifs fall.

Two views, for two questions.

:func:`tracks` draws one row per motif, scaled across the whole record, which
answers *where along the sequence* - a cluster of sites in a promoter looks
like a cluster, and a motif that peppers the whole record looks like that
instead. Where several hits land on one column of a narrow terminal they are
marked as crowded rather than overprinted, because a map that silently loses a
hit is worse than one that admits it ran out of room.

:func:`blocks` prints the residues themselves, wrapped, with the span of each
hit underlined beneath the line it falls on. That is only readable for a short
record, so the caller decides when to ask for it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from motif_scan.search import Hit
from motif_scan.text import thousands

#: Character marking a single hit in a column of a track.
MARK = "|"

#: Character marking a column holding more than one hit.
CROWDED = "*"

#: Symbols given out to motifs, in order, for the block view. Letters first so
#: that a legend reads as a list rather than as punctuation.
SYMBOLS = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"

#: Used when there are more motifs than symbols.
SPARE = "+"


@dataclass(frozen=True)
class Legend:
    """Which symbol stands for which motif in the block view."""

    symbols: tuple[tuple[str, str], ...]

    def of(self, motif: str) -> str:
        """The symbol standing for one motif."""
        for name, symbol in self.symbols:
            if name == motif:
                return symbol
        return SPARE

    def lines(self) -> list[str]:
        """The legend, one motif per line."""
        return [f"  {symbol}  {name}" for name, symbol in self.symbols]


def legend(motifs: Sequence[str]) -> Legend:
    """Give each motif a symbol, in the order the motifs were applied."""
    return Legend(
        tuple(
            (motif, SYMBOLS[index] if index < len(SYMBOLS) else SPARE)
            for index, motif in enumerate(motifs)
        )
    )


def ruler(length: int, width: int, circular: bool) -> list[str]:
    """Draw the record as a scaled line with its end positions.

    The line is never drawn wider than the width asked for, which at a width of
    one or two leaves room for the end caps and nothing between them.
    """
    caps, fill = ("o", "=") if circular else ("[", "-")
    if width < 2:
        line = caps * max(0, width)
    else:
        line = caps + fill * (width - 2) + caps
    labels = f"1{' ' * max(1, width - len(thousands(length)) - 1)}{thousands(length)}"
    return [line, labels]


def column_of(position: int, length: int, width: int) -> int:
    """Which column of a scaled track a position falls in."""
    if length <= 0 or width <= 0:
        return 0
    return min(width - 1, (position - 1) * width // length)


def track(hits: Sequence[Hit], length: int, width: int) -> str:
    """One row of marks, scaled across the record."""
    if width <= 0 or length <= 0:
        return ""
    counts = [0] * width
    for hit in hits:
        counts[column_of(hit.start, length, width)] += 1
    return "".join(
        " " if not count else (MARK if count == 1 else CROWDED) for count in counts
    )


def tracks(
    hits: Sequence[Hit],
    motifs: Sequence[str],
    length: int,
    width: int = 48,
    circular: bool = False,
    label_width: int = 16,
) -> list[str]:
    """Draw one scaled row per motif, above a ruler of the whole record.

    Motifs with no hits are drawn too, as an empty row. Leaving them out would
    make an absent element look like a motif nobody searched for, and the
    difference between those two matters more than the space the row costs.
    """
    if length <= 0 or width <= 0:
        return []
    lines = []
    for motif in motifs:
        found = [hit for hit in hits if hit.motif == motif]
        name = motif if len(motif) <= label_width else motif[: label_width - 1] + "\u2026"
        lines.append(
            f"  {name:<{label_width}} {track(found, length, width)}"
            f"  {thousands(len(found))}"
        )
    pad = " " * (label_width + 3)
    lines.extend(f"{pad}{line}" for line in ruler(length, width, circular))
    return lines


def _runs(hit: Hit, length: int) -> list[tuple[int, int]]:
    """The 0-based half-open spans a hit covers, splitting one that wraps."""
    first = hit.start - 1
    if not hit.spans_origin:
        return [(first, first + hit.length)]
    return [(first, length), (0, hit.end)]


def blocks(
    sequence: str,
    hits: Sequence[Hit],
    width: int = 60,
    names: Legend | None = None,
    gutter: int = 9,
) -> list[str]:
    """Print the residues, wrapped, with each hit underlined beneath its line.

    Each motif with a hit on a line gets its own annotation line, so overlapping
    motifs stay legible rather than fighting for the same column. A hit on the
    reverse strand opens with ``<`` and one on the forward strand closes with
    ``>``, which is the only place the direction appears in this view.
    """
    length = len(sequence)
    if not length or width <= 0:
        return []
    marks = names or legend(sorted({hit.motif for hit in hits}))

    spans: list[tuple[int, int, Hit]] = [
        (start, end, hit) for hit in hits for start, end in _runs(hit, length)
    ]
    lines: list[str] = []
    for start in range(0, length, width):
        stop = min(start + width, length)
        lines.append(f"{thousands(start + 1):>{gutter}} {sequence[start:stop]}")
        here = [entry for entry in spans if entry[0] < stop and entry[1] > start]
        for motif, _ in marks.symbols:
            drawn = [entry for entry in here if entry[2].motif == motif]
            if not drawn:
                continue
            row = [" "] * (stop - start)
            symbol = marks.of(motif)
            for first, last, hit in drawn:
                for index in range(max(first, start), min(last, stop)):
                    row[index - start] = symbol
                if hit.strand == "-" and start <= first < stop:
                    row[first - start] = "<"
                if hit.strand != "-" and start < last <= stop:
                    row[last - 1 - start] = ">"
            lines.append(f"{' ' * gutter} {''.join(row).rstrip()}")
    return lines
