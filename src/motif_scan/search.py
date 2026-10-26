"""Finding motifs in a sequence.

Coordinates are 1-based and inclusive, and always given on the forward strand
whichever strand the motif was matched on. The matched text, by contrast, is
given as the motif reads it: a hit on the reverse strand shows the motif, not
the forward-strand text, because the point of reporting it is to see the
element that was found.

Four things here are easy to get wrong and are handled explicitly:

* A motif that is not its own reverse complement occurs on either strand, so
  both are searched. One that reads identically on both - ``CANNTG``, ``GGWCC``
  - is searched once and labelled ``=``, because a second pass would only
  rediscover the same positions. Protein has no second strand at all, and its
  hits are labelled ``.``.
* On a circular molecule a motif may span the join, so the search runs over the
  sequence with as much of its start appended as the motif could need - and no
  hit may be longer than the molecule itself, which a pattern of variable
  length will otherwise reach for, matching its way around the circle and past
  its own start.
* An ambiguity code in the *sequence* means the residue is not known. Such a
  position cannot be called a hit, but neither is it certainly not one, so it
  is reported separately as uncertain rather than silently dropped.
* Overlapping occurrences are all reported. A motif such as ``CG`` or
  ``CANNTG`` tiles along the sequence, and a scanner that resumed after each
  match would report every other one.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, Sequence

from motif_scan.alphabet import has_ambiguity, reverse_complement
from motif_scan.patterns import GROUP, Compiled, PatternError

#: Strand labels. ``=`` marks a motif that reads identically on both strands,
#: so that a palindromic motif is not reported twice; ``.`` marks a sequence
#: that has no second strand, which is every protein.
STRANDS: tuple[str, ...] = ("+", "-", "=", ".")

#: What ``--strand`` may ask for.
STRAND_CHOICES: tuple[str, ...] = ("both", "forward", "reverse")

#: How a position of a fixed-length pattern can meet a sequence residue.
_MATCH, _UNCERTAIN, _MISMATCH = 0, 1, 2


@dataclass(frozen=True)
class Hit:
    """One occurrence of a motif.

    Attributes:
        record: identifier of the record the motif was found in.
        motif: name of the motif.
        start: 1-based first residue of the hit on the forward strand.
        end: 1-based last residue, inclusive. On a circular molecule ``end``
            may be smaller than ``start``, which means the hit spans the join.
        strand: one of :data:`STRANDS`.
        matched: the text matched, read on the strand named.
        mismatches: positions that did not match. Always 0 for a
            variable-length pattern, which has no positions to count.
        uncertain: positions that could not be called either way, because the
            sequence carries an ambiguity code there.
    """

    record: str
    motif: str
    start: int
    end: int
    strand: str
    matched: str
    mismatches: int = 0
    uncertain: int = 0

    @property
    def length(self) -> int:
        """How many residues the hit covers."""
        return len(self.matched)

    @property
    def spans_origin(self) -> bool:
        """Whether the hit wraps past the end of a circular sequence."""
        return self.end < self.start

    @property
    def exact(self) -> bool:
        """Whether every position matched."""
        return self.mismatches == 0

    @property
    def certain(self) -> bool:
        """Whether no position depended on an ambiguity code in the sequence."""
        return self.uncertain == 0


def _lookup(
    compiled: Compiled,
) -> tuple[dict[str, int], ...]:
    """For each position, how every sequence code meets it.

    Built once per motif and read once per residue examined, which is what
    makes the mismatch-tolerant scan affordable in Python: the alternative is
    two set operations per residue per position.
    """
    table = []
    for allowed in compiled.positions:
        verdicts: dict[str, int] = {}
        for code, residues in compiled.alphabet.codes.items():
            if residues <= allowed:
                verdict = _MATCH
            elif residues & allowed:
                verdict = _UNCERTAIN
            else:
                verdict = _MISMATCH
            verdicts[code] = verdict
            verdicts[code.lower()] = verdict
        table.append(verdicts)
    return tuple(table)


def _window(sequence: str, span: int, circular: bool) -> str:
    """The text to search, extended so a circular hit can span the join.

    A motif of length ``span`` can begin at most ``span - 1`` residues before
    the join and still cross it, so that much of the start is appended.
    """
    if not circular or span < 2 or not sequence:
        return sequence
    overlap = min(span - 1, len(sequence))
    return sequence + sequence[:overlap]


def _scan_exact(compiled: Compiled, text: str) -> Iterator[tuple[int, str, int, int]]:
    """Find matches with no mismatches allowed, in one pass.

    Used for the patterns this package builds itself, which carry a zero-width
    lookahead so that overlapping occurrences are all found in a single scan.
    """
    for match in compiled.regex.finditer(text):
        matched = match.group(GROUP)
        if not matched:
            raise PatternError(
                f"{compiled.name}: pattern matched an empty string at position "
                f"{match.start() + 1}, so it cannot be counted"
            )
        yield match.start(), matched, 0, 0


def _scan_anchored(
    compiled: Compiled, text: str, stop: int, reach: int = 0
) -> Iterator[tuple[int, str, int, int]]:
    """Match at every start position in turn, optionally capping the reach.

    Slower than one pass of :func:`_scan_exact` by a constant factor, and
    needed in two cases that pass cannot serve.

    A user-written regular expression is compiled untouched, so it has no
    lookahead to make ``finditer`` return overlapping matches; testing each
    start position gives the same answer without rewriting somebody else's
    pattern.

    A pattern of variable length searched on a circle has to be stopped from
    matching more than once around the molecule. ``reach`` caps how far from
    the start position a match may extend, which the regular expression engine
    honours directly: without it, ``A(2,8)`` on a circle of four adenines
    matches eight residues and reports a four-residue span, the coordinates and
    the matched text then describing different things.
    """
    size = len(text)
    for offset in range(min(stop, size)):
        endpos = min(offset + reach, size) if reach else size
        match = compiled.regex.match(text, offset, endpos)
        if match is None:
            continue
        matched = compiled.matched(match)
        if not matched:
            raise PatternError(
                f"{compiled.name}: pattern matched an empty string at position "
                f"{offset + 1}, so it cannot be counted"
            )
        yield offset, matched, 0, 0


def _scan_tolerant(
    compiled: Compiled, text: str, mismatches: int
) -> Iterator[tuple[int, str, int, int]]:
    """Find matches position by position, allowing mismatches and uncertainty.

    A position where the sequence carries an ambiguity code that the motif does
    not fully cover is counted as uncertain, not as a mismatch. The two are
    different claims: a mismatch says the sequence is known and does not match,
    while an uncertain position says the sequence is not known. Spending the
    mismatch budget on missing information would turn every run of Ns into a
    hit for any motif short enough.
    """
    table = _lookup(compiled)
    span = len(table)
    budget = mismatches
    for offset in range(len(text) - span + 1):
        wrong = unknown = 0
        for index, verdicts in enumerate(table):
            verdict = verdicts.get(text[offset + index], _MISMATCH)
            if verdict == _MATCH:
                continue
            if verdict == _UNCERTAIN:
                unknown += 1
                continue
            wrong += 1
            if wrong > budget:
                break
        else:
            yield offset, text[offset: offset + span], wrong, unknown


def _strand_passes(
    compiled: Compiled, strand: str
) -> tuple[tuple[str, bool], ...]:
    """Which strands to search, and whether each is read reverse-complemented.

    Raises:
        PatternError: if the reverse strand is asked for in an alphabet that
            has none. A protein sequence has one strand, and silently ignoring
            ``--strand reverse`` would answer a different question from the one
            asked.
    """
    if not compiled.alphabet.has_strands:
        if strand == "reverse":
            raise PatternError(
                f"{compiled.name}: {compiled.alphabet.name} sequence has no "
                "reverse strand to search"
            )
        return ((".", False),)
    if strand == "forward":
        return (("+", False),)
    if strand == "reverse":
        return (("-", True),)
    if compiled.palindromic:
        return (("=", False),)
    return (("+", False), ("-", True))


def _occurrences(
    sequence: str,
    compiled: Compiled,
    mismatches: int,
    strand: str,
    circular: bool,
) -> Iterator[tuple[str, int, int, str, int, int]]:
    """Every window the motif matches, as forward-strand coordinates.

    Yields ``(strand, start, end, matched, mismatches, uncertain)`` with start
    and end 1-based and inclusive, end smaller than start for a hit spanning
    the join of a circle. Both :func:`find` and :func:`count_uncertain` are
    built on this, so the coordinate arithmetic - which is the part that is
    easy to get wrong - exists once.

    Raises:
        PatternError: if the search asked for cannot be performed.
    """
    if mismatches and not compiled.supports_mismatches:
        raise PatternError(
            f"{compiled.name}: mismatches cannot be allowed against a "
            f"{compiled.motif.syntax} pattern of variable length, because it has "
            "no fixed positions to count them against"
        )
    if mismatches < 0:
        raise PatternError(f"{compiled.name}: mismatches cannot be negative")

    length = len(sequence)
    if not length or length < compiled.minimum:
        return
    if circular and compiled.anchored:
        raise PatternError(
            f"{compiled.name}: a pattern anchored to a terminus cannot be "
            "searched on a circular sequence, which has no termini"
        )
    if circular and compiled.maximum is None:
        raise PatternError(
            f"{compiled.name}: a pattern of unbounded length cannot be searched "
            "across the join of a circular sequence; bound its repeats or scan "
            "the sequence as linear"
        )

    span = compiled.maximum or 0
    window = _window(sequence, span, circular)
    tolerant = bool(mismatches) or (
        compiled.supports_mismatches and has_ambiguity(sequence, compiled.alphabet)
    )

    for label, reversed_strand in _strand_passes(compiled, strand):
        text = reverse_complement(window) if reversed_strand else window
        if tolerant:
            matches = _scan_tolerant(compiled, text, mismatches)
        elif not compiled.wrapped:
            # Somebody else's regular expression: tested at each start position
            # rather than rewritten to find overlaps.
            matches = _scan_anchored(compiled, text, len(text))
        elif circular and compiled.fixed_length is None:
            # Variable length on a circle: capped so no hit laps the molecule.
            matches = _scan_anchored(compiled, text, length, min(span, length))
        else:
            matches = _scan_exact(compiled, text)
        for offset, matched, wrong, unknown in matches:
            size = len(matched)
            if size > length:
                # A hit cannot cover more residues than the molecule holds; the
                # scan above caps this, and this is the guard that it did.
                raise PatternError(
                    f"{compiled.name}: matched {size} residues in a sequence of "
                    f"{length}, which cannot be reported as one occurrence"
                )
            if compiled.anchor_start and offset != 0:
                continue
            if compiled.anchor_end and offset + size != length:
                continue
            if reversed_strand:
                # A hit found in the reverse complement of the window covers
                # the mirror-image span of the window itself.
                first = len(text) - offset - size
            else:
                first = offset
            if first >= length:
                # Already reported from the unwrapped part of the window.
                continue
            last = first + size
            if last > length:
                if not circular:
                    continue
                last -= length
            yield label, first + 1, last, matched, wrong, unknown


def find(
    sequence: str,
    compiled: Compiled,
    *,
    mismatches: int = 0,
    strand: str = "both",
    circular: bool = False,
    overlapping: bool = True,
    record: str = "",
    include_uncertain: bool = False,
) -> list[Hit]:
    """Every occurrence of one motif in one sequence.

    Args:
        sequence: the residues to search, in the motif's alphabet.
        compiled: the motif, already compiled.
        mismatches: how many positions may fail to match. Only a fixed-length
            pattern can allow any.
        strand: one of :data:`STRAND_CHOICES`.
        circular: treat the sequence as a closed circle, so a motif may span
            the join between its end and its start.
        overlapping: report every occurrence, including those overlapping one
            another. When false, the leftmost of a run of overlapping hits is
            kept, applied separately to each strand.
        record: identifier to record on each hit.
        include_uncertain: also return the windows that could not be called
            because the sequence carries an ambiguity code in them. They are
            left out by default: a run of Ns is compatible with any short motif
            and would otherwise fill the output with positions that are not
            findings. They are counted either way, by
            :func:`count_uncertain`.

    Returns:
        Hits sorted by forward-strand start, then strand, then motif.

    Raises:
        PatternError: if mismatches are asked of a pattern that cannot support
            them, if the reverse strand is asked for in protein, or if a
            circular search is asked for with a pattern that is anchored or of
            unbounded length.
    """
    found = [
        Hit(
            record=record,
            motif=compiled.name,
            start=start,
            end=end,
            strand=label,
            matched=matched,
            mismatches=wrong,
            uncertain=unknown,
        )
        for label, start, end, matched, wrong, unknown in _occurrences(
            sequence, compiled, mismatches, strand, circular
        )
        if include_uncertain or not unknown
    ]
    if not overlapping:
        found = _without_overlaps(found, len(sequence), circular)
    found.sort(key=lambda hit: (hit.start, hit.strand, hit.motif))
    return found


def count_uncertain(
    sequence: str,
    compiled: Compiled,
    *,
    mismatches: int = 0,
    strand: str = "both",
    circular: bool = False,
) -> int:
    """How many windows match but cannot be called, for want of information.

    A window counts when it would be a hit except that at least one of its
    residues is an ambiguity code in the *sequence*, so whether the motif is
    there depends on information the file does not carry. Counting rather than
    listing keeps the cost of a sequence with long runs of Ns proportional to
    its length, not to the number of motifs times the length of the runs.
    """
    if not compiled.supports_mismatches:
        # Uncertainty is assessed position by position, which a pattern of
        # variable length does not have. Reporting zero would be a claim; the
        # caller is told the assessment did not happen instead.
        return 0
    if not has_ambiguity(sequence, compiled.alphabet):
        return 0
    return sum(
        1
        for *_, unknown in _occurrences(
            sequence, compiled, mismatches, strand, circular
        )
        if unknown
    )


def covered_by(hit: Hit, length: int) -> set[int]:
    """The 0-based residues a hit covers, taking in the join of a circle."""
    first = hit.start - 1
    if not hit.spans_origin:
        return set(range(first, min(first + hit.length, length)))
    return set(range(first, length)) | set(range(0, hit.end))


def _without_overlaps(
    hits: Sequence[Hit], length: int, circular: bool = False
) -> list[Hit]:
    """Keep the leftmost of each run of overlapping hits, strand by strand.

    Strands are thinned separately: a motif on the forward strand and one on
    the reverse that cover the same residues are two findings, not one, and
    collapsing them would hide half of an inverted repeat.

    Overlap is decided over the residues each hit actually covers rather than
    by comparing a start against the furthest position reached. On a circular
    sequence the two differ: a hit that spans the join starts after everything
    kept so far and still shares residues with the first hit of the molecule,
    so position arithmetic alone lets it through.
    """
    kept: list[Hit] = []
    for label in STRANDS:
        occupied: set[int] = set()
        for hit in sorted(
            (hit for hit in hits if hit.strand == label),
            key=lambda hit: (hit.start, hit.length, hit.end),
        ):
            residues = covered_by(hit, length)
            if residues & occupied:
                continue
            kept.append(hit)
            occupied |= residues
    return kept
