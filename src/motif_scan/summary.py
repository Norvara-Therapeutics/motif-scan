"""Counting hits, expecting hits, and placing them along the sequence.

Three questions are answered here, and the second is the one that needs care.

*How many.* Hits per motif, split by strand, by whether every position matched,
and by whether any position could not be called.

*How many would a sequence like this have anyway.* A motif six bases long turns
up about once every four kilobases by chance, and in a GC-rich payload a
GC-rich motif turns up far more often than that. The expected count here is
computed from the record's own residue composition, position by position over
the pattern, with the mismatch budget taken into account exactly rather than
approximated. It assumes residues are independent of one another, which real
sequence is not: the expectation is a yardstick for reading a count, not a
test, and nothing here reports a p-value.

For the yardstick to be read beside the count, both have to be over the same
windows. A window carrying an ambiguity code cannot yield a hit, so it is left
out of the expectation too - otherwise a record of one adenine followed by nine
Ns reports one hit against ten expected and reads as a tenfold depletion of a
motif that matched every position it could be assessed at.

*Where.* Counts per bin along the sequence, which is what the distribution
chart draws.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Mapping, Sequence

from motif_scan.alphabet import Alphabet
from motif_scan.patterns import Compiled, complement_positions
from motif_scan.search import Hit

#: Default number of bins the sequence is divided into for the distribution.
DEFAULT_BINS = 40


@dataclass(frozen=True)
class MotifCount:
    """What one motif found in one record.

    Attributes:
        motif: name of the motif.
        hits: how many occurrences were reported.
        forward: hits on the forward strand.
        reverse: hits on the reverse strand.
        palindromic: hits of a motif that reads the same on both strands, which
            are reported once rather than twice.
        unstranded: hits in a sequence with only one strand, which is protein.
        exact: hits where every position matched definitely - no mismatch, and
            no position resting on an ambiguity code in the sequence.
        with_mismatches: hits where at least one position did not match.
        listed_uncertain: hits carrying at least one position that could not be
            called, which is zero unless uncertain windows were asked for. A
            hit may be both inexact and uncertain, so these three do not
            partition the hits.
        uncertain: windows that could not be called either way because the
            sequence carries an ambiguity code in them. Not included in
            ``hits`` unless they were explicitly asked for.
        assessed: whether uncertain windows could be looked for at all. False
            for a variable-length pattern, where ``uncertain`` is 0 because
            nothing was assessed, not because nothing was found.
        expected: how many hits a sequence of this length and composition would
            carry by chance, or None when the pattern has no fixed length and
            so no expectation can be computed.
        residues: length of the record, for density.
    """

    motif: str
    hits: int
    forward: int
    reverse: int
    palindromic: int
    unstranded: int
    exact: int
    with_mismatches: int
    listed_uncertain: int
    uncertain: int
    assessed: bool
    expected: float | None
    residues: int

    @property
    def per_kb(self) -> float:
        """Hits per thousand residues."""
        if not self.residues:
            return 0.0
        return self.hits * 1000 / self.residues

    @property
    def enrichment(self) -> float | None:
        """Observed hits divided by expected hits.

        None when there is no expectation to divide by, and None rather than
        infinity when the expectation rounds to zero: a ratio against zero
        expectation says only that the motif is rare, which the counts already
        say more honestly.
        """
        if self.expected is None or self.expected <= 0:
            return None
        return self.hits / self.expected


@dataclass(frozen=True)
class Distribution:
    """Where the hits of each motif fall along one record.

    Attributes:
        bins: how many equal parts the record was divided into.
        residues: length of the record.
        counts: for each motif, the hits in each bin.
    """

    bins: int
    residues: int
    counts: tuple[tuple[str, tuple[int, ...]], ...]

    @property
    def width(self) -> float:
        """How many residues one bin covers."""
        return self.residues / self.bins if self.bins else 0.0

    def span(self, index: int) -> tuple[int, int]:
        """The 1-based residue range one bin covers, inclusive."""
        start = int(index * self.width) + 1
        end = int((index + 1) * self.width) if index + 1 < self.bins else self.residues
        return start, max(start, end)

    @property
    def busiest(self) -> int:
        """The largest count in any bin of any motif, for scaling a chart."""
        return max((max(counts) for _, counts in self.counts if counts), default=0)


@dataclass(frozen=True)
class RecordScan:
    """One record, scanned for every motif.

    Attributes:
        identifier: the record's identifier.
        description: the rest of its header line.
        source: the file it came from.
        residues: how long it is.
        alphabet: which alphabet it was read in.
        circular: whether it was treated as a closed circle.
        ambiguous: how many residues stand for more than one residue.
        hits: every hit reported, sorted by position.
        counts: one entry per motif, in the order the motifs were given.
        distribution: hits per bin along the record.
    """

    identifier: str
    description: str
    source: str
    residues: int
    alphabet: str
    circular: bool
    ambiguous: int
    hits: tuple[Hit, ...]
    counts: tuple[MotifCount, ...]
    distribution: Distribution

    @property
    def label(self) -> str:
        """Identifier, qualified by file when several files were read."""
        return f"{self.source}:{self.identifier}" if self.source else self.identifier

    @property
    def topology(self) -> str:
        """``circular`` or ``linear``."""
        return "circular" if self.circular else "linear"

    def hits_for(self, motif: str) -> tuple[Hit, ...]:
        """Every hit of one motif, in position order."""
        return tuple(hit for hit in self.hits if hit.motif == motif)

    def count_for(self, motif: str) -> MotifCount | None:
        """The counts for one motif, or None if it was not applied."""
        for count in self.counts:
            if count.motif == motif:
                return count
        return None


@dataclass(frozen=True)
class Scan:
    """A whole run: every motif against every record.

    Attributes:
        name: what the scan is of, for titling a report.
        records: one entry per record scanned.
        motifs: the motifs applied, compiled.
        parameters: the settings the scan ran under, for reporting.
        mismatches: the mismatch budget every motif was searched under.
        overlapping: whether overlapping hits were all reported. When false,
            no expectation is reported, because the two counts would no longer
            be answering the same question.
        include_uncertain: whether windows that could not be called were listed
            as hits. When they were, they are inside ``hits`` as well as being
            counted separately, and the expectation is taken over every window
            rather than only the assessable ones.
    """

    name: str
    records: tuple[RecordScan, ...]
    motifs: tuple[Compiled, ...]
    parameters: tuple[tuple[str, str], ...] = ()
    mismatches: int = 0
    overlapping: bool = True
    include_uncertain: bool = False

    @property
    def residues(self) -> int:
        """Total residues scanned."""
        return sum(record.residues for record in self.records)

    @property
    def hits(self) -> int:
        """Total hits reported."""
        return sum(len(record.hits) for record in self.records)

    @property
    def uncertain(self) -> int:
        """Total windows that could not be called."""
        return sum(
            count.uncertain for record in self.records for count in record.counts
        )

    @property
    def alphabet(self) -> str:
        """The alphabet the records were read in, or ``mixed``."""
        names = {record.alphabet for record in self.records}
        if len(names) == 1:
            return names.pop()
        return "mixed"

    def totals(self) -> tuple[MotifCount, ...]:
        """Counts for each motif, summed over every record.

        Density is over the residues of every record together, and an
        expectation is reported only when it could be computed for every record
        the motif was applied to: a total that silently omitted one record
        would read as a smaller expectation rather than an incomplete one.
        """
        totals: list[MotifCount] = []
        for compiled in self.motifs:
            parts = [
                count
                for record in self.records
                if (count := record.count_for(compiled.name)) is not None
            ]
            if not parts:
                continue
            expectations = [part.expected for part in parts]
            totals.append(
                MotifCount(
                    motif=compiled.name,
                    hits=sum(part.hits for part in parts),
                    forward=sum(part.forward for part in parts),
                    reverse=sum(part.reverse for part in parts),
                    palindromic=sum(part.palindromic for part in parts),
                    unstranded=sum(part.unstranded for part in parts),
                    exact=sum(part.exact for part in parts),
                    with_mismatches=sum(part.with_mismatches for part in parts),
                    listed_uncertain=sum(part.listed_uncertain for part in parts),
                    uncertain=sum(part.uncertain for part in parts),
                    assessed=all(part.assessed for part in parts),
                    expected=(
                        None
                        if any(value is None for value in expectations)
                        else sum(value for value in expectations if value is not None)
                    ),
                    residues=sum(part.residues for part in parts),
                )
            )
        return tuple(totals)


def frequencies(
    composition: Mapping[str, int], compiled: Compiled
) -> dict[str, float] | None:
    """Residue frequencies of a record, over the residues that are known.

    Ambiguity codes are left out of the denominator rather than spread over the
    residues they cover. Spreading them would be a guess about the very
    positions the scan refuses to guess about, and a record that is mostly Ns
    would then be given a composition it was never shown to have.

    Returns None when the record holds no unambiguous residue at all, which
    leaves nothing to compute an expectation from.
    """
    definite = compiled.alphabet.definite
    counted = {
        residue: count
        for residue, count in composition.items()
        if residue in definite
    }
    total = sum(counted.values())
    if not total:
        return None
    # A code standing for one residue is counted as that residue, so that U in
    # an RNA record contributes to T and selenocysteine to itself.
    spread: Counter[str] = Counter()
    for residue, count in counted.items():
        (single,) = compiled.alphabet.expand(residue)
        spread[single] += count
    return {residue: count / total for residue, count in spread.items()}


def assessable_windows(
    sequence: str,
    span: int,
    alphabet: Alphabet,
    circular: bool = False,
) -> int:
    """How many windows of a given width carry no ambiguity code.

    These are the windows a hit can come from, and so the windows an
    expectation should be computed over. A window holding an ambiguity code is
    reported as uncertain rather than as a hit, and counting it on the expected
    side but not the observed side is what turns a run of Ns into an apparent
    depletion.
    """
    length = len(sequence)
    if span <= 0 or length < span:
        return 0
    definite = alphabet.definite
    text = sequence + sequence[: span - 1] if circular and span > 1 else sequence
    total = length if circular else length - span + 1
    found = 0
    run = 0
    for index, residue in enumerate(text):
        run = run + 1 if residue in definite else 0
        if run >= span and index - span + 1 < total:
            found += 1
    return found


def _match_probability(
    positions: Sequence[frozenset[str]],
    frequency: Mapping[str, float],
    mismatches: int,
) -> float:
    """The chance that a window of independent residues matches a pattern.

    A dynamic programme over the positions, carrying the distribution of the
    number of mismatches so far. Enumerating the words a degenerate pattern
    with a mismatch budget accepts would be hopeless - a ten-position pattern
    with two mismatches allowed covers millions - while this is linear in the
    length of the pattern times the budget, and exact under the model.
    """
    distribution = [1.0] + [0.0] * mismatches
    for allowed in positions:
        hit = sum(frequency.get(residue, 0.0) for residue in allowed)
        miss = 1.0 - hit
        updated = [0.0] * (mismatches + 1)
        for wrong, chance in enumerate(distribution):
            if not chance:
                continue
            updated[wrong] += chance * hit
            if wrong + 1 <= mismatches:
                updated[wrong + 1] += chance * miss
        distribution = updated
    return sum(distribution)


def expected_count(
    compiled: Compiled,
    composition: Mapping[str, int],
    length: int,
    *,
    mismatches: int = 0,
    strand: str = "both",
    circular: bool = False,
    windows: int | None = None,
) -> float | None:
    """How many hits a record of this length and composition would carry by chance.

    Args:
        windows: how many windows to count the expectation over. Defaults to
            every window of the record, which is right only when every window
            could have yielded a hit; pass the count from
            :func:`assessable_windows` for a record carrying ambiguity codes.

    Returns None when no expectation can be computed: a variable-length
    pattern has no positions to compute one over, and a record with no
    unambiguous residues gives nothing to compute it from.

    The model is the simplest one that can be stated honestly: residues drawn
    independently from the record's own composition. It therefore ignores
    everything real sequence does - dinucleotide bias, repeats, the CpG
    depletion of a mammalian genome - and a motif can sit well above or below
    this expectation for reasons that have nothing to do with selection.
    """
    if not compiled.positions:
        return None
    frequency = frequencies(composition, compiled)
    if frequency is None:
        return None

    span = len(compiled.positions)
    if length < span:
        return 0.0
    if windows is None:
        if compiled.anchored:
            windows = 1
        elif circular:
            windows = length
        else:
            windows = length - span + 1
    elif compiled.anchored:
        windows = min(windows, 1)

    passes = [compiled.positions]
    if compiled.alphabet.has_strands and strand == "both" and not compiled.palindromic:
        passes.append(complement_positions(compiled.positions, compiled.alphabet))
    elif compiled.alphabet.has_strands and strand == "reverse":
        passes = [complement_positions(compiled.positions, compiled.alphabet)]

    counted = windows
    return sum(
        _match_probability(positions, frequency, mismatches) * counted
        for positions in passes
    )


def tally(
    compiled: Compiled,
    hits: Sequence[Hit],
    length: int,
    composition: Mapping[str, int],
    *,
    mismatches: int = 0,
    strand: str = "both",
    circular: bool = False,
    uncertain: int = 0,
    overlapping: bool = True,
    windows: int | None = None,
) -> MotifCount:
    """Count one motif's hits in one record.

    No expectation is reported when ``overlapping`` is false. The expectation
    counts every window the motif matches, while suppressing overlaps reports
    only the leftmost of each run, so the ratio of the two would move with a
    reporting choice rather than with the sequence: a run of five alanines
    holds four AA windows and is reported as two hits, which would read as a
    twofold depletion of a motif that is in fact exactly as common as chance
    predicts. Rather than divide one question by another, the figure is
    withheld and the reason given.
    """
    strands = Counter(hit.strand for hit in hits)
    return MotifCount(
        motif=compiled.name,
        hits=len(hits),
        forward=strands["+"],
        reverse=strands["-"],
        palindromic=strands["="],
        unstranded=strands["."],
        exact=sum(1 for hit in hits if hit.exact and hit.certain),
        with_mismatches=sum(1 for hit in hits if not hit.exact),
        listed_uncertain=sum(1 for hit in hits if not hit.certain),
        uncertain=uncertain,
        assessed=compiled.supports_mismatches,
        expected=(
            expected_count(
                compiled,
                composition,
                length,
                mismatches=mismatches,
                strand=strand,
                circular=circular,
                windows=windows,
            )
            if overlapping
            else None
        ),
        residues=length,
    )


def distribute(
    hits: Sequence[Hit],
    motifs: Sequence[str],
    length: int,
    bins: int = DEFAULT_BINS,
) -> Distribution:
    """Count each motif's hits in equal bins along a record.

    A hit is placed by its start, so a motif spanning the join of a circle
    counts in the last bin rather than the first. Bins are never wider than one
    residue: a record shorter than the number of bins asked for gets one bin
    per residue, so the chart never implies a resolution the sequence has not
    got.
    """
    width = max(1, bins)
    if length <= 0:
        return Distribution(0, 0, tuple((motif, ()) for motif in motifs))
    width = min(width, length)
    counts: dict[str, list[int]] = {motif: [0] * width for motif in motifs}
    for hit in hits:
        if hit.motif not in counts:
            continue
        index = min(width - 1, (hit.start - 1) * width // length)
        counts[hit.motif][index] += 1
    return Distribution(
        bins=width,
        residues=length,
        counts=tuple((motif, tuple(counts[motif])) for motif in motifs),
    )
