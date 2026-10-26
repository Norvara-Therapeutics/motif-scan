"""Command-line interface."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import IO, Sequence

from motif_scan import __version__, library
from motif_scan.alphabet import (
    ALPHABETS,
    NUCLEOTIDE,
    NUCLEOTIDE_CODES,
    PROTEIN,
    PROTEIN_CODES,
    detect,
    get as get_alphabet,
)
from motif_scan.mapview import blocks, legend, ruler, track
from motif_scan.patterns import (
    Compiled,
    Motif,
    PatternError,
    compile_motif,
    named,
)
from motif_scan.reader import Record, SequenceError, check_alphabet, read_records
from motif_scan.report import render_html
from motif_scan.search import STRAND_CHOICES, Hit, count_uncertain, find
from motif_scan.summary import (
    DEFAULT_BINS,
    MotifCount,
    RecordScan,
    Scan,
    assessable_windows,
    distribute,
    tally,
)
from motif_scan.text import number, plural, ratio, thousands, truncate

STDOUT_FORMATS = ("table", "tsv", "json", "none")

_HIT_COLUMNS = (
    "record",
    "motif",
    "start",
    "end",
    "strand",
    "length",
    "match",
    "mismatches",
    "uncertain",
    "pattern",
    "syntax",
)

_COUNT_COLUMNS = (
    "record",
    "motif",
    "pattern",
    "syntax",
    "residues",
    "hits",
    "forward",
    "reverse",
    "palindromic",
    "unstranded",
    "exact",
    "with_mismatches",
    "listed_uncertain",
    "uncertain",
    "uncertainty_assessed",
    "per_kb",
    "expected",
    "observed_over_expected",
)

_MOTIF_COLUMNS = (
    "name",
    "pattern",
    "syntax",
    "alphabet",
    "category",
    "length",
    "description",
    "source",
)


def _non_negative(text: str) -> int:
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not an integer: {text!r}") from None
    if value < 0:
        raise argparse.ArgumentTypeError(f"must not be negative: {value}")
    return value


def _positive(text: str) -> int:
    value = _non_negative(text)
    if value < 1:
        raise argparse.ArgumentTypeError(f"must be 1 or greater: {value}")
    return value


def _motif_names(text: str) -> list[str]:
    names = [part.strip() for part in text.replace(" ", ",").split(",")]
    chosen = [name for name in names if name]
    if not chosen:
        raise argparse.ArgumentTypeError("no motif names given")
    return chosen


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser."""
    parser = argparse.ArgumentParser(
        prog="motif-scan",
        description=(
            "Find motifs in nucleotide or protein sequence: IUPAC degenerate "
            "patterns, PROSITE patterns, or regular expressions, with "
            f"mismatches allowed and both strands searched. "
            f"{len(library.all_motifs())} regulatory and protein elements are "
            "built in."
        ),
    )
    parser.add_argument(
        "--version", action="version", version=f"motif-scan {__version__}"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_scan_parser(subparsers)
    _add_motifs_parser(subparsers)
    _add_codes_parser(subparsers)
    return parser


def _add_scan_parser(subparsers: argparse._SubParsersAction) -> None:
    scan = subparsers.add_parser(
        "scan",
        help="find motifs in sequence",
        description=(
            "Scan one or more FASTA files for the motifs named, for patterns "
            "given on the command line, or for the whole built-in catalogue. "
            "Positions are 1-based and inclusive on the forward strand; the "
            "matched text is read on the strand reported."
        ),
    )
    scan.add_argument(
        "-i",
        "--input",
        required=True,
        action="append",
        metavar="PATH",
        help="FASTA file to scan, repeatable, or '-' for standard input",
    )
    chosen = scan.add_argument_group("motifs to look for")
    chosen.add_argument(
        "-m",
        "--motif",
        type=_motif_names,
        action="append",
        metavar="NAMES",
        help="catalogue motif to apply, comma-separated or repeated",
    )
    chosen.add_argument(
        "--all",
        action="store_true",
        help=(
            "apply every catalogue motif written in the alphabet of the input, "
            "subject to --category"
        ),
    )
    chosen.add_argument(
        "-p",
        "--pattern",
        action="append",
        metavar="[NAME=]IUPAC",
        help=(
            "IUPAC pattern to apply, as GCCRCCATGG or kozak=GCCRCCATGG; "
            "repeatable"
        ),
    )
    chosen.add_argument(
        "--prosite",
        action="append",
        metavar="[NAME=]PATTERN",
        help="PROSITE pattern to apply, as 'N-{P}-[ST]-{P}'; repeatable",
    )
    chosen.add_argument(
        "--regex",
        action="append",
        metavar="[NAME=]REGEX",
        help=(
            "regular expression to apply; repeatable. Mismatches cannot be "
            "allowed against one, and a pattern that can match an empty string "
            "is refused"
        ),
    )
    chosen.add_argument(
        "--motifs-file",
        metavar="PATH",
        help=(
            "tab-separated file of motifs, with at least the columns 'name' and "
            f"'pattern' (also accepted: {', '.join(library.COLUMNS[2:])})"
        ),
    )
    chosen.add_argument(
        "--category",
        metavar="NAME",
        help=(
            "with --all, apply only this category "
            f"({', '.join(library.categories())})"
        ),
    )

    how = scan.add_argument_group("how to search")
    how.add_argument(
        "--mismatches",
        type=_non_negative,
        default=0,
        metavar="N",
        help=(
            "allow up to N positions to differ (default: 0). Only a pattern of "
            "fixed length can allow any"
        ),
    )
    how.add_argument(
        "--strand",
        choices=STRAND_CHOICES,
        default="both",
        help=(
            "which strand of a nucleotide sequence to search (default: both). A "
            "motif reading identically on both strands is reported once"
        ),
    )
    how.add_argument(
        "-c",
        "--circular",
        action="store_true",
        help="treat each sequence as a closed circle, as a plasmid is",
    )
    how.add_argument(
        "--no-overlapping",
        action="store_true",
        help=(
            "report only the leftmost of each run of overlapping hits of the "
            "same motif on the same strand"
        ),
    )
    how.add_argument(
        "--type",
        choices=("auto", *ALPHABETS),
        default="auto",
        help=(
            "alphabet of the input (default: auto, judged from composition and "
            "always reported)"
        ),
    )
    how.add_argument(
        "--include-uncertain",
        action="store_true",
        help=(
            "list the windows that match except where the sequence carries an "
            "ambiguity code; they are counted either way"
        ),
    )

    out = scan.add_argument_group("output")
    out.add_argument(
        "-f",
        "--format",
        choices=STDOUT_FORMATS,
        default="table",
        help="what to write to standard output (default: table)",
    )
    out.add_argument(
        "-o",
        "--output",
        metavar="PATH",
        help="write every hit to this file as TSV",
    )
    out.add_argument(
        "--counts",
        metavar="PATH",
        help="write the per-motif counts to this file as TSV",
    )
    out.add_argument(
        "--report",
        metavar="PATH",
        help="write a self-contained HTML report to this file",
    )
    out.add_argument(
        "--limit",
        type=_non_negative,
        default=50,
        metavar="N",
        help=(
            "most hits to list per record in the table and the report; 0 lists "
            "every one. The TSV and JSON outputs are never truncated "
            "(default: 50)"
        ),
    )
    out.add_argument(
        "--sort",
        choices=("position", "motif"),
        default="position",
        help="order hits by position or by motif (default: position)",
    )
    out.add_argument(
        "--sequence",
        action="store_true",
        help="print the residues with every hit marked beneath them",
    )
    out.add_argument(
        "--no-map",
        action="store_true",
        help="leave out the scaled map of where the hits fall",
    )
    out.add_argument(
        "--bins",
        type=_positive,
        default=DEFAULT_BINS,
        metavar="N",
        help=(
            f"bins to divide each record into for the distribution chart in the "
            f"report (default: {DEFAULT_BINS})"
        ),
    )
    out.add_argument(
        "--map-width",
        type=_positive,
        default=48,
        metavar="N",
        help="width of the scaled map, in characters (default: 48)",
    )
    out.add_argument(
        "--line-width",
        type=_positive,
        default=60,
        metavar="N",
        help="residues per line where the sequence is printed (default: 60)",
    )
    out.add_argument(
        "--highlight-limit",
        type=_non_negative,
        default=20000,
        metavar="N",
        help=(
            "longest record whose residues are shown in full, in the report and "
            "with --sequence (default: 20000)"
        ),
    )


def _add_motifs_parser(subparsers: argparse._SubParsersAction) -> None:
    motifs = subparsers.add_parser(
        "motifs",
        help="list the built-in motif catalogue",
        description=(
            "List the motifs built in, with the pattern each is written as and "
            "the paper its consensus comes from."
        ),
    )
    motifs.add_argument(
        "--category",
        metavar="NAME",
        help=f"list only this category ({', '.join(library.categories())})",
    )
    motifs.add_argument(
        "--type",
        choices=tuple(ALPHABETS),
        help="list only motifs written in this alphabet",
    )
    motifs.add_argument(
        "--search",
        metavar="TEXT",
        help="list only motifs whose name, description or pattern contains this",
    )
    motifs.add_argument(
        "-f",
        "--format",
        choices=STDOUT_FORMATS,
        default="table",
        help="what to write to standard output (default: table)",
    )
    motifs.add_argument(
        "-o",
        "--output",
        metavar="PATH",
        help="write the catalogue to this file as TSV, ready to edit and reuse",
    )


def _add_codes_parser(subparsers: argparse._SubParsersAction) -> None:
    codes = subparsers.add_parser(
        "codes",
        help="show the residue codes and the pattern notations",
        description=(
            "Show which residues each IUPAC nucleotide code and each amino-acid "
            "ambiguity code stands for, and summarise the pattern notations."
        ),
    )
    codes.add_argument(
        "--type",
        choices=tuple(ALPHABETS),
        help="show only one alphabet",
    )
    codes.add_argument(
        "-f",
        "--format",
        choices=STDOUT_FORMATS,
        default="table",
        help="what to write to standard output (default: table)",
    )


def _same_file(first: str, second: str) -> bool:
    """Whether two paths name one file.

    Comparing the paths as written is not enough: ``out.tsv`` and ``./out.tsv``
    are the same file and do not compare equal, and a hard link is the same
    file under two different names that resolve differently. So existing files
    are compared by identity, through :meth:`~pathlib.Path.samefile`, and
    anything not yet created is compared by its resolved path - which is the
    only thing there is to compare before a file exists.
    """
    one, two = Path(first), Path(second)
    try:
        if one.exists() and two.exists():
            return one.samefile(two)
    except OSError:  # pragma: no cover - an unreadable path is not a match
        pass
    return one.resolve() == two.resolve()


def _inputs_of(args: argparse.Namespace) -> list[str]:
    """Every file this run reads.

    The motif file belongs here as much as the sequence does. It is read
    before the outputs are written, so overwriting it does not corrupt the
    scan - it destroys a file the user wrote by hand, which is worse.
    """
    reads = [path for path in args.input if path and path != "-"]
    motifs = getattr(args, "motifs_file", None)
    if motifs and motifs != "-":
        reads.append(motifs)
    return reads


def _guard_outputs(args: argparse.Namespace) -> None:
    """Refuse to write an output over an input, or two outputs over each other.

    Without the first, ``--output`` aimed at the file being scanned destroys
    the sequence and the input is gone. Without the second, the last write
    silently wins and the user is left with one file where two were asked for.
    """
    outputs = [
        (flag, value)
        for flag, value in (
            ("--output", args.output),
            ("--counts", args.counts),
            ("--report", args.report),
        )
        if value and value != "-"
    ]
    for flag, value in outputs:
        for source in _inputs_of(args):
            if _same_file(source, value):
                also = (
                    ""
                    if Path(source) == Path(value)
                    else f" (the same file as {source})"
                )
                raise ValueError(f"{flag} would overwrite an input file: {value}{also}")
    for index, (flag, value) in enumerate(outputs):
        for other_flag, other in outputs[index + 1:]:
            if _same_file(value, other):
                also = "" if Path(other) == Path(value) else f" and {other}"
                raise ValueError(
                    f"{flag} and {other_flag} would be written to the same file: "
                    f"{value}{also}"
                )


def ineffective_options(args: argparse.Namespace) -> list[str]:
    """Options that will have no effect on what this run writes.

    Reported rather than ignored: a flag that does nothing is a question the
    user asked and did not get an answer to.
    """
    notes: list[str] = []
    shows_table = args.format == "table"
    writes_report = bool(args.report)
    if not shows_table and not writes_report:
        for flag, value, default in (
            ("--limit", args.limit, 50),
            ("--map-width", args.map_width, 48),
        ):
            if value != default:
                notes.append(f"{flag} affects only the table and the report")
        if args.no_map:
            notes.append("--no-map affects only the table and the report")
    if not writes_report and args.bins != DEFAULT_BINS:
        notes.append("--bins affects only the report")
    if not writes_report and not args.sequence:
        if args.line_width != 60:
            notes.append("--line-width affects only the report and --sequence")
        if args.highlight_limit != 20000:
            notes.append("--highlight-limit affects only the report and --sequence")
    if args.format == "none" and not (args.output or args.counts or args.report):
        notes.append(
            "--format none with no --output, --counts or --report writes nothing"
        )
    return notes


def _chosen_motifs(args: argparse.Namespace, alphabet: str) -> list[Motif]:
    """Every motif this run should apply, in the order they were asked for.

    Raises:
        PatternError: if no motif was asked for, if a name is unknown, or if
            two motifs would share a name.
    """
    motifs: list[Motif] = []
    if args.motif:
        names = [name for group in args.motif for name in group]
        try:
            motifs.extend(library.resolve(names))
        except KeyError as error:
            raise PatternError(str(error).strip("\"'")) from None
    if args.all:
        category = args.category or ""
        if category and category not in library.categories():
            raise PatternError(
                f"not a category: {category!r} "
                f"(known: {', '.join(library.categories())})"
            )
        chosen = library.select(category=category, alphabet=alphabet)
        if not chosen:
            where = f" in category {category}" if category else ""
            raise PatternError(
                f"the catalogue holds no {alphabet} motif{where}"
            )
        motifs.extend(chosen)
    for index, text in enumerate(args.pattern or (), start=1):
        motifs.append(named(text, "iupac", alphabet, index))
    for index, text in enumerate(args.prosite or (), start=1):
        motifs.append(named(text, "prosite", alphabet, index))
    for index, text in enumerate(args.regex or (), start=1):
        motifs.append(named(text, "regex", alphabet, index))
    if args.motifs_file:
        motifs.extend(library.load(args.motifs_file))
    if not motifs:
        raise PatternError(
            "no motif was asked for; name one with --motif, give a pattern with "
            "--pattern, --prosite or --regex, read a file with --motifs-file, or "
            "apply the whole catalogue with --all"
        )
    if args.category and not args.all:
        raise PatternError("--category applies to --all, which was not given")
    return library.unique(motifs)


def _read_inputs(paths: Sequence[str]) -> list[Record]:
    """Read every record of every input, labelling each with its file."""
    records: list[Record] = []
    several = len(paths) > 1
    for path in paths:
        label = path if several else ""
        records.extend(read_records(path, label=label))
    if not records:
        raise SequenceError("the input holds no records")
    return records


def _decide_alphabet(records: Sequence[Record], asked: str) -> str:
    """Settle which alphabet the records are written in.

    Raises:
        SequenceError: if the records disagree and none was forced. Scanning a
            mixed file against one alphabet would silently read half of it as
            the wrong kind of sequence.
    """
    if asked != "auto":
        return asked
    judged = {detect(record.sequence) for record in records} - {None}
    if len(judged) > 1:
        shown = ", ".join(sorted(str(name) for name in judged))
        raise SequenceError(
            f"the records are not all the same kind of sequence ({shown}); scan "
            "them separately, or force one alphabet with --type"
        )
    return judged.pop() if judged else NUCLEOTIDE.name


def _compile_all(motifs: Sequence[Motif], alphabet: str) -> list[Compiled]:
    """Compile every motif against the alphabet of the input.

    Raises:
        PatternError: naming the motifs written in the other alphabet, which is
            the one mistake that is worth catching before any sequence is read.
        """
    target = get_alphabet(alphabet)
    wrong = [motif.name for motif in motifs if motif.alphabet != alphabet]
    if wrong:
        shown = ", ".join(wrong)
        raise PatternError(
            f"the input is {alphabet} sequence but these motifs are not: {shown}"
        )
    return [compile_motif(motif, target) for motif in motifs]


def run_scan(
    records: Sequence[Record],
    compiled: Sequence[Compiled],
    *,
    mismatches: int = 0,
    strand: str = "both",
    circular: bool = False,
    overlapping: bool = True,
    include_uncertain: bool = False,
    bins: int = DEFAULT_BINS,
    sort: str = "position",
    name: str = "input",
    parameters: Sequence[tuple[str, str]] = (),
) -> Scan:
    """Scan every record for every motif."""
    names = [motif.name for motif in compiled]
    alphabet = compiled[0].alphabet if compiled else NUCLEOTIDE
    scanned: list[RecordScan] = []
    for record in records:
        composition = record.composition()
        ambiguous = record.ambiguous(alphabet)
        hits: list[Hit] = []
        counts: list[MotifCount] = []
        for motif in compiled:
            found = find(
                record.sequence,
                motif,
                mismatches=mismatches,
                strand=strand,
                circular=circular,
                overlapping=overlapping,
                record=record.identifier,
                include_uncertain=include_uncertain,
            )
            uncertain = count_uncertain(
                record.sequence,
                motif,
                mismatches=mismatches,
                strand=strand,
                circular=circular,
            )
            # A window carrying an ambiguity code cannot yield a hit, so the
            # expectation is taken over the windows that could - unless those
            # windows were asked for as hits, in which case both sides count
            # them again.
            windows = None
            span = motif.fixed_length
            if span and ambiguous and not include_uncertain:
                windows = assessable_windows(
                    record.sequence, span, motif.alphabet, circular
                )
            hits.extend(found)
            counts.append(
                tally(
                    motif,
                    found,
                    len(record.sequence),
                    composition,
                    mismatches=mismatches,
                    strand=strand,
                    circular=circular,
                    uncertain=uncertain,
                    overlapping=overlapping,
                    windows=windows,
                )
            )
        if sort == "motif":
            hits.sort(key=lambda hit: (names.index(hit.motif), hit.start, hit.strand))
        else:
            hits.sort(key=lambda hit: (hit.start, hit.strand, hit.motif))
        scanned.append(
            RecordScan(
                identifier=record.identifier,
                description=record.description,
                source=record.source,
                residues=len(record.sequence),
                alphabet=alphabet.name,
                circular=circular,
                ambiguous=ambiguous,
                hits=tuple(hits),
                counts=tuple(counts),
                distribution=distribute(hits, names, len(record.sequence), bins),
            )
        )
    return Scan(
        name=name,
        records=tuple(scanned),
        motifs=tuple(compiled),
        parameters=tuple(parameters),
        mismatches=mismatches,
        overlapping=overlapping,
        include_uncertain=include_uncertain,
    )


def _parameters(args: argparse.Namespace, alphabet: str, motifs: int) -> list[tuple[str, str]]:
    """The settings the scan ran under, for the report."""
    rows = [
        ("Input", ", ".join(args.input)),
        ("Alphabet", alphabet + (" (judged from composition)" if args.type == "auto" else " (given)")),
        ("Motifs", thousands(motifs)),
        ("Mismatches allowed", thousands(args.mismatches)),
    ]
    if get_alphabet(alphabet).has_strands:
        rows.append(("Strand", args.strand))
        rows.append(("Topology", "circular" if args.circular else "linear"))
    rows.append(
        (
            "Overlapping hits",
            "leftmost of each run only, so no expectation is reported"
            if args.no_overlapping
            else "all reported",
        )
    )
    if args.include_uncertain:
        rows.append(("Uncertain windows", "listed as hits"))
    return rows


def _hit_rows(scan: Scan) -> list[list[str]]:
    """Every hit as a row of the hits file."""
    patterns = {
        motif.name: (motif.motif.pattern, motif.motif.syntax) for motif in scan.motifs
    }
    rows = []
    for record in scan.records:
        for hit in record.hits:
            pattern, syntax = patterns.get(hit.motif, ("", ""))
            rows.append(
                [
                    record.label,
                    hit.motif,
                    str(hit.start),
                    str(hit.end),
                    hit.strand,
                    str(hit.length),
                    hit.matched,
                    str(hit.mismatches),
                    str(hit.uncertain),
                    pattern,
                    syntax,
                ]
            )
    return rows


def _count_rows(scan: Scan) -> list[list[str]]:
    """Every per-motif count as a row of the counts file."""
    patterns = {
        motif.name: (motif.motif.pattern, motif.motif.syntax) for motif in scan.motifs
    }
    rows = []
    for record in scan.records:
        for count in record.counts:
            pattern, syntax = patterns.get(count.motif, ("", ""))
            rows.append(
                [
                    record.label,
                    count.motif,
                    pattern,
                    syntax,
                    str(count.residues),
                    str(count.hits),
                    str(count.forward),
                    str(count.reverse),
                    str(count.palindromic),
                    str(count.unstranded),
                    str(count.exact),
                    str(count.with_mismatches),
                    str(count.listed_uncertain),
                    str(count.uncertain),
                    "yes" if count.assessed else "no",
                    f"{count.per_kb:.6g}",
                    "" if count.expected is None else f"{count.expected:.6g}",
                    "" if count.enrichment is None else f"{count.enrichment:.6g}",
                ]
            )
    return rows


def _write_tsv(
    handle: IO[str], columns: Sequence[str], rows: Sequence[Sequence[str]]
) -> None:
    """Write a header and rows, tab-separated."""
    handle.write("\t".join(columns) + "\n")
    for row in rows:
        handle.write("\t".join(row) + "\n")


def _as_json(scan: Scan) -> dict[str, object]:
    """The whole scan as a JSON-ready structure."""
    return {
        "motif_scan": __version__,
        "name": scan.name,
        "alphabet": scan.alphabet,
        "mismatches": scan.mismatches,
        "parameters": dict(scan.parameters),
        "motifs": [
            {
                "name": motif.name,
                "pattern": motif.motif.pattern,
                "syntax": motif.motif.syntax,
                "alphabet": motif.motif.alphabet,
                "category": motif.motif.category,
                "description": motif.motif.description,
                "source": motif.motif.source,
                "length": motif.fixed_length,
                "palindromic": motif.palindromic,
                "supports_mismatches": motif.supports_mismatches,
            }
            for motif in scan.motifs
        ],
        "records": [
            {
                "identifier": record.identifier,
                "description": record.description,
                "source": record.source,
                "residues": record.residues,
                "alphabet": record.alphabet,
                "topology": record.topology,
                "ambiguous": record.ambiguous,
                "counts": [
                    {
                        "motif": count.motif,
                        "hits": count.hits,
                        "forward": count.forward,
                        "reverse": count.reverse,
                        "palindromic": count.palindromic,
                        "unstranded": count.unstranded,
                        "exact": count.exact,
                        "with_mismatches": count.with_mismatches,
                        "listed_uncertain": count.listed_uncertain,
                        "uncertain": count.uncertain,
                        "uncertainty_assessed": count.assessed,
                        "per_kb": count.per_kb,
                        "expected": count.expected,
                        "observed_over_expected": count.enrichment,
                    }
                    for count in record.counts
                ],
                "hits": [
                    {
                        "motif": hit.motif,
                        "start": hit.start,
                        "end": hit.end,
                        "strand": hit.strand,
                        "length": hit.length,
                        "match": hit.matched,
                        "mismatches": hit.mismatches,
                        "uncertain": hit.uncertain,
                        "spans_origin": hit.spans_origin,
                    }
                    for hit in record.hits
                ],
                "distribution": {
                    "bins": record.distribution.bins,
                    "counts": {
                        motif: list(counts)
                        for motif, counts in record.distribution.counts
                    },
                },
            }
            for record in scan.records
        ],
    }


def _count_table(
    counts: Sequence[MotifCount],
    motifs: Sequence[Compiled],
    strands: bool = True,
) -> list[str]:
    """The per-motif counts, as aligned columns.

    The strand columns are left out for a sequence that has only one strand.
    Printing a column of zeroes headed by a minus sign for a protein would
    imply a reverse strand that does not exist.
    """
    patterns = {motif.name: motif.motif.pattern for motif in motifs}
    columns = f"{'+':>5} {'-':>5} " if strands else ""
    lines = [
        f"  {'MOTIF':<16} {'PATTERN':<18} {'HITS':>6} {columns}"
        f"{'EXACT':>6} {'PER KB':>8} {'EXPECTED':>9} {'OBS/EXP':>8}"
    ]
    for count in counts:
        forward = count.forward + count.palindromic + count.unstranded
        split = (
            f"{thousands(forward):>5} {thousands(count.reverse):>5} "
            if strands
            else ""
        )
        lines.append(
            f"  {truncate(count.motif, 16):<16} "
            f"{truncate(patterns.get(count.motif, ''), 18):<18} "
            f"{thousands(count.hits):>6} {split}"
            f"{thousands(count.exact):>6} "
            f"{number(count.per_kb):>8} {number(count.expected, 1):>9} "
            f"{ratio(count.enrichment):>8}"
        )
    return lines


def _hit_table(hits: Sequence[Hit], limit: int) -> list[str]:
    """The hits themselves, as aligned columns."""
    shown = hits[:limit] if limit else hits
    lines = [
        f"  {'START':>9} {'END':>9} {'STRAND':<7} {'MOTIF':<16} {'MATCH':<24} MM"
    ]
    for hit in shown:
        note = " (spans join)" if hit.spans_origin else ""
        lines.append(
            f"  {thousands(hit.start):>9} {thousands(hit.end):>9} "
            f"{hit.strand:<7} {truncate(hit.motif, 16):<16} "
            f"{truncate(hit.matched, 24):<24} {hit.mismatches}{note}"
        )
    if limit and len(hits) > limit:
        lines.append(
            f"  ... {plural(len(hits) - limit, 'further hit')} not shown; "
            "raise --limit or write them with --output"
        )
    return lines


def _render_scan(
    scan: Scan,
    sequences: Sequence[str],
    args: argparse.Namespace,
) -> list[str]:
    """The whole scan, as the terminal table.

    The residues are taken by position rather than by record name, because two
    records of one file may share an identifier.
    """
    lines: list[str] = []
    if not scan.overlapping and scan.motifs:
        lines.append(
            "  note       no expected count is reported: suppressing overlapping"
        )
        lines.append(
            "             hits reports the leftmost of each run, while the"
        )
        lines.append(
            "             expectation counts every matching window"
        )
        lines.append("")
    for index, record in enumerate(scan.records):
        unit = "aa" if record.alphabet == "protein" else "bp"
        head = (
            f"> {record.label}  {thousands(record.residues)} {unit}  "
            f"{record.alphabet}  {record.topology}"
        )
        lines.append(head)
        if record.description:
            lines.append(f"  {truncate(record.description, 72)}")
        lines.append("")
        lines.extend(
            _count_table(
                record.counts, scan.motifs, strands=record.alphabet != "protein"
            )
        )

        uncertain = [count for count in record.counts if count.uncertain]
        if uncertain:
            listed = ", ".join(
                f"{count.motif} ({thousands(count.uncertain)})" for count in uncertain
            )
            lines.append("")
            where = (
                "listed among the hits below"
                if scan.include_uncertain
                else "not counted as hits"
            )
            lines.append(
                f"  uncertain  positions that could not be called, because the "
                f"sequence is ambiguous there: {listed} ({where})"
            )

        if not args.no_map and record.residues and record.hits:
            lines.append("")
            width = args.map_width
            for motif in scan.motifs:
                found = record.hits_for(motif.name)
                lines.append(
                    f"  {truncate(motif.name, 16):<16} "
                    f"{track(found, record.residues, width)}  {thousands(len(found))}"
                )
            lines.extend(
                f"  {' ' * 17}{line}"
                for line in ruler(record.residues, width, record.circular)
            )

        if args.sequence and record.residues <= args.highlight_limit:
            sequence = sequences[index] if index < len(sequences) else ""
            if sequence:
                lines.append("")
                marks = legend([motif.name for motif in scan.motifs])
                lines.extend(
                    blocks(sequence, record.hits, args.line_width, marks)
                )
                lines.append("")
                lines.extend(
                    line
                    for name, _ in marks.symbols
                    for line in [f"  {marks.of(name)}  {name}"]
                    if record.hits_for(name)
                )
        elif args.sequence:
            lines.append("")
            lines.append(
                f"  sequence   not shown: {thousands(record.residues)} {unit} is "
                f"past the --highlight-limit of {thousands(args.highlight_limit)}"
            )

        if record.hits:
            lines.append("")
            lines.extend(_hit_table(record.hits, args.limit))
        else:
            lines.append("")
            lines.append("  no hits")
        lines.append("")
    return lines


def _motif_table(motifs: Sequence[Motif]) -> list[str]:
    """The catalogue, as aligned columns."""
    lines = [
        f"{'MOTIF':<16} {'PATTERN':<26} {'SYNTAX':<8} {'ALPHABET':<11} "
        f"{'CATEGORY':<21} DESCRIPTION"
    ]
    for motif in motifs:
        lines.append(
            f"{truncate(motif.name, 16):<16} {truncate(motif.pattern, 26):<26} "
            f"{motif.syntax:<8} {motif.alphabet:<11} {motif.category:<21} "
            f"{motif.description}"
        )
    return lines


def _motif_rows(motifs: Sequence[Motif]) -> list[list[str]]:
    """The catalogue as rows of a motif file, ready to edit and read back."""
    rows = []
    for motif in motifs:
        try:
            length = compile_motif(motif).fixed_length
        except PatternError:  # pragma: no cover - the catalogue compiles
            length = None
        rows.append(
            [
                motif.name,
                motif.pattern,
                motif.syntax,
                motif.alphabet,
                motif.category,
                "" if length is None else str(length),
                motif.description,
                motif.source,
            ]
        )
    return rows


def _run_motifs(args: argparse.Namespace, out: IO[str]) -> int:
    """List the catalogue."""
    if args.category and args.category not in library.categories():
        raise PatternError(
            f"not a category: {args.category!r} "
            f"(known: {', '.join(library.categories())})"
        )
    motifs = library.select(
        category=args.category or "",
        alphabet=args.type or "",
        search=args.search or "",
    )
    if not motifs:
        raise PatternError("no catalogue motif matches those filters")
    if args.format == "table":
        out.write("\n".join(_motif_table(motifs)) + "\n")
    elif args.format == "tsv":
        _write_tsv(out, _MOTIF_COLUMNS, _motif_rows(motifs))
    elif args.format == "json":
        json.dump(
            [
                dict(zip(_MOTIF_COLUMNS, row))
                for row in _motif_rows(motifs)
            ],
            out,
            indent=2,
        )
        out.write("\n")
    if args.output:
        with open(args.output, "w", encoding="utf-8", newline="\n") as handle:
            _write_tsv(handle, _MOTIF_COLUMNS, _motif_rows(motifs))
    return 0


def _code_rows(alphabet: str | None) -> list[list[str]]:
    """Every code of the alphabets asked for, and what it stands for."""
    rows = []
    for name, codes in (
        (NUCLEOTIDE.name, NUCLEOTIDE_CODES),
        (PROTEIN.name, PROTEIN_CODES),
    ):
        if alphabet and name != alphabet:
            continue
        for code, residues in codes.items():
            rows.append(
                [name, code, "".join(sorted(residues)), str(len(residues))]
            )
    return rows


_NOTATION = """
Pattern notations

  iupac     a string of alphabet codes, one per position: TATAWAWR, CANNTG.
            Fixed length, so mismatches can be allowed and both strands
            searched position by position.

  prosite   PROSITE syntax, tokens joined by hyphens:
              A        this residue            [ST]     any of these
              {P}      any except these        x        any residue
              x(2)     exactly twice           x(2,4)   two to four times
              <        anchored to the start    >       anchored to the end
            A pattern whose repeats are all exact counts has a fixed length and
            behaves like an iupac one. One with a range does not, so it cannot
            be searched with mismatches.

  regex     a Python regular expression, applied as written. No mismatches, and
            a pattern that can match an empty string is refused.

Matching is asymmetric by design. An ambiguity code in the pattern widens what
the pattern accepts; an ambiguity code in the sequence means the residue is not
known, and those positions are reported as uncertain rather than as hits.
""".strip()


def _run_codes(args: argparse.Namespace, out: IO[str]) -> int:
    """Show the residue codes and the pattern notations."""
    rows = _code_rows(args.type)
    if args.format == "table":
        lines = [f"{'ALPHABET':<11} {'CODE':<5} {'STANDS FOR':<21} COUNT"]
        for alphabet, code, residues, count in rows:
            lines.append(
                f"{alphabet:<11} {code:<5} {truncate(residues, 21):<21} {count}"
            )
        out.write("\n".join(lines) + "\n\n")
        out.write(_NOTATION + "\n")
    elif args.format == "tsv":
        _write_tsv(out, ("alphabet", "code", "stands_for", "residues"), rows)
    elif args.format == "json":
        json.dump(
            [
                dict(zip(("alphabet", "code", "stands_for", "residues"), row))
                for row in rows
            ],
            out,
            indent=2,
        )
        out.write("\n")
    return 0


def _run_scan(args: argparse.Namespace, out: IO[str], err: IO[str]) -> int:
    """Scan the inputs."""
    _guard_outputs(args)
    for note in ineffective_options(args):
        err.write(f"motif-scan: note: {note}\n")

    records = _read_inputs(args.input)
    alphabet = _decide_alphabet(records, args.type)
    target = get_alphabet(alphabet)
    for record in records:
        check_alphabet(record, target)
    if args.circular and not target.has_strands:
        raise PatternError(
            f"--circular describes a closed molecule of nucleic acid, and the "
            f"input is {alphabet}"
        )

    motifs = _chosen_motifs(args, alphabet)
    compiled = _compile_all(motifs, alphabet)

    sources = ", ".join(dict.fromkeys(args.input))
    name = sources if len(args.input) == 1 else f"{plural(len(args.input), 'input')}"
    scan = run_scan(
        records,
        compiled,
        mismatches=args.mismatches,
        strand=args.strand,
        circular=args.circular,
        overlapping=not args.no_overlapping,
        include_uncertain=args.include_uncertain,
        bins=args.bins,
        sort=args.sort,
        name=name,
        parameters=_parameters(args, alphabet, len(compiled)),
    )
    sequences = [record.sequence for record in records]

    if args.format == "table":
        out.write("\n".join(_render_scan(scan, sequences, args)) + "\n")
    elif args.format == "tsv":
        _write_tsv(out, _HIT_COLUMNS, _hit_rows(scan))
    elif args.format == "json":
        json.dump(_as_json(scan), out, indent=2)
        out.write("\n")

    if args.output:
        with open(args.output, "w", encoding="utf-8", newline="\n") as handle:
            _write_tsv(handle, _HIT_COLUMNS, _hit_rows(scan))
    if args.counts:
        with open(args.counts, "w", encoding="utf-8", newline="\n") as handle:
            _write_tsv(handle, _COUNT_COLUMNS, _count_rows(scan))
    if args.report:
        html = render_html(
            scan,
            sequences,
            __version__,
            width=args.line_width,
            highlight_limit=args.highlight_limit,
            limit=args.limit,
            map_width=args.map_width,
        )
        with open(args.report, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(html)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command line.

    Returns the process exit status: 0 on success, 1 for anything the user can
    fix - an unreadable file, a pattern that cannot be compiled, an option that
    contradicts the input.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    out = sys.stdout
    err = sys.stderr
    try:
        if args.command == "scan":
            return _run_scan(args, out, err)
        if args.command == "motifs":
            return _run_motifs(args, out)
        return _run_codes(args, out)
    except (SequenceError, PatternError, ValueError, KeyError) as error:
        err.write(f"motif-scan: {str(error).strip(chr(39)).strip(chr(34))}\n")
        return 1
    except FileNotFoundError as error:
        err.write(f"motif-scan: {error}\n")
        return 1
    except IsADirectoryError as error:
        err.write(f"motif-scan: {error}\n")
        return 1
    except BrokenPipeError:  # pragma: no cover - depends on the consumer
        return 0
