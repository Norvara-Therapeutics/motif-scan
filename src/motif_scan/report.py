"""Self-contained HTML reporting.

The document embeds its own style rules, uses no JavaScript and loads nothing
from the network, so a report stays readable from a local file, an email
attachment or a strict content-security policy years later. Output is
byte-deterministic for a given scan: nothing in it depends on the clock or the
filesystem beyond the names passed in.

The hits table is filtered without scripting. Each control is a checkbox
standing before the table, and each row carries a class for its motif, its
strand and whether every position matched; a rule per control hides the rows it
switches off. The checkbox itself is moved off-screen and drawn by its label,
which is why every input is immediately followed by the label that belongs to
it and both sit as direct children of the section: the adjacent-sibling rule
that draws the box, and the general-sibling rule that reaches the table, both
depend on that order. What this buys is filtering that works offline and
prints; what it costs is that the counts above the table are the counts of the
whole scan and do not follow the filter, which the report says in as many
words.
"""

from __future__ import annotations

from html import escape
from typing import Mapping, Sequence

from motif_scan.mapview import ruler, track
from motif_scan.patterns import Compiled
from motif_scan.search import Hit
from motif_scan.summary import Distribution, MotifCount, RecordScan, Scan
from motif_scan.text import number, plural, ratio, thousands, truncate

#: How many motifs get their own highlight colour. Beyond this the colours
#: repeat, and the legend is what distinguishes them.
_COLOURS = 8

_STYLE = """
html { -webkit-text-size-adjust: 100%; }
body { font: 16px/1.55 Georgia, 'Times New Roman', serif; color: #1a1a1a;
       background: #fff; margin: 0 auto; padding: 3rem 1.25rem 4rem;
       max-width: 52rem; }
h1 { font-size: 1.5rem; font-weight: normal; margin: 0 0 .25rem; }
h2, h3 { font-size: .8rem; font-variant: small-caps; letter-spacing: .08em;
     font-weight: normal; margin: 2.5rem 0 .6rem; border-bottom: 1px solid #1a1a1a;
     padding-bottom: .2rem; }
h3 { margin-top: 1.8rem; border-bottom: 1px solid #c8c8c8; }
p { margin: .5rem 0; }
.lede { color: #555; font-size: .9rem; margin-bottom: .5rem; }
.note { color: #555; font-size: .85rem; }
code, .seq, td.num, th.num, .mono { font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace; }
table { border-collapse: collapse; width: 100%; margin: .6rem 0 .2rem;
        font-size: .9rem; }
th { text-align: left; font-weight: normal; font-variant: small-caps;
     letter-spacing: .04em; border-bottom: 1px solid #1a1a1a; padding: .25rem .5rem .25rem 0; }
td { padding: .22rem .5rem .22rem 0; border-bottom: 1px solid #e6e6e6;
     vertical-align: baseline; }
td.num, th.num { text-align: right; padding-right: 1rem; font-size: .85rem; }
tr:last-child td { border-bottom: none; }
dl.metrics { display: grid; grid-template-columns: max-content 1fr; gap: .15rem 1.2rem;
             margin: .6rem 0; font-size: .9rem; }
dl.metrics dt { color: #555; }
dl.metrics dd { margin: 0; font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace; }
.bar { display: block; height: 2px; background: #1a1a1a; margin: .3rem 0 .1rem; }
.tag { font-variant: small-caps; letter-spacing: .04em; font-size: .8rem; color: #0f6f6a; }
.tag.warn { color: #8a5a00; }
details { margin: 1.2rem 0; }
summary { cursor: pointer; font-size: .9rem; }
summary h2 { display: inline-block; margin: .4rem 0 .2rem; }
pre.map { font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace;
     font-size: .78rem; line-height: 1.35; overflow-x: auto; background: #fafafa;
     border: 1px solid #e6e6e6; padding: .7rem .8rem; margin: .6rem 0; }
.chart { margin: .6rem 0; }
.chart .row { display: flex; align-items: flex-end; gap: 1px; height: 2.2rem;
     border-bottom: 1px solid #c8c8c8; }
.chart .row span { flex: 1 1 0; background: #1a1a1a; min-height: 0; }
.chart .row span.zero { background: #e6e6e6; height: 1px; }
.chart .label { font-size: .8rem; color: #555; margin: .5rem 0 .1rem; }
.chart .scale { display: flex; justify-content: space-between; font-size: .75rem;
     color: #555; font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace; }
.seqmap { background: #fafafa; border: 1px solid #e6e6e6; padding: .7rem .8rem;
     margin: .6rem 0; overflow-x: auto; }
.seqline { font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace;
     font-size: .78rem; line-height: 1.5; white-space: pre; }
.seqline .at { color: #999; user-select: none; }
.seqline:target { background: #fff3cd; }
mark { background: none; padding: 0; }
mark.m0 { background: #d8e9f7; }
mark.m1 { background: #f7e3d8; }
mark.m2 { background: #dcf0da; }
mark.m3 { background: #efdcf0; }
mark.m4 { background: #f7f1cf; }
mark.m5 { background: #d8f0ee; }
mark.m6 { background: #f0dade; }
mark.m7 { background: #e4e2f5; }
mark.multi { box-shadow: inset 0 -2px 0 #8a5a00; }
.key { font-size: .8rem; color: #555; margin: .4rem 0; }
.key span { margin-right: .9rem; white-space: nowrap; }
.hits { margin: .6rem 0; }
.hits > .controls { font-size: .8rem; color: #555; font-variant: small-caps;
     letter-spacing: .04em; margin-right: .4rem; }
.hits > label { font-size: .8rem; color: #555; margin-right: .8rem;
     white-space: nowrap; cursor: pointer; }
.filter { position: absolute; left: -9999px; }
.filter + label::before { content: "\\2611\\00a0"; }
.filter:not(:checked) + label::before { content: "\\2610\\00a0"; }
.filter:not(:checked) + label { color: #aaa; }
footer { margin-top: 3rem; padding-top: .5rem; border-top: 1px solid #e6e6e6;
         color: #555; font-size: .8rem; }
@media print { body { padding: 0; max-width: none; }
    h2, h3 { page-break-after: avoid; }
    .hits > label, .hits > .controls { display: none; }
    details { display: block; } summary { display: none; } }
"""

#: Row classes by strand label, so that a CSS rule can name a strand without
#: putting a symbol in a selector.
_STRAND_CLASS = {"+": "s-plus", "-": "s-minus", "=": "s-both", ".": "s-single"}

#: How a strand is written in the hits table.
_STRAND_NAME = {"+": "forward", "-": "reverse", "=": "both", ".": "single"}


def _definitions(rows: Sequence[tuple[str, str]]) -> list[str]:
    parts = ['<dl class="metrics">']
    for name, value in rows:
        parts.append(f"<dt>{escape(name)}</dt><dd>{escape(value)}</dd>")
    parts.append("</dl>")
    return parts


def _motif_index(scan: Scan) -> dict[str, int]:
    """Position of each motif in the order they were applied."""
    return {compiled.name: index for index, compiled in enumerate(scan.motifs)}


def _summary_table(
    counts: Sequence[MotifCount],
    motifs: Sequence[Compiled],
    unit: str,
    overlapping: bool = True,
    ambiguity: bool = False,
) -> list[str]:
    """Hits per motif, with a bar for the count: the frequency breakdown."""
    by_name = {compiled.name: compiled for compiled in motifs}
    busiest = max((count.hits for count in counts), default=0)
    parts = [
        "<table>",
        "<tr><th>Motif</th><th>Pattern</th><th>Syntax</th>"
        '<th class="num">Hits</th><th></th><th class="num">+</th>'
        '<th class="num">&minus;</th><th class="num">Exact</th>'
        '<th class="num">Per kb</th><th class="num">Expected</th>'
        '<th class="num">Obs/exp</th></tr>',
    ]
    for count in counts:
        compiled = by_name.get(count.motif)
        pattern = compiled.motif.pattern if compiled else ""
        syntax = compiled.motif.syntax if compiled else ""
        width = 100.0 * count.hits / busiest if busiest else 0.0
        bar = f'<span class="bar" style="width: {width:.1f}%"></span>'
        forward = count.forward + count.palindromic + count.unstranded
        parts.append(
            f"<tr><td>{escape(count.motif)}</td>"
            f'<td class="mono">{escape(pattern)}</td>'
            f"<td>{escape(syntax)}</td>"
            f'<td class="num">{thousands(count.hits)}</td>'
            f"<td>{bar}</td>"
            f'<td class="num">{thousands(forward)}</td>'
            f'<td class="num">{thousands(count.reverse)}</td>'
            f'<td class="num">{thousands(count.exact)}</td>'
            f'<td class="num">{number(count.per_kb)}</td>'
            f'<td class="num">{number(count.expected, 1)}</td>'
            f'<td class="num">{ratio(count.enrichment)}</td></tr>'
        )
    parts.append("</table>")
    parts.append(
        '<p class="note">The + column counts hits on the forward strand, and '
        "also the hits of a motif that reads identically on both strands, which "
        f"are reported once rather than twice. Per kb is hits per thousand {escape(unit)}. "
        "Expected is how many hits a record of this length and residue "
        "composition would carry by chance, with the mismatch budget taken into "
        "account, under a model of independent residues: it is a yardstick for "
        "reading the count beside it, not a test, and real sequence breaks the "
        "model in both directions. It is left blank for a pattern with no fixed "
        "length, which has no positions to compute it over.</p>"
    )
    if ambiguity:
        parts.append(
            '<p class="note">This record carries ambiguity codes. A window '
            "holding one cannot yield a hit, so it is left out of the "
            "expectation as well as out of the count, and both figures are over "
            "the windows that could be called.</p>"
        )
    if not overlapping:
        parts.append(
            '<p class="note">No expectation is reported at all for this scan, '
            "because overlapping hits were suppressed. The expectation counts "
            "every window a motif matches, while this scan reports only the "
            "leftmost of each run of overlapping hits, so the two figures would "
            "no longer be answering the same question and their ratio would "
            "move with that choice rather than with the sequence.</p>"
        )
    return parts


def _distribution_chart(
    distribution: Distribution, motifs: Sequence[str], unit: str
) -> list[str]:
    """One row of bars per motif, binned along the record."""
    if not distribution.bins or not distribution.busiest:
        return []
    tallest = distribution.busiest
    parts = ['<div class="chart">']
    present = {motif: counts for motif, counts in distribution.counts}
    for motif in motifs:
        counts = present.get(motif, ())
        if not counts or not any(counts):
            continue
        parts.append(f'<div class="label">{escape(motif)}</div>')
        parts.append('<div class="row">')
        for index, count in enumerate(counts):
            first, last = distribution.span(index)
            title = (
                f"{thousands(count)} in {thousands(first)}..{thousands(last)}"
            )
            if not count:
                parts.append(f'<span class="zero" title="{escape(title)}"></span>')
            else:
                height = 100.0 * count / tallest
                parts.append(
                    f'<span style="height: {height:.1f}%" title="{escape(title)}"></span>'
                )
        parts.append("</div>")
    parts.append(
        f'<div class="scale"><span>1</span><span>{thousands(distribution.residues)}'
        f" {escape(unit)}</span></div>"
    )
    parts.append("</div>")
    parts.append(
        f'<p class="note">Each row is one motif, divided into '
        f"{thousands(distribution.bins)} equal bins along the record, and every "
        "row is drawn to the same scale, so rows can be compared with one "
        f"another. The tallest bar is {plural(tallest, 'hit')}.</p>"
    )
    return parts


def _runs(hit: Hit, length: int) -> list[tuple[int, int]]:
    """The 0-based half-open spans a hit covers, splitting one that wraps."""
    first = hit.start - 1
    if not hit.spans_origin:
        return [(first, min(first + hit.length, length))]
    return [(first, length), (0, hit.end)]


def _cover(
    sequence: str, hits: Sequence[Hit], order: Mapping[str, int]
) -> list[tuple[int, ...]]:
    """For each residue, the motifs covering it, as indices in motif order."""
    covered: list[tuple[int, ...]] = [()] * len(sequence)
    for hit in hits:
        index = order.get(hit.motif)
        if index is None:
            continue
        for start, stop in _runs(hit, len(sequence)):
            for position in range(start, stop):
                covered[position] = covered[position] + (index,)
    return covered


def _sequence_map(
    sequence: str,
    hits: Sequence[Hit],
    order: Mapping[str, int],
    prefix: str,
    width: int,
) -> list[str]:
    """The residues themselves, with every hit highlighted in place.

    Residues covered by more than one motif are drawn in the colour of the
    first and underlined, rather than being given a colour of their own: an
    overlap is common enough that inventing a colour for every combination
    would exhaust both the palette and the reader.
    """
    covered = _cover(sequence, hits, order)
    gutter = len(thousands(len(sequence)))
    lines = ['<div class="seqmap">']
    for start in range(0, len(sequence), width):
        stop = min(start + width, len(sequence))
        pieces: list[str] = []
        index = start
        while index < stop:
            here = covered[index]
            run = index
            while run < stop and covered[run] == here:
                run += 1
            text = escape(sequence[index:run])
            if here:
                classes = f"m{here[0] % _COLOURS}"
                if len(here) > 1:
                    classes += " multi"
                pieces.append(f'<mark class="{classes}">{text}</mark>')
            else:
                pieces.append(text)
            index = run
        label = thousands(start + 1).rjust(gutter)
        lines.append(
            f'<div class="seqline" id="{prefix}-{start + 1}">'
            f'<span class="at">{label}</span> {"".join(pieces)}</div>'
        )
    lines.append("</div>")
    return lines


def _key(motifs: Sequence[str], order: Mapping[str, int]) -> list[str]:
    """Which colour stands for which motif."""
    if not motifs:
        return []
    parts = ['<p class="key">']
    for motif in motifs:
        index = order.get(motif, 0) % _COLOURS
        parts.append(f'<span><mark class="m{index}">&nbsp;&nbsp;</mark> {escape(motif)}</span>')
    parts.append(
        '<span><mark class="multi">&nbsp;&nbsp;</mark> more than one motif</span></p>'
    )
    return parts


def _controls(
    prefix: str,
    motifs: Sequence[str],
    order: Mapping[str, int],
    strands: Sequence[str],
    uncertain: bool = False,
) -> tuple[list[str], list[str]]:
    """Filter checkboxes for one hits table, and the rules that drive them."""
    controls: list[str] = ['<span class="controls">Show</span>']
    rules: list[str] = []

    def switch(box: str, caption: str, row: str) -> None:
        # The label follows its input directly: the rule that draws the box is
        # an adjacent-sibling one, and separating them would leave a control
        # that cannot be seen or clicked.
        controls.append(
            f'<input type="checkbox" class="filter" id="{box}" checked>'
            f'<label for="{box}">{caption}</label>'
        )
        rules.append(f"#{box}:not(:checked) ~ table tr.{row} {{ display: none; }}")

    for motif in motifs:
        index = order[motif]
        switch(f"{prefix}-m{index}", escape(motif), f"m-{index}")
    for label in strands:
        css = _STRAND_CLASS[label]
        switch(f"{prefix}-{css}", escape(_STRAND_NAME[label]), css)
    switch(f"{prefix}-inexact", "with mismatches", "inexact")
    if uncertain:
        switch(f"{prefix}-uncertain", "uncertain", "uncertain")
    return controls, rules


def _hits_table(
    hits: Sequence[Hit],
    order: Mapping[str, int],
    prefix: str,
    width: int,
    limit: int,
    linked: bool,
) -> list[str]:
    """Every hit, one row each, classed so that CSS can filter them.

    A start position is a link into the sequence map when that map was drawn.
    For a record too long to show in full there is nothing to link to, and a
    link that goes nowhere is worse than plain text.
    """
    parts = [
        "<table>",
        '<tr><th class="num">Start</th><th class="num">End</th><th>Strand</th>'
        "<th>Motif</th><th>Match</th><th class=\"num\">Mismatches</th></tr>",
    ]
    for hit in hits[:limit] if limit else hits:
        index = order.get(hit.motif, 0)
        classes = [f"m-{index}", _STRAND_CLASS[hit.strand]]
        if not hit.exact:
            classes.append("inexact")
        if not hit.certain:
            classes.append("uncertain")
        anchor = ((hit.start - 1) // width) * width + 1
        note = " \u00b7 spans join" if hit.spans_origin else ""
        if not hit.certain:
            note += " \u00b7 uncertain"
        position = (
            f'<a href="#{prefix}-{anchor}">{thousands(hit.start)}</a>'
            if linked
            else thousands(hit.start)
        )
        parts.append(
            f'<tr class="{" ".join(classes)}">'
            f'<td class="num">{position}</td>'
            f'<td class="num">{thousands(hit.end)}{escape(note)}</td>'
            f"<td>{escape(_STRAND_NAME[hit.strand])}</td>"
            f"<td>{escape(hit.motif)}</td>"
            f'<td class="mono">{escape(truncate(hit.matched, 48))}</td>'
            f'<td class="num">{thousands(hit.mismatches)}</td></tr>'
        )
    parts.append("</table>")
    if limit and len(hits) > limit:
        parts.append(
            f'<p class="note">Showing the first {plural(limit, "hit")} of '
            f"{thousands(len(hits))}. The rest are in the hits file, which "
            "--output writes in full.</p>"
        )
    return parts


def _record_section(
    record: RecordScan,
    scan: Scan,
    sequence: str,
    index: int,
    order: Mapping[str, int],
    width: int,
    highlight_limit: int,
    limit: int,
    map_width: int,
) -> tuple[list[str], list[str]]:
    """One record: its metrics, charts, sequence map and hits."""
    prefix = f"r{index}"
    unit = "aa" if record.alphabet == "protein" else "bp"
    parts = [
        "<details open>",
        f"<summary><h2>{escape(record.identifier)}</h2></summary>",
    ]
    if record.description:
        parts.append(f'<p class="lede">{escape(record.description)}</p>')
    rows = [
        ("Length", f"{thousands(record.residues)} {unit}"),
        ("Alphabet", record.alphabet),
        ("Topology", record.topology),
        ("Hits", thousands(len(record.hits))),
    ]
    if record.ambiguous:
        rows.append(("Ambiguous residues", thousands(record.ambiguous)))
    if record.source:
        rows.append(("File", record.source))
    parts.extend(_definitions(rows))

    applied = [compiled.name for compiled in scan.motifs]
    parts.extend(
        _summary_table(
            record.counts,
            scan.motifs,
            unit,
            scan.overlapping,
            ambiguity=bool(record.ambiguous) and not scan.include_uncertain,
        )
    )

    uncertain = [count for count in record.counts if count.uncertain]
    if uncertain:
        listed = ", ".join(
            f"{count.motif} ({plural(count.uncertain, 'position')})"
            for count in uncertain
        )
        fate = (
            "They are listed among the hits below, marked uncertain, and are "
            "counted neither as exact nor as mismatched: a position resting on "
            "an ambiguity code did not match and did not fail to match. The "
            "expectation is taken over every window, as the hits now are."
            if scan.include_uncertain
            else "They are not counted as hits, and the windows holding them "
            "are left out of the expectation as well, so that both figures are "
            "over the windows that could be called."
        )
        parts.append(
            f'<p class="note">Positions that could not be called either way, '
            f"because the sequence carries an ambiguity code there: "
            f"{escape(listed)}. {fate}</p>"
        )
    unassessed = [count.motif for count in record.counts if not count.assessed]
    if unassessed and record.ambiguous:
        parts.append(
            f'<p class="note">Uncertain positions were not assessed for '
            f"{escape(', '.join(unassessed))}: a pattern of variable length has no "
            "fixed positions to assess them over. A zero there means nothing was "
            "looked for, not that nothing was found.</p>"
        )

    chart = _distribution_chart(record.distribution, applied, unit)
    if chart:
        parts.append("<h3>Distribution</h3>")
        parts.extend(chart)

    scaled = [
        f"  {truncate(compiled.name, 16):<16} "
        f"{track(record.hits_for(compiled.name), record.residues, map_width)}"
        f"  {thousands(len(record.hits_for(compiled.name)))}"
        for compiled in scan.motifs
    ]
    if scaled and record.residues:
        scaled.extend(
            f"{' ' * 19}{line}"
            for line in ruler(record.residues, map_width, record.circular)
        )
        parts.append("<h3>Map</h3>")
        parts.append(f'<pre class="map">{escape(chr(10).join(scaled))}</pre>')

    rules: list[str] = []
    if record.hits:
        parts.append("<h3>Sequence</h3>")
        drawn = record.residues <= highlight_limit and bool(sequence)
        if drawn:
            found = [
                compiled.name
                for compiled in scan.motifs
                if record.hits_for(compiled.name)
            ]
            parts.extend(_key(found, order))
            parts.extend(_sequence_map(sequence, record.hits, order, prefix, width))
        else:
            parts.append(
                f'<p class="note">The residues are not shown: this record is '
                f"{thousands(record.residues)} {unit} long, past the "
                f"{thousands(highlight_limit)} {unit} this report highlights in "
                "full. Raise the limit to see them, or scan a region.</p>"
            )

        parts.append("<h3>Hits</h3>")
        strands = [
            label
            for label in ("+", "-", "=", ".")
            if any(hit.strand == label for hit in record.hits)
        ]
        controls, control_rules = _controls(
            prefix,
            applied,
            order,
            strands,
            uncertain=any(not hit.certain for hit in record.hits),
        )
        rules.extend(control_rules)
        parts.append('<section class="hits">')
        parts.extend(controls)
        parts.extend(_hits_table(record.hits, order, prefix, width, limit, drawn))
        parts.append("</section>")
        linked_note = (
            "A start position links to the line of the sequence above that holds "
            "it. "
            if drawn
            else ""
        )
        parts.append(
            f'<p class="note">{linked_note}The counts above the table are the '
            "counts of the whole scan: hiding rows here changes what is listed, "
            "not what was found.</p>"
        )
    else:
        parts.append("<p>No motif was found in this record.</p>")
    parts.append("</details>")
    return parts, rules


def render_html(
    scan: Scan,
    sequences: Sequence[str],
    version: str,
    *,
    width: int = 60,
    highlight_limit: int = 20000,
    limit: int = 500,
    map_width: int = 48,
) -> str:
    """Render a complete, self-contained HTML document for one scan.

    Args:
        scan: the scan to report.
        sequences: the residues of each record, in the same order as
            ``scan.records``, for the sequence map. Pass an empty sequence to
            report without the residues. Addressing them by position rather
            than by name is deliberate: two records of a FASTA file may carry
            the same identifier, and keying by it drew one record's hits over
            another record's residues.
        version: the version of the tool, for the footer.
        width: residues per line of the sequence map.
        highlight_limit: longest record whose residues are shown in full.
        limit: most hits listed per record; 0 lists every one.
        map_width: width of the scaled map, in characters.

    Raises:
        ValueError: if the residues given do not correspond one to one with the
            records of the scan.
    """
    if sequences and len(sequences) != len(scan.records):
        raise ValueError(
            f"{len(sequences)} sequences were given for "
            f"{len(scan.records)} records; they are matched by position"
        )
    order = _motif_index(scan)
    unit = "aa" if scan.alphabet == "protein" else "bp"
    sections: list[str] = []
    rules: list[str] = []
    for index, record in enumerate(scan.records):
        part, rule = _record_section(
            record,
            scan,
            sequences[index] if index < len(sequences) else "",
            index,
            order,
            width,
            highlight_limit,
            limit,
            map_width,
        )
        sections.extend(part)
        rules.extend(rule)

    parts = [
        "<!DOCTYPE html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>Motif scan of {escape(scan.name)}</title>",
        f"<style>{_STYLE}{chr(10).join(rules)}{chr(10)}</style>",
        "</head>",
        "<body>",
        f"<h1>Motif scan of {escape(scan.name)}</h1>",
        f'<p class="lede">{escape(plural(len(scan.motifs), "motif"))} &middot; '
        f'{escape(plural(len(scan.records), "record"))} &middot; '
        f"{thousands(scan.residues)} {escape(unit)} &middot; "
        f'{escape(plural(scan.hits, "hit"))}</p>',
    ]
    if scan.parameters:
        parts.append("<h2>Parameters</h2>")
        parts.extend(_definitions(list(scan.parameters)))

    if len(scan.records) > 1:
        parts.append("<h2>Motifs</h2>")
        parts.extend(
            _summary_table(
                scan.totals(),
                scan.motifs,
                unit,
                scan.overlapping,
                ambiguity=any(record.ambiguous for record in scan.records)
                and not scan.include_uncertain,
            )
        )

    parts.extend(sections)

    catalogued = sorted(
        {
            compiled.motif.source
            for compiled in scan.motifs
            if compiled.motif.source
        }
    )
    if catalogued:
        parts.append("<h2>Sources</h2>")
        parts.append('<ul class="note">')
        for compiled in scan.motifs:
            if compiled.motif.source:
                parts.append(
                    f"<li>{escape(compiled.name)}: {escape(compiled.motif.description)}"
                    f" &mdash; {escape(compiled.motif.source)}</li>"
                )
        parts.append("</ul>")
        parts.append(
            '<p class="note">A consensus is a blunt instrument: it over-reports, '
            "because a short consensus occurs by chance, and under-reports, "
            "because a real site that differs from the consensus at one position "
            "is still a real site. Read each count against its expectation, and "
            "for transcription factor binding prefer a matrix model over any "
            "consensus string.</p>"
        )

    parts.extend(
        [
            "<footer>",
            "<p>Positions are 1-based and inclusive, and always given on the "
            "forward strand; the matched text is read on the strand named. A hit "
            "whose end is smaller than its start spans the join of a circular "
            "sequence.</p>",
            f"<p>motif-scan {escape(version)}.</p>",
            "</footer>",
            "</body>",
            "</html>",
            "",
        ]
    )
    return "\n".join(parts)
