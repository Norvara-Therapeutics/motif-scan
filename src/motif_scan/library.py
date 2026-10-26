"""The built-in motif catalogue.

A small, traceable set of elements that come up when designing an expression
cassette or reading a capsid sequence, each written as its published consensus
and each carrying the paper that consensus comes from. Nothing here is derived
from a position weight matrix, fitted to anything, or scored: a consensus is a
blunt instrument, and the catalogue says which blunt instrument it is.

That is also its limit. A consensus string over-reports - ``CANNTG`` occurs
every few hundred bases by chance - and under-reports, because a real binding
site that differs from the consensus at one position is a real binding site.
For the first, read the expected count the scan reports beside the observed
one; for the second, allow a mismatch. For transcription factor binding
specifically, a matrix model from JASPAR or HOCOMOCO is a better instrument
than anything in this file, and this catalogue is not a substitute for one.

Pattern syntax is PROSITE for every protein motif and plain IUPAC codes for
every nucleotide one, except where a nucleotide motif needs a repeat count.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import IO, Iterable, Mapping, Sequence

from motif_scan.patterns import SYNTAXES, Motif, PatternError

#: The catalogue. Kept as a flat tuple in a deliberate order - core promoter,
#: then the factors, then processing, translation, vector elements, and the
#: protein motifs - so that output reads in the order a cassette is built.
CATALOGUE: tuple[Motif, ...] = (
    Motif(
        name="TATA",
        pattern="TATAWAWR",
        category="core-promoter",
        description="TATA box, about 25-30 bp upstream of the start site",
        source="Bucher 1990, J Mol Biol 212:563",
    ),
    Motif(
        name="BRE",
        pattern="SSRCGCC",
        category="core-promoter",
        description="TFIIB recognition element, immediately upstream of the TATA box",
        source="Lagrange et al. 1998, Genes Dev 12:34",
    ),
    Motif(
        name="Inr",
        pattern="YYANWYY",
        category="core-promoter",
        description="initiator element, overlapping the transcription start site",
        source="Javahery et al. 1994, Mol Cell Biol 14:116",
    ),
    Motif(
        name="DPE",
        pattern="RGWYVT",
        category="core-promoter",
        description="downstream promoter element, about +28 to +33",
        source="Kutach and Kadonaga 2000, Mol Cell Biol 20:4754",
    ),
    Motif(
        name="CCAAT",
        pattern="CCAAT",
        category="core-promoter",
        description="CCAAT box, bound by NF-Y",
        source="Mantovani 1998, Nucleic Acids Res 26:1135",
    ),
    Motif(
        name="GC-box",
        pattern="GGGGCGGGGC",
        category="core-promoter",
        description="GC box, bound by Sp1",
        source="Briggs et al. 1986, Science 234:47",
    ),
    Motif(
        name="CACCC",
        pattern="CACCC",
        category="core-promoter",
        description="CACCC box, bound by the KLF family; essential in the "
        "beta-globin promoter at -90",
        source="Dierks et al. 1983, Cell 32:695",
    ),
    Motif(
        name="AP-1",
        pattern="TGASTCA",
        category="transcription-factor",
        description="AP-1 site (TRE), bound by Jun/Fos dimers",
        source="Angel et al. 1987, Cell 49:729",
    ),
    Motif(
        name="CRE",
        pattern="TGACGTCA",
        category="transcription-factor",
        description="cyclic AMP response element, bound by CREB",
        source="Montminy et al. 1986, Proc Natl Acad Sci USA 83:6682",
    ),
    Motif(
        name="NF-kB",
        pattern="GGGRNNYYCC",
        category="transcription-factor",
        description="NF-kappaB site",
        source="Kunsch et al. 1992, Mol Cell Biol 12:4412",
    ),
    Motif(
        name="E-box",
        pattern="CANNTG",
        category="transcription-factor",
        description="E box, bound by basic helix-loop-helix factors",
        source="Ephrussi et al. 1985, Science 227:134",
    ),
    Motif(
        name="GATA",
        pattern="WGATAR",
        category="transcription-factor",
        description="GATA factor site",
        source="Merika and Orkin 1993, Mol Cell Biol 13:3999",
    ),
    Motif(
        name="octamer",
        pattern="ATGCAAAT",
        category="transcription-factor",
        description="octamer element, bound by POU factors such as Oct-1",
        source="Falkner and Zachau 1984, Nature 310:71",
    ),
    Motif(
        name="HSE",
        pattern="TTCNNGAANNTTC",
        category="transcription-factor",
        description="heat shock element: three inverted nGAAn repeats",
        source="Amin et al. 1988, Mol Cell Biol 8:3761",
    ),
    Motif(
        name="CArG",
        pattern="CCWWWWWWGG",
        category="transcription-factor",
        description="CArG box of the serum response element",
        source="Minty and Kedes 1986, Mol Cell Biol 6:2125",
    ),
    Motif(
        name="GRE",
        pattern="GGTACANNNTGTTCT",
        category="transcription-factor",
        description="glucocorticoid response element, an imperfect palindrome",
        source="Beato 1989, Cell 56:335",
    ),
    Motif(
        name="polyA",
        pattern="AATAAA",
        category="processing",
        description="polyadenylation signal",
        source="Proudfoot and Brownlee 1976, Nature 263:211",
    ),
    Motif(
        name="polyA-variant",
        pattern="ATTAAA",
        category="processing",
        description="the commonest single-base variant of the polyadenylation signal",
        source="Beaudoing et al. 2000, Genome Res 10:1001",
    ),
    Motif(
        name="splice-donor",
        pattern="MAGGTRAGT",
        category="processing",
        description="5' splice site; the intron begins at the GT",
        source="Shapiro and Senapathy 1987, Nucleic Acids Res 15:7155",
    ),
    Motif(
        name="splice-acceptor",
        pattern="YYYYYYYYYYNCAGG",
        category="processing",
        description="3' splice site: polypyrimidine tract, then the AG the intron ends at",
        source="Shapiro and Senapathy 1987, Nucleic Acids Res 15:7155",
    ),
    Motif(
        name="branch-point",
        pattern="YTRAY",
        category="processing",
        description="intron branch point, 20-50 bases upstream of the 3' splice site",
        source="Gao et al. 2008, Nucleic Acids Res 36:2257",
    ),
    Motif(
        name="Kozak",
        pattern="GCCRCCATGG",
        category="translation",
        description="Kozak consensus around the initiator ATG",
        source="Kozak 1987, Nucleic Acids Res 15:8125",
    ),
    Motif(
        name="Kozak-core",
        pattern="RNNATGG",
        category="translation",
        description="the two context positions that matter most: a purine at -3 "
        "and G at +4",
        source="Kozak 1986, Cell 44:283",
    ),
    Motif(
        name="CpG",
        pattern="CG",
        category="composition",
        description="CpG dinucleotide; unmethylated CpG in a vector genome is a "
        "TLR9 ligand, and depleting it lowers the response to an AAV payload",
        source="Faust et al. 2013, J Clin Invest 123:2994",
    ),
    Motif(
        name="AAV-RBE",
        pattern="GCTCGCTC",
        category="vector",
        description="Rep binding element: two of the four imperfect GCTC repeats "
        "of the AAV inverted terminal repeat",
        source="McCarty et al. 1994, J Virol 68:4988",
    ),
    Motif(
        name="AAV-trs",
        pattern="GTTGGCC",
        category="vector",
        description="terminal resolution site; Rep nicks between the two "
        "thymidines, and only this strand is required",
        source="Brister and Muzyczka 1999, J Virol 73:9325",
    ),
    Motif(
        name="N-glycosylation",
        pattern="N-{P}-[ST]-{P}",
        syntax="prosite",
        alphabet="protein",
        category="protein-modification",
        description="N-glycosylation sequon; occupancy also depends on structure",
        source="PROSITE PS00001; Gavel and von Heijne 1990, Protein Eng 3:433",
    ),
    Motif(
        name="CK2-site",
        pattern="[ST]-x(2)-[DE]",
        syntax="prosite",
        alphabet="protein",
        category="protein-modification",
        description="casein kinase II phosphorylation site",
        source="PROSITE PS00006",
    ),
    Motif(
        name="PKA-site",
        pattern="[RK](2)-x-[ST]",
        syntax="prosite",
        alphabet="protein",
        category="protein-modification",
        description="cAMP-dependent protein kinase phosphorylation site",
        source="PROSITE PS00004",
    ),
    Motif(
        name="myristoylation",
        pattern="G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}",
        syntax="prosite",
        alphabet="protein",
        category="protein-modification",
        description="N-myristoylation site at a glycine exposed by cleavage of "
        "the initiator methionine",
        source="PROSITE PS00008",
    ),
    Motif(
        name="CAAX",
        pattern="C-[AGILMV](2)-x>",
        syntax="prosite",
        alphabet="protein",
        category="protein-modification",
        description="C-terminal prenylation box",
        source="Zhang and Casey 1996, Annu Rev Biochem 65:241",
    ),
    Motif(
        name="RGD",
        pattern="R-G-D",
        syntax="prosite",
        alphabet="protein",
        category="protein-interaction",
        description="integrin-binding cell attachment sequence",
        source="PROSITE PS00016; Ruoslahti and Pierschbacher 1986, Cell 44:517",
    ),
    Motif(
        name="heparin-binding",
        pattern="x-[RK](2)-x-[RK]-x",
        syntax="prosite",
        alphabet="protein",
        category="protein-interaction",
        description="Cardin-Weintraub XBBXBX heparin binding motif",
        source="Cardin and Weintraub 1989, Arteriosclerosis 9:21",
    ),
    Motif(
        name="PLA2",
        pattern="H-D-x(2)-Y",
        syntax="prosite",
        alphabet="protein",
        category="protein-interaction",
        description="catalytic motif of the phospholipase A2 domain in the "
        "unique region of AAV VP1, required for endosomal escape",
        source="Zadori et al. 2001, Dev Cell 1:291",
    ),
    Motif(
        name="NLS",
        pattern="K-[KR]-x-[KR]",
        syntax="prosite",
        alphabet="protein",
        category="protein-trafficking",
        description="classical monopartite nuclear localisation signal",
        source="Chelsky et al. 1989, Mol Cell Biol 9:2487",
    ),
    Motif(
        name="NLS-bipartite",
        pattern="[KR](2)-x(10,12)-[KR](3,5)",
        syntax="prosite",
        alphabet="protein",
        category="protein-trafficking",
        description="bipartite nuclear localisation signal; written strictly, "
        "as a run of basic residues rather than three basics within five",
        source="Dingwall and Laskey 1991, Trends Biochem Sci 16:478",
    ),
    Motif(
        name="KDEL",
        pattern="[KR]-D-E-L>",
        syntax="prosite",
        alphabet="protein",
        category="protein-trafficking",
        description="C-terminal endoplasmic reticulum retention signal",
        source="Munro and Pelham 1987, Cell 48:899",
    ),
    Motif(
        name="furin-site",
        pattern="R-x-[KR]-R",
        syntax="prosite",
        alphabet="protein",
        category="protein-cleavage",
        description="furin cleavage site; cleavage follows the last arginine",
        source="Hosaka et al. 1991, J Biol Chem 266:12127",
    ),
    Motif(
        name="TEV-site",
        pattern="E-x(2)-Y-x-Q-[GS]",
        syntax="prosite",
        alphabet="protein",
        category="protein-cleavage",
        description="tobacco etch virus protease site; cleavage follows the glutamine",
        source="Dougherty et al. 1989, Virology 172:302",
    ),
    Motif(
        name="His-tag",
        pattern="H(6)",
        syntax="prosite",
        alphabet="protein",
        category="protein-tag",
        description="hexahistidine affinity tag",
        source="Hochuli et al. 1988, Nat Biotechnol 6:1321",
    ),
)

#: Expected columns of a motif file, in order. ``name`` and ``pattern`` are
#: required; the rest may be left out or left empty.
COLUMNS: tuple[str, ...] = (
    "name",
    "pattern",
    "syntax",
    "alphabet",
    "category",
    "description",
    "source",
)

_BY_NAME: dict[str, Motif] = {motif.name.lower(): motif for motif in CATALOGUE}


def all_motifs() -> tuple[Motif, ...]:
    """Every motif in the catalogue."""
    return CATALOGUE


def categories() -> tuple[str, ...]:
    """Every category present, in catalogue order."""
    seen: dict[str, None] = {}
    for motif in CATALOGUE:
        seen.setdefault(motif.category, None)
    return tuple(seen)


def get(name: str) -> Motif:
    """Look one motif up by name, ignoring case.

    Raises:
        KeyError: if no motif has that name, with the closest names offered.
    """
    try:
        return _BY_NAME[name.strip().lower()]
    except KeyError:
        near = [
            motif.name
            for motif in CATALOGUE
            if name.strip().lower() in motif.name.lower()
        ]
        hint = f"; did you mean {', '.join(near)}?" if near else ""
        raise KeyError(f"not a motif in the catalogue: {name!r}{hint}") from None


def resolve(names: Iterable[str]) -> list[Motif]:
    """Look up several motifs by name, keeping the order asked for.

    A name repeated is returned once: applying a motif twice would double every
    count it contributes to.
    """
    chosen: dict[str, Motif] = {}
    for name in names:
        motif = get(name)
        chosen.setdefault(motif.name, motif)
    return list(chosen.values())


def select(
    category: str = "",
    alphabet: str = "",
    search: str = "",
) -> list[Motif]:
    """Every catalogue motif matching all of the filters given.

    Args:
        category: keep only this category.
        alphabet: keep only motifs written in this alphabet.
        search: keep only motifs whose name, description or pattern contains
            this text, ignoring case.
    """
    needle = search.strip().lower()
    return [
        motif
        for motif in CATALOGUE
        if (not category or motif.category == category)
        and (not alphabet or motif.alphabet == alphabet)
        and (
            not needle
            or needle in motif.name.lower()
            or needle in motif.description.lower()
            or needle in motif.pattern.lower()
        )
    ]


def _row_motif(row: Mapping[str, str], where: str) -> Motif:
    """Build one motif from a row of a motif file.

    The columns themselves are checked once, against the header, by
    :func:`read_motifs`; this reads a row that has already been matched to it.
    """
    fields = {key: (value or "").strip() for key, value in row.items()}
    name = fields.get("name", "")
    pattern = fields.get("pattern", "")
    if not name or not pattern:
        raise PatternError(f"{where}: both a name and a pattern are required")
    syntax = fields.get("syntax") or "iupac"
    if syntax not in SYNTAXES:
        raise PatternError(
            f"{where}: not a pattern syntax: {syntax!r} "
            f"(known: {', '.join(SYNTAXES)})"
        )
    alphabet = fields.get("alphabet") or "nucleotide"
    if alphabet not in ("nucleotide", "protein"):
        raise PatternError(
            f"{where}: not an alphabet: {alphabet!r} (known: nucleotide, protein)"
        )
    return Motif(
        name=name,
        pattern=pattern,
        syntax=syntax,
        alphabet=alphabet,
        category=fields.get("category", "file"),
        description=fields.get("description", ""),
        source=fields.get("source", ""),
    )


def read_motifs(handle: IO[str], where: str = "motif file") -> list[Motif]:
    """Read motifs from an open tab-separated file.

    The first non-comment line is the header, and must name its columns from
    :data:`COLUMNS`; ``name`` and ``pattern`` are the two that are required.
    Lines beginning with ``#`` and empty lines are ignored, so a file can carry
    its own provenance.

    The file has to be rectangular, and its columns have to be distinct. A
    dictionary reader would accept neither problem noisily: a repeated column
    name silently keeps the last of the two, so a header of
    ``name pattern name`` quietly renames every motif in the file, and a row
    with a field too many has the extra one filed under no column at all and
    dropped. Both are read here as what they are - a malformed file - because
    the alternative is a scan that runs, succeeds, and reports a motif nobody
    wrote.

    Raises:
        PatternError: if the file has no header, a repeated or unknown column,
            a row whose length does not match the header, or a row that cannot
            be read as a motif.
    """
    lines = [
        line
        for line in handle
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not lines:
        raise PatternError(f"{where}: holds no motifs")
    rows = list(csv.reader(lines, delimiter="\t"))
    header = [name.strip().lower() for name in rows[0]]
    if not header or not any(header):
        raise PatternError(f"{where}: has no header line")

    repeated = sorted({name for name in header if header.count(name) > 1})
    if repeated:
        shown = ", ".join(repr(name) for name in repeated)
        raise PatternError(
            f"{where}: the header names a column more than once: {shown}; "
            "which of the two a row meant cannot be known"
        )
    unknown = sorted(set(header) - set(COLUMNS))
    if unknown:
        shown = ", ".join(repr(column) for column in unknown)
        raise PatternError(
            f"{where}: unknown column(s): {shown} (expected {', '.join(COLUMNS)})"
        )
    if "name" not in header or "pattern" not in header:
        raise PatternError(
            f"{where}: the header must name at least the columns 'name' and "
            f"'pattern'; found {', '.join(rows[0])}"
        )

    motifs: list[Motif] = []
    for number, row in enumerate(rows[1:], start=2):
        if len(row) != len(header):
            raise PatternError(
                f"{where} line {number}: has {len(row)} field(s) against "
                f"{len(header)} column(s) in the header"
            )
        motifs.append(_row_motif(dict(zip(header, row)), f"{where} line {number}"))
    if not motifs:
        raise PatternError(f"{where}: header only, no motifs")
    return motifs


def load(path: str | Path) -> list[Motif]:
    """Read motifs from a tab-separated file.

    Raises:
        PatternError: if the file cannot be read as motifs.
        OSError: if it cannot be opened.
    """
    where = "-" if str(path) == "-" else str(path)
    if where == "-":
        import sys

        return read_motifs(sys.stdin, "standard input")
    with Path(path).open(encoding="utf-8-sig") as handle:
        return read_motifs(handle, where)


def unique(motifs: Sequence[Motif]) -> list[Motif]:
    """Reject a set of motifs with a repeated name.

    Two motifs of the same name would be indistinguishable in every table and
    in the hits file, and one would silently overwrite the other wherever
    results are kept by name.

    Raises:
        PatternError: naming the motifs that collide.
    """
    seen: dict[str, Motif] = {}
    clashes: list[str] = []
    for motif in motifs:
        if motif.name in seen and seen[motif.name] != motif:
            clashes.append(motif.name)
        seen.setdefault(motif.name, motif)
    if clashes:
        shown = ", ".join(sorted(set(clashes)))
        raise PatternError(
            f"more than one motif is named: {shown}; rename one with NAME=PATTERN"
        )
    return list(seen.values())
