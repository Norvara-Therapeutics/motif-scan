"""Motif patterns: what one can be written in, and what it compiles to.

Three notations are accepted, and they are not interchangeable.

``iupac``
    A string of alphabet codes, one per position: ``TATAWAWR``, ``CANNTG``.
    Fixed length, so it can be searched with mismatches allowed and
    reverse-complemented position by position.

``prosite``
    PROSITE pattern syntax: ``N-{P}-[ST]-{P}``, ``[ST]-x(2)-[DE]``,
    ``C-[AVILM](2)-x>``. Residue sets, exclusions, wildcards, repeat counts,
    and the ``<`` and ``>`` terminal anchors. A pattern whose repeats are all
    exact counts has a fixed length and behaves like an ``iupac`` one; a
    pattern with a range such as ``x(2,4)`` does not.

``regex``
    A Python regular expression, applied to the sequence as written. The most
    expressive and the least constrained: no mismatch tolerance, no length,
    and no guarantee the author of the pattern meant what the engine does. It
    is compiled exactly as given, with nothing wrapped around it, so its own
    groups keep their numbers and its own inline flags stay where they were
    written; the overlapping scan it needs is done by matching at each start
    position instead.

Everything that follows depends on the distinction between a pattern of fixed
length, which is held as one set of permitted residues per position, and one
of variable length, which is held only as a compiled regular expression. The
first can be searched with mismatches, reverse-complemented exactly, and have
its expected frequency computed; the second cannot, and says so rather than
quietly returning zero.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable, Sequence

from motif_scan.alphabet import Alphabet, get as get_alphabet

#: The notations a pattern can be written in.
SYNTAXES: tuple[str, ...] = ("iupac", "prosite", "regex")

#: Name of the capture group wrapped around a pattern this module builds
#: itself, which is every pattern except a user-written regular expression.
#: Those are compiled untouched: wrapping one in a group renumbers the groups
#: inside it, so ``([AC])\\1`` becomes a reference to the wrapper and fails to
#: compile, and prefixing anything at all moves an inline ``(?i)`` off the
#: start of the expression, where Python requires it to be.
GROUP = "motif_scan_hit"

#: Longest repeat count a PROSITE token may ask for. A pattern is a motif, not
#: a homopolymer run finder, and ``x(100000)`` compiles to a regular expression
#: large enough to matter.
MAX_REPEAT = 1000


class PatternError(ValueError):
    """Raised when a pattern cannot be read, or cannot do what is asked of it."""


@dataclass(frozen=True)
class Motif:
    """A named pattern, before it is compiled.

    Attributes:
        name: what the motif is called in output. Unique within one scan.
        pattern: the pattern as written, in the notation given by ``syntax``.
        syntax: one of :data:`SYNTAXES`.
        alphabet: ``nucleotide`` or ``protein``.
        description: one line on what the motif is.
        category: grouping used to filter the catalogue.
        source: where the consensus comes from, for anything from the
            catalogue. Empty for a pattern given on the command line.
    """

    name: str
    pattern: str
    syntax: str = "iupac"
    alphabet: str = "nucleotide"
    description: str = ""
    category: str = ""
    source: str = ""

    def __post_init__(self) -> None:
        if not self.name:
            raise PatternError("a motif must have a name")
        if self.syntax not in SYNTAXES:
            known = ", ".join(SYNTAXES)
            raise PatternError(
                f"{self.name}: not a pattern syntax: {self.syntax!r} (known: {known})"
            )
        if not self.pattern:
            raise PatternError(f"{self.name}: pattern is empty")


@dataclass(frozen=True)
class Compiled:
    """A motif prepared for searching.

    Attributes:
        motif: the motif this came from.
        alphabet: the alphabet it is searched in.
        positions: permitted residues at each position, for a fixed-length
            pattern. Empty for a pattern of variable length.
        regex: the pattern as a regular expression that admits overlapping
            matches, used for the exact-match search path.
        anchor_start: the motif must begin at the first residue.
        anchor_end: the motif must end at the last residue.
        minimum: shortest text the pattern can match.
        maximum: longest text the pattern can match, or None when the pattern
            sets no upper bound. A pattern with no bound cannot be searched
            across the join of a circular sequence, because there is no way to
            know how much of the start to append.
    """

    motif: Motif
    alphabet: Alphabet
    regex: re.Pattern[str]
    positions: tuple[frozenset[str], ...] = ()
    anchor_start: bool = False
    anchor_end: bool = False
    minimum: int = 0
    maximum: int | None = None
    _complemented: tuple[frozenset[str], ...] = field(default=(), repr=False)

    @property
    def name(self) -> str:
        """The motif's name."""
        return self.motif.name

    @property
    def fixed_length(self) -> int | None:
        """How long every match must be, or None for a variable-length pattern."""
        return len(self.positions) if self.positions else None

    @property
    def supports_mismatches(self) -> bool:
        """Whether mismatches can be allowed against this pattern.

        Only a fixed-length pattern can: a mismatch is a position that does not
        match, and a pattern of variable length does not have a fixed set of
        positions to count against. Asking for mismatches on a regular
        expression is refused rather than ignored.
        """
        return bool(self.positions)

    @property
    def palindromic(self) -> bool:
        """Whether the pattern reads identically on the reverse strand.

        Only meaningful for a fixed-length nucleotide pattern, where it is
        decided position by position over the permitted residue sets, so that
        ``GGWCC`` and ``CANNTG`` are recognised as palindromes even though
        neither is a literal one. A variable-length pattern is never claimed to
        be palindromic: whether a regular expression is its own reverse
        complement is not a question the engine can answer, so both strands are
        searched and a palindromic pattern will report each site twice.
        """
        return bool(self.positions) and self._complemented == self.positions

    @property
    def anchored(self) -> bool:
        """Whether either terminus of the sequence is required."""
        return self.anchor_start or self.anchor_end

    @property
    def wrapped(self) -> bool:
        """Whether the compiled pattern carries this module's capture group.

        True for everything built here from residue sets, false for a
        user-written regular expression, which is compiled as it was written.
        The scanner needs to know which group holds the matched text.
        """
        return self.motif.syntax != "regex"

    def matched(self, match: re.Match[str]) -> str:
        """The text one match covers."""
        return match.group(GROUP) if self.wrapped else match.group(0)


_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def complement_positions(
    positions: Sequence[frozenset[str]], alphabet: Alphabet
) -> tuple[frozenset[str], ...]:
    """Reverse-complement a positional pattern, set by set.

    Returns an empty tuple for an alphabet with one strand, so that a protein
    pattern has no reverse reading rather than a meaningless one.
    """
    if not alphabet.has_strands:
        return ()
    return tuple(
        frozenset(_COMPLEMENT[residue] for residue in allowed)
        for allowed in reversed(positions)
    )


def _character_class(allowed: Iterable[str], alphabet: Alphabet) -> str:
    """Every sequence code a set of permitted residues admits, as a regex class.

    The class is built by asking which codes the position accepts rather than
    by listing residues, because the sequence may itself be written in
    ambiguity codes: a position permitting every base must admit a sequence
    ``N`` explicitly, since ``N`` is not one of A, C, G or T. Deriving the
    class this way keeps the fast path and
    :func:`motif_scan.alphabet.matches` from drifting apart.
    """
    permitted = frozenset(allowed)
    codes = sorted(
        code for code, residues in alphabet.codes.items() if residues <= permitted
    )
    if not codes:
        # No code can satisfy this position, so nothing can match it. A class
        # that matches nothing is spelled as a negated class of everything.
        return "[^\\s\\S]"
    both_cases = "".join(codes) + "".join(codes).lower()
    return f"[{re.escape(both_cases)}]"


def _positional_regex(
    positions: Sequence[frozenset[str]], alphabet: Alphabet
) -> re.Pattern[str]:
    """Compile a fixed-length pattern, keeping overlapping matches."""
    body = "".join(_character_class(allowed, alphabet) for allowed in positions)
    return re.compile(f"(?=(?P<{GROUP}>{body}))")


def _wrap(body: str) -> re.Pattern[str]:
    """Wrap a regular expression so that overlapping matches are all found.

    ``finditer`` resumes after the end of each match, so a motif that tiles
    along a sequence - ``CANNTG`` in a CA-rich region, ``CG`` in an island -
    would report every other occurrence. A zero-width lookahead consumes
    nothing, so every start position is tested.
    """
    try:
        return re.compile(f"(?=(?P<{GROUP}>{body}))")
    except re.error as error:
        raise PatternError(f"not a valid regular expression: {error}") from None


def parse_iupac(pattern: str, alphabet: Alphabet) -> tuple[frozenset[str], ...]:
    """Read a string of alphabet codes as one permitted set per position."""
    text = pattern.strip().upper()
    if not text:
        raise PatternError("pattern is empty")
    unknown = sorted(set(text) - set(alphabet.codes))
    if unknown:
        shown = ", ".join(repr(code) for code in unknown)
        hint = ""
        if any(character in "[]{}()<>.-" for character in unknown):
            hint = "; bracket and repeat notation needs --prosite"
        raise PatternError(
            f"not {alphabet.name} codes: {shown}{hint}"
        )
    return tuple(alphabet.expand(code) for code in text)


_TOKEN = re.compile(
    r"""
    (?P<body> \[ (?P<any> [^\]]* ) \]      # [ST]   any of these
            | \{ (?P<not> [^}]*  ) \}      # {P}    any except these
            | (?P<wild> [xX] )             # x      any residue
            | (?P<one> [A-Za-z*] )         # A      this residue
    )
    (?: \( \s* (?P<low> \d+ ) \s*          # (2)    exactly twice
        (?: , \s* (?P<high> \d+ ) \s* )?   # (2,4)  two to four times
    \) )?
    """,
    re.VERBOSE,
)


def _token_set(match: re.Match[str], alphabet: Alphabet, token: str) -> frozenset[str]:
    """The residues one PROSITE token permits."""
    if match.group("wild"):
        return alphabet.anything
    if match.group("one"):
        code = match.group("one").upper()
        if code not in alphabet.codes:
            raise PatternError(f"{token!r}: not a {alphabet.name} code: {code!r}")
        return alphabet.expand(code)

    listed = (match.group("any") or match.group("not") or "").upper()
    if not listed:
        raise PatternError(f"{token!r}: empty residue set")
    unknown = sorted(set(listed) - set(alphabet.codes))
    if unknown:
        shown = ", ".join(repr(code) for code in unknown)
        raise PatternError(f"{token!r}: not {alphabet.name} codes: {shown}")
    residues = frozenset().union(*(alphabet.expand(code) for code in listed))
    if match.group("not") is not None:
        residues = alphabet.anything - residues
        if not residues:
            raise PatternError(f"{token!r}: excludes every residue, so cannot match")
    return residues


def parse_prosite(
    pattern: str, alphabet: Alphabet
) -> tuple[tuple[frozenset[str], ...], str, bool, bool, int]:
    """Read a PROSITE pattern.

    Returns the permitted set at each position - empty when the pattern has a
    repeat range and so no fixed length - the body of the equivalent regular
    expression, whether each terminus is anchored, and the longest text the
    pattern can match.

    Raises:
        PatternError: if the pattern cannot be read.
    """
    text = pattern.strip().rstrip(".")
    anchor_start = text.startswith("<")
    text = text[1:] if anchor_start else text
    anchor_end = text.endswith(">")
    text = text[:-1] if anchor_end else text
    if not text:
        raise PatternError("pattern is empty")

    tokens = [token.strip() for token in text.split("-")]
    if any(not token for token in tokens):
        raise PatternError(
            f"{pattern!r}: empty token; PROSITE tokens are joined by single hyphens"
        )

    positions: list[frozenset[str]] = []
    body: list[str] = []
    fixed = True
    longest = 0
    for token in tokens:
        match = _TOKEN.fullmatch(token)
        if match is None:
            raise PatternError(
                f"{token!r}: not a PROSITE token; expected a residue, [ABC], "
                "{ABC} or x, with an optional (n) or (n,m) repeat"
            )
        allowed = _token_set(match, alphabet, token)
        low = int(match.group("low")) if match.group("low") else 1
        high = int(match.group("high")) if match.group("high") else low
        if high < low:
            raise PatternError(f"{token!r}: repeat range runs backwards")
        if low < 1:
            raise PatternError(f"{token!r}: a repeat of zero cannot be searched for")
        if high > MAX_REPEAT:
            raise PatternError(
                f"{token!r}: repeat of {high} exceeds the limit of {MAX_REPEAT}"
            )
        klass = _character_class(allowed, alphabet)
        longest += high
        if low == high:
            positions.extend([allowed] * low)
            body.append(klass if low == 1 else f"{klass}{{{low}}}")
        else:
            fixed = False
            body.append(f"{klass}{{{low},{high}}}")

    return (
        tuple(positions) if fixed else (),
        "".join(body),
        anchor_start,
        anchor_end,
        longest,
    )


def compile_motif(motif: Motif, alphabet: Alphabet | None = None) -> Compiled:
    """Prepare a motif for searching.

    Args:
        motif: the motif to compile.
        alphabet: the alphabet to compile it against. Defaults to the one the
            motif declares, and must agree with it: a protein motif cannot be
            searched in a nucleotide sequence.

    Raises:
        PatternError: if the pattern cannot be read, if it is written in
            another alphabet, or if it can match an empty string.
    """
    declared = get_alphabet(motif.alphabet)
    target = alphabet or declared
    if target.name != declared.name:
        raise PatternError(
            f"{motif.name}: a {declared.name} motif cannot be searched in "
            f"{target.name} sequence"
        )

    positions: tuple[frozenset[str], ...] = ()
    anchor_start = anchor_end = False
    maximum: int | None = None
    try:
        if motif.syntax == "iupac":
            positions = parse_iupac(motif.pattern, target)
            regex = _positional_regex(positions, target)
            minimum = maximum = len(positions)
        elif motif.syntax == "prosite":
            positions, body, anchor_start, anchor_end, maximum = parse_prosite(
                motif.pattern, target
            )
            regex = _wrap(body)
            minimum = len(positions) if positions else _shortest(body)
        else:
            # Compiled untouched: see GROUP for why nothing may be wrapped
            # around a pattern somebody else wrote.
            try:
                regex = re.compile(motif.pattern)
            except re.error as error:
                raise PatternError(
                    f"not a valid regular expression: {error}"
                ) from None
            minimum = _shortest(motif.pattern)
    except PatternError as error:
        raise PatternError(f"{motif.name}: {error}") from None

    # A pattern that matches nothing at all would report a hit at every
    # position, so it is refused here rather than at every one of them.
    empty = regex.match("") if motif.syntax == "regex" else regex.search("")
    if empty is not None and not (
        empty.group(GROUP) if motif.syntax != "regex" else empty.group(0)
    ):
        raise PatternError(
            f"{motif.name}: pattern can match an empty string, so it would "
            "report a hit at every position"
        )

    return Compiled(
        motif=motif,
        alphabet=target,
        regex=regex,
        positions=positions,
        anchor_start=anchor_start,
        anchor_end=anchor_end,
        minimum=minimum,
        maximum=maximum,
        _complemented=complement_positions(positions, target),
    )


def _shortest(body: str) -> int:
    """A cheap lower bound on the length of text a regex body can match.

    Only used to skip sequences that are too short to hold the motif at all, so
    an underestimate is safe and zero is the safe default. Nothing is parsed:
    the bound comes from the smallest repeat counts written in the pattern,
    which is why a pattern with no repeat notation yields zero rather than its
    literal length.
    """
    counts = [int(low) for low in re.findall(r"\{\s*(\d+)", body)]
    return min(counts) if counts else 0


#: What may stand before the ``=`` of ``NAME=PATTERN``. A regular expression
#: may itself contain an equals sign - in a lookahead, or a named
#: backreference such as ``(?P=x)`` - so text that cannot be a name is taken
#: as the whole pattern rather than split in two.
_NAME = re.compile(r"[\w.\-]+")


def named(text: str, syntax: str, alphabet: str, index: int) -> Motif:
    """Read a motif given on the command line as ``NAME=PATTERN`` or ``PATTERN``.

    An unnamed pattern is named after itself, which keeps output readable and
    self-explanatory. A pattern long enough to crowd a table is named after its
    position on the command line instead.

    Raises:
        PatternError: if the name is empty or the pattern is missing.
    """
    name, separator, pattern = text.partition("=")
    if separator and name.strip() and not _NAME.fullmatch(name.strip()):
        # The equals sign belongs to the pattern, not to a name. An *empty*
        # name is still an error: nobody means a pattern beginning with '='.
        name, separator, pattern = "", "", text
    if not separator:
        pattern = text.strip()
        if not pattern:
            raise PatternError("pattern is empty")
        name = pattern if len(pattern) <= 24 else f"{syntax}-{index}"
    else:
        name = name.strip()
        pattern = pattern.strip()
        if not name:
            raise PatternError(f"{text!r}: name before '=' is empty")
        if not pattern:
            raise PatternError(f"{text!r}: pattern after '=' is empty")
    return Motif(
        name=name,
        pattern=pattern,
        syntax=syntax,
        alphabet=alphabet,
        description="given on the command line",
    )
