"""Residue alphabets: IUPAC nucleotide codes and amino-acid ambiguity codes.

A motif is written in one alphabet and searched in a sequence written in the
same one, and the two meet under a deliberately asymmetric rule set out in
:func:`matches`: an ambiguous code in the *motif* widens what the motif
accepts, while an ambiguous code in the *sequence* means the residue is not
known and cannot be reported as a match.

The two alphabets overlap completely. Every IUPAC nucleotide code is also an
amino-acid letter or ambiguity code - T is threonine, R arginine, N asparagine,
V valine - so no sequence can be assigned to an alphabet by looking for a
letter unique to it. Which alphabet a file is written in is therefore decided
by composition, in :func:`detect`, and the decision is always reported rather
than made silently.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

#: The bases each IUPAC nucleotide code stands for. U is uracil, and stands for
#: the same base as T: an RNA sequence is not an uncertain one.
NUCLEOTIDE_CODES: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "U": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

#: The twenty standard amino acids, plus the two encoded by recoded stop
#: codons, the three IUPAC-IUBMB ambiguity codes and X for an unknown residue.
#: ``*`` is admitted because a translated sequence carries stop codons, and a
#: motif may need to say where they are; it stands only for itself.
PROTEIN_CODES: dict[str, frozenset[str]] = {
    **{letter: frozenset(letter) for letter in "ACDEFGHIKLMNPQRSTVWY"},
    "U": frozenset("U"),  # selenocysteine
    "O": frozenset("O"),  # pyrrolysine
    "B": frozenset("DN"),
    "Z": frozenset("EQ"),
    "J": frozenset("IL"),
    "X": frozenset("ACDEFGHIKLMNPQRSTVWY"),
    "*": frozenset("*"),
}

#: Complement of every nucleotide code, including the ambiguity codes.
_COMPLEMENT = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn",
)


@dataclass(frozen=True)
class Alphabet:
    """One residue alphabet.

    Attributes:
        name: ``nucleotide`` or ``protein``.
        codes: the residues each code stands for.
        unit: what one residue is called, for counting in messages.
        has_strands: whether a sequence in this alphabet has a second strand to
            search. Protein does not, and asking for one is an error rather
            than a silently ignored flag.
    """

    name: str
    codes: Mapping[str, frozenset[str]]
    unit: str
    has_strands: bool
    wildcard: str

    @property
    def anything(self) -> frozenset[str]:
        """The residues a pattern wildcard stands for.

        Protein wildcards stand for the twenty standard amino acids, as they do
        in PROSITE, and so match neither a stop codon nor one of the two
        residues encoded by recoded stop codons. A motif that should match
        those names them.
        """
        return self.expand(self.wildcard)

    @property
    def definite(self) -> frozenset[str]:
        """Codes standing for exactly one residue."""
        return frozenset(
            code for code, residues in self.codes.items() if len(residues) == 1
        )

    @property
    def residues(self) -> frozenset[str]:
        """Every residue any code of this alphabet can stand for."""
        return frozenset().union(*self.codes.values())

    def expand(self, code: str) -> frozenset[str]:
        """The residues a code stands for.

        Raises:
            KeyError: if the character is not a code of this alphabet.
        """
        try:
            return self.codes[code.upper()]
        except KeyError:
            raise KeyError(f"not a {self.name} code: {code!r}") from None


#: Nucleotide sequence, which has two strands.
NUCLEOTIDE = Alphabet("nucleotide", NUCLEOTIDE_CODES, "bp", True, "N")

#: Protein sequence, which has one.
PROTEIN = Alphabet("protein", PROTEIN_CODES, "aa", False, "X")

#: Both alphabets by name.
ALPHABETS: dict[str, Alphabet] = {
    NUCLEOTIDE.name: NUCLEOTIDE,
    PROTEIN.name: PROTEIN,
}

#: Deleting every code that stands for a single residue leaves only the
#: ambiguous ones, which makes :func:`has_ambiguity` a C-speed scan rather than
#: a Python loop. The difference matters: the check decides whether the slower
#: of two search paths is needed at all, and it is asked once per sequence.
_SIEVES: dict[str, dict[int, int | None]] = {
    name: str.maketrans(
        "",
        "",
        "".join(sorted(alphabet.definite)) + "".join(sorted(alphabet.definite)).lower(),
    )
    for name, alphabet in ALPHABETS.items()
}

#: Residues that settle the alphabet of a sequence on sight, because they are
#: amino acids with no reading as a nucleotide code. X is included: a file
#: writing unknown *bases* as X rather than N is rejected with a message saying
#: so, which is better than scanning it as a peptide.
PROTEIN_ONLY: frozenset[str] = frozenset(PROTEIN_CODES) - frozenset(NUCLEOTIDE_CODES)

#: Codes counted as nucleotide when judging composition. The ambiguity codes
#: other than N are left out on purpose: R, Y, S and W are common amino acids
#: and rare bases, so counting them as nucleotide would pull peptides across
#: the line.
_NUCLEOTIDE_CORE: frozenset[str] = frozenset("ACGTUN")

#: Share of residues that must read as unambiguous nucleotide for a sequence to
#: be called nucleotide. The same 85% rule the EMBOSS and Biopython guessers
#: use; a coding sequence is ~100%, a peptide of this composition would have to
#: be made almost entirely of alanine, cysteine, glycine, threonine and
#: asparagine.
NUCLEOTIDE_SHARE = 0.85


def has_ambiguity(sequence: str, alphabet: Alphabet) -> bool:
    """Whether any residue stands for more than one residue.

    Used to decide whether uncertain positions need to be looked for at all, so
    it must never report False for a sequence that does contain an ambiguity
    code: that would silently switch off the uncertainty reporting rather than
    change its speed.
    """
    return bool(sequence.translate(_SIEVES[alphabet.name]))


def complement(sequence: str) -> str:
    """Complement a nucleotide sequence, preserving case and ambiguity codes."""
    return sequence.translate(_COMPLEMENT)


def reverse_complement(sequence: str) -> str:
    """Reverse complement a nucleotide sequence, preserving case and codes."""
    return sequence.translate(_COMPLEMENT)[::-1]


def is_palindrome(pattern: str) -> bool:
    """Whether a nucleotide pattern reads the same on both strands."""
    return bool(pattern) and pattern.upper() == reverse_complement(pattern.upper())


def matches(motif_code: str, sequence_code: str, alphabet: Alphabet) -> bool:
    """Whether a sequence residue certainly satisfies a motif position.

    The test is containment, not overlap: the residues the sequence might be
    must all be residues the motif accepts. So ``N`` in a nucleotide motif is
    satisfied by any base, while ``N`` in the sequence satisfies only ``N``.

    The asymmetry is the point. A motif written with an ambiguous code really
    does accept every residue that code covers, so widening there is correct.
    An ambiguous code in the sequence is missing information, and reporting a
    hit that may not be there is worse than reporting none: it sends someone
    to redesign a cassette around an element that was never found. Those
    positions are counted and reported separately as uncertain rather than
    silently dropped.
    """
    return alphabet.expand(sequence_code) <= alphabet.expand(motif_code)


def could_match(motif_code: str, sequence_code: str, alphabet: Alphabet) -> bool:
    """Whether a sequence residue might satisfy a motif position.

    True when the two codes have any residue in common, so an ``N`` in a
    nucleotide sequence might match anything. Used to find the positions that
    cannot be called either way.
    """
    return bool(alphabet.expand(motif_code) & alphabet.expand(sequence_code))


def detect(sequence: str) -> str | None:
    """Judge which alphabet a sequence is written in.

    Returns ``nucleotide``, ``protein``, or None when there is nothing to
    judge, which happens only for an empty sequence. A single amino acid with
    no nucleotide reading settles it outright; otherwise the decision is the
    share of residues that are unambiguous bases, against
    :data:`NUCLEOTIDE_SHARE`.

    The caller is expected to report the answer, not act on it silently: a
    twelve-residue peptide of glycines and alanines reads as DNA by any
    composition rule, and the only cure for that is for the user to see the
    guess and override it.
    """
    residues = sequence.upper()
    if not residues:
        return None
    if any(residue in PROTEIN_ONLY for residue in residues):
        return PROTEIN.name
    core = sum(1 for residue in residues if residue in _NUCLEOTIDE_CORE)
    if core / len(residues) >= NUCLEOTIDE_SHARE:
        return NUCLEOTIDE.name
    return PROTEIN.name


def get(name: str) -> Alphabet:
    """Look an alphabet up by name.

    Raises:
        KeyError: if the name is not an alphabet.
    """
    try:
        return ALPHABETS[name]
    except KeyError:
        known = ", ".join(sorted(ALPHABETS))
        raise KeyError(f"not an alphabet: {name!r} (known: {known})") from None
