"""Pattern and regulatory element finding in nucleotide and protein sequence.

Positions are 1-based and inclusive, and always reported on the forward strand;
the matched text is reported as the motif reads it, so a hit on the reverse
strand shows the element rather than its complement. An ambiguity code in a
pattern widens what the pattern accepts, while an ambiguity code in the
sequence means the residue is not known: those positions are counted as
uncertain rather than reported as hits.
"""

from motif_scan.alphabet import (
    NUCLEOTIDE,
    PROTEIN,
    Alphabet,
    detect,
    reverse_complement,
)
from motif_scan.library import all_motifs, categories, load, select
from motif_scan.patterns import (
    Compiled,
    Motif,
    PatternError,
    compile_motif,
    parse_iupac,
    parse_prosite,
)
from motif_scan.reader import Record, SequenceError, read_records
from motif_scan.search import Hit, count_uncertain, find
from motif_scan.summary import (
    Distribution,
    MotifCount,
    RecordScan,
    Scan,
    distribute,
    expected_count,
    tally,
)

__version__ = "1.0.0rc3"

__all__ = [
    "Alphabet",
    "Compiled",
    "Distribution",
    "Hit",
    "Motif",
    "MotifCount",
    "NUCLEOTIDE",
    "PROTEIN",
    "PatternError",
    "Record",
    "RecordScan",
    "Scan",
    "SequenceError",
    "__version__",
    "all_motifs",
    "categories",
    "compile_motif",
    "count_uncertain",
    "detect",
    "distribute",
    "expected_count",
    "find",
    "load",
    "parse_iupac",
    "parse_prosite",
    "read_records",
    "reverse_complement",
    "select",
    "tally",
]
