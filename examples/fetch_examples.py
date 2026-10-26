"""Fetch the bundled example sequences from their public sources.

The files in this directory are committed, so nothing here needs to run for the
examples in the README to work. It exists so that the provenance of each file
is executable rather than asserted: run it and the committed files are
reproduced byte for byte, or the difference is a change in the source database.

    python examples/fetch_examples.py

Three sequences are fetched, each chosen because what it should contain is
documented elsewhere and can be checked:

* ``aav2.fasta`` - the AAV2 genome. Its inverted terminal repeats carry the Rep
  binding element and the terminal resolution site, and the p5, p19 and p40
  promoters carry TATA boxes.
* ``aav2_vp1.fasta`` - the AAV2 VP1 capsid protein. Its unique N-terminal region
  carries the HDxxY catalytic motif of a phospholipase A2 domain.
* ``hbb.fasta`` - the human beta-globin gene with 5 kb of flank on each side.
  Transcription starts at position 5001, which puts the CCAAT box at -76 and
  the duplicated CACCC box at -90 and -105. Its TATA box is the variant
  CATAAAA, which an exact search for TATAWAWR does not find and a search
  allowing one mismatch does.
"""

from __future__ import annotations

import sys
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent

#: A courtesy identifier, as the NCBI guidelines ask for.
TOOL = "motif-scan"

#: Line width of the written files, which is what makes the output comparable
#: with the committed one.
WIDTH = 60

NUCLEOTIDE = (
    (
        "aav2.fasta",
        "NC_001401.2",
        "NC_001401.2 Adeno-associated virus 2, complete genome",
    ),
    (
        "hbb.fasta",
        "NG_059281.1",
        "NG_059281.1 Homo sapiens hemoglobin subunit beta (HBB), RefSeqGene; "
        "the gene occupies 5001..6608",
    ),
)

PROTEIN = (
    (
        "aav2_vp1.fasta",
        "P03135",
        "P03135 Capsid protein VP1, adeno-associated virus 2 "
        "(isolate Srivastava/1982)",
    ),
)


def fetch(url: str) -> str:
    """Read a URL as text."""
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read().decode()


def residues(fasta: str) -> str:
    """The sequence of a single-record FASTA document, uppercased."""
    lines = [line.strip() for line in fasta.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise SystemExit("the download is not FASTA; the source may have moved")
    return "".join(lines[1:]).upper()


def write(name: str, header: str, sequence: str) -> None:
    """Write one record, wrapped, with Unix line endings."""
    body = "\n".join(
        sequence[index: index + WIDTH] for index in range(0, len(sequence), WIDTH)
    )
    (HERE / name).write_text(f">{header}\n{body}\n", encoding="utf-8", newline="\n")
    print(f"{name}: {len(sequence)} residues")


def main() -> int:
    for name, accession, header in NUCLEOTIDE:
        query = urllib.parse.urlencode(
            {
                "db": "nuccore",
                "id": accession,
                "rettype": "fasta",
                "retmode": "text",
                "tool": TOOL,
            }
        )
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}"
        write(name, header, residues(fetch(url)))

    for name, accession, header in PROTEIN:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
        write(name, header, residues(fetch(url)))

    print("\nSequences are from NCBI Nucleotide and UniProtKB.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
