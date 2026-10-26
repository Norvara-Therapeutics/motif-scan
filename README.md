# motif-scan

Find sequence motifs from the command line: Kozak contexts and TATA boxes,
cryptic splice sites and poly(A) signals in a payload, Rep binding elements in
an ITR, glycosylation sequons and protease sites in a capsid. Patterns are
written as IUPAC degenerate codes, as PROSITE patterns, or as regular
expressions; mismatches are allowed where that means something; both strands
are searched; and the report is a single HTML file that opens from a local disk
with no network connection.

Pure Python, no third-party dependencies, 40 regulatory and protein elements
built in - 26 nucleotide, 14 protein - each carrying the paper its consensus
comes from.

## Install

Install the prerelease from PyPI:

```bash
python -m pip install --pre norvara-motif-scan
```

The Python package and command keep the shorter project name:

```python
import motif_scan
```

```bash
motif-scan --version
```

## Quickstart

```bash
motif-scan scan -i examples/aav2.fasta -m AAV-RBE,AAV-trs,TATA,polyA --limit 6
```

```
> NC_001401.2  4,679 bp  nucleotide  linear
  Adeno-associated virus 2, complete genome

  MOTIF            PATTERN              HITS     +     -  EXACT   PER KB  EXPECTED  OBS/EXP
  AAV-RBE          GCTCGCTC                8     4     4      8     1.71       0.2    41.93
  AAV-trs          GTTGGCC                 2     1     1      2     0.43       0.7     2.81
  TATA             TATAWAWR                0     0     0      0     0.00       0.6     0.00
  polyA            AATAAA                  2     2     0      2     0.43       1.5     1.33

  AAV-RBE          *                                             **  8
  AAV-trs           |                                            |   2
  TATA                                                               0
  polyA                                  |                      |    2
                   [----------------------------------------------[
                   1                                          4,679

      START       END STRAND  MOTIF            MATCH                    MM
         22        29 +       AAV-RBE          GCTCGCTC                 0
         26        33 +       AAV-RBE          GCTCGCTC                 0
         93       100 -       AAV-RBE          GCTCGCTC                 0
         97       104 -       AAV-RBE          GCTCGCTC                 0
        120       126 -       AAV-trs          GTTGGCC                  0
      2,182     2,187 +       polyA            AATAAA                   0
  ... 6 further hits not shown; raise --limit or write them with --output
```

The Rep binding element appears four times at each end of the genome, twice on
each strand, because the inverted terminal repeat carries two GCTC repeats and
is itself a palindrome. The terminal resolution site appears once per end,
because Rep needs it on the nicked strand only. Both are where the literature
puts them, and the bundled tests assert exactly that.

The TATA box count is zero, and that is the interesting part.

## Mismatches, and why the expected column is there

AAV2 has three promoters, named p5, p19 and p40 for their positions in map
units. Each has a TATA box. None of them matches the TATAWAWR consensus
exactly:

```bash
motif-scan scan -i examples/aav2.fasta -m TATA --mismatches 1 --no-map -f tsv
```

```
record        motif  start  end    strand  length  match     mismatches
NC_001401.2   TATA   255    262    +       8       TATTTAAG  1
NC_001401.2   TATA   843    850    +       8       TATTTAAG  1
NC_001401.2   TATA   1823   1830   +       8       TATAAGTG  1
...
```

255, 843 and 1,823 are 5%, 18% and 39% of the way along a 4,679 base genome.
That is p5, p19 and p40.

Allowing a mismatch is also how a scan starts finding things that are not
there. A seven-base motif with one mismatch allowed matches about one window in
750, so on a few kilobases of either strand it will match a dozen times
whatever the sequence is. That is what the `EXPECTED` and `OBS/EXP` columns are
for:

```
  MOTIF            PATTERN              HITS  ...  EXPECTED  OBS/EXP
  AAV-RBE          GCTCGCTC               16            4.6     3.51
  AAV-trs          GTTGGCC                13           15.0     0.87
```

With one mismatch allowed the terminal resolution site is found thirteen times
and expected fifteen, so twelve of those thirteen are noise and the count
carries no information. The expectation is computed from the record's own
residue composition, position by position over the pattern, with the mismatch
budget taken into account exactly. It assumes residues are independent of one
another, which real sequence is not, so it is a yardstick for reading a count
and not a test. No p-value is reported, because none would be honest.

Both figures are over the same windows. A window carrying an ambiguity code
cannot yield a hit, so it is left out of the expectation as well as out of the
count - otherwise one adenine followed by nine Ns would report one hit against
ten expected and read as a tenfold depletion of a motif that matched every
position it could be assessed at.

## Protein

The same scanner, the same report, PROSITE syntax:

```bash
motif-scan scan -i examples/aav2_vp1.fasta -m PLA2,N-glycosylation,furin-site \
    --sequence --line-width 70 --no-map
```

```
> P03135  735 aa  protein  linear
  Capsid protein VP1, adeno-associated virus 2 (isolate Srivastava/1982)

  MOTIF            PATTERN              HITS  EXACT   PER KB  EXPECTED  OBS/EXP
  PLA2             H-D-x(2)-Y              1      1     1.36       0.0    21.43
  N-glycosylation  N-{P}-[ST]-{P}          6      6     8.16       6.4     0.93
  furin-site       R-x-[KR]-R              1      1     1.36       0.2     6.42

        1 MAADGYLPDWLEDTLSEGIRQWWKLKPGPPPPKPAERHKDDSRGLVLPGYKYLGPFNGLDKGEPVNEADA
       71 AALEHDKAYDRQLDSGDNPYLKYNHADAEFQERLKEDTSFGGNLGRAVFQAKKRVLEPLGLVEEPVKTAP
              aaaa>
```

The phospholipase A2 catalytic motif sits at residue 75, inside the VP1 unique
region that VP2 and VP3 do not carry - which is the whole reason VP1 is needed
for endosomal escape. The six glycosylation sequons sit at an entirely
unremarkable six for a protein of this length and composition, which the
`OBS/EXP` column of 0.93 says plainly: a sequon is a three-residue pattern, and
finding six in 735 residues is what chance does.

A protein has one strand, so `--strand reverse` and `--circular` are refused
rather than ignored, and the strand columns are left out of the table.

## Patterns

Three notations, and `motif-scan codes` prints this summary along with every
residue code:

| Syntax | Example | Fixed length | Mismatches | Strands |
| --- | --- | --- | --- | --- |
| `iupac` | `TATAWAWR`, `CANNTG` | yes | yes | both |
| `prosite` | `N-{P}-[ST]-{P}`, `[ST]-x(2)-[DE]` | when repeats are exact | when fixed | both |
| `regex` | `A{8,}`, `(GGC){3,}` | no | no | both |

PROSITE syntax is the full notation: `A` for a residue, `[ST]` for any of a
set, `{P}` for any except a set, `x` for any residue, `x(2)` and `x(2,4)` for
repeats, `<` and `>` for the termini. A pattern whose repeats are all exact
counts has a fixed length, and is searched like an IUPAC one.

A regular expression is compiled exactly as written, with nothing wrapped
around it, so its own groups keep their numbers and its own inline flags stay
where you put them: `([AC])\1`, `(?P<x>A)(?P=x)`, `(?i)acgt` and `AC(?=G)` all
work as they do in Python. Overlapping occurrences are found by testing each
start position rather than by rewriting your pattern. A pattern containing `=`
is not mistaken for `NAME=PATTERN`, because only text that could be a name is
read as one.

Matching is deliberately asymmetric. An ambiguity code in the **pattern**
widens what the pattern accepts, so `N` matches any base. An ambiguity code in
the **sequence** means the residue is not known, so it matches only a position
that covers everything it might be. Those windows are not hits; they are
counted separately as uncertain, and `--include-uncertain` lists them. A run of
Ns is compatible with any short motif, and reporting it as a finding would send
somebody to redesign a cassette around an element that was never found.

Mismatches never come out of that budget either. A position where the sequence
is ambiguous is uncertain, not mismatched - otherwise every run of Ns would be
a hit for any motif short enough.

## Circular sequences

`--circular` searches across the join, as a plasmid requires:

```bash
motif-scan scan -i plasmid.fasta --all --circular
```

A hit that spans the join is reported with its end smaller than its start, and
never covers more residues than the molecule holds: a pattern of variable
length such as `A(2,8)` is capped so that it cannot match its way around the
circle and past its own start, which would leave the coordinates and the
matched text describing different things.

A pattern anchored to a terminus, or one of unbounded length such as `A+`,
cannot be searched on a circle and is refused with a reason rather than quietly
skipped.

## The report

```bash
motif-scan scan -i examples/aav2.fasta --all --mismatches 1 --circular \
    --report aav2.html
```

One file, no JavaScript, nothing loaded from the network, byte-deterministic
for a given scan. It carries the parameters, a per-motif table with the counts
and their expectations, a binned distribution chart per motif, a scaled map, the
residues themselves with every hit highlighted in place, and the hits table.

The hits table filters without scripting: the controls are checkboxes and the
rules are CSS, so filtering by motif, by strand, or by whether every position
matched works offline, prints, and survives a strict content-security policy.
The one thing that costs is sorting, which belongs to `--sort` on the command
line. The counts above the table are the counts of the whole scan and do not
follow the filter, which the report says in as many words.

## Overlapping hits

Every occurrence is reported by default, including occurrences that overlap one
another, because a motif such as `CG` or `CANNTG` tiles along real sequence and
a scanner that resumed after each match would report every other one.

`--no-overlapping` keeps only the leftmost of each run instead, strand by
strand, so an inverted repeat still reports both of its arms. Overlap is
decided over the residues a hit actually covers, which on a circular sequence
is not the same as comparing positions: a hit spanning the join starts after
everything kept so far and still shares residues with the first hit of the
molecule.

Asking for it also withholds the expected counts, and the output says so. The
expectation counts every window a motif matches, while this reports the
leftmost of each run, so the ratio of the two would move with the reporting
choice rather than with the sequence - five alanines hold four `AA` windows and
would be reported as two hits, reading as a twofold depletion of something
exactly as common as chance predicts.

## Your own motifs

```bash
motif-scan scan -i examples/aav2.fasta --motifs-file examples/payload_motifs.tsv
```

A motif file is a TSV needing only `name` and `pattern`; `syntax`, `alphabet`,
`category`, `description` and `source` are optional. `motif-scan motifs -f tsv`
writes the built-in catalogue in that format, ready to edit and read back.

The file has to be rectangular and its column names distinct. Neither
restriction is pedantry: a repeated column name silently keeps the last of the
two, so a header of `name pattern name` renames every motif in the file, and a
row with one field too many has the extra one filed under no column and
dropped. Both are refused by line number.

## Output

| Flag | What it writes |
| --- | --- |
| `-f table` | the aligned view above (the default) |
| `-f tsv` | every hit, tab-separated, on standard output |
| `-f json` | motifs, counts, hits and the distribution, in full |
| `-o PATH` | every hit as TSV, never truncated |
| `--counts PATH` | the per-motif counts, with expectations |
| `--report PATH` | the self-contained HTML report |

`--limit` truncates the table and the report, never the TSV or the JSON. Any
option that cannot affect what a run writes is reported on standard error
rather than ignored.

An output that would overwrite an input is refused before anything is read -
and the motif file counts as an input, because it is a file you wrote by hand.
So are two outputs aimed at one file. Paths are compared as files rather than
as strings: `out.tsv` and `./out.tsv` are the same file, and so is a hard link
under another name.

## What this is not

It is a consensus matcher, and a consensus is a blunt instrument. It
over-reports, because `CANNTG` occurs every few hundred bases by chance, and it
under-reports, because a real binding site that differs from the consensus at
one position is still a real binding site. Read every count against its
expectation.

For transcription factor binding specifically, a position weight matrix from
[JASPAR](https://jaspar.elixir.no) or HOCOMOCO is a better instrument than any
string in this catalogue, and nothing here is a substitute for one: there is no
matrix model, no score, no threshold, and no background model beyond the
composition of the record in front of it. For protein domains,
[InterProScan](https://www.ebi.ac.uk/interpro/) and the full PROSITE release
search thousands of patterns and profiles rather than the handful bundled here.

What this does instead is run locally over unpublished or patient-derived
sequence with no network and no dependencies, in both alphabets, with
mismatches, across the join of a circle, and emit a report that will still open
in ten years.

## Examples

The bundled sequences are fetched from public databases by
`examples/fetch_examples.py`, which is executable provenance rather than an
assertion:

| File | Source | Why |
| --- | --- | --- |
| `aav2.fasta` | NCBI `NC_001401.2` | ITRs carry the RBE and the *trs*; p5, p19 and p40 carry TATA boxes |
| `aav2_vp1.fasta` | UniProt `P03135` | VP1u carries the phospholipase A2 HDxxY motif |
| `hbb.fasta` | NCBI `NG_059281.1` | CCAAT at -76, CACCC at -90 and -105, and a variant TATA box at -32 |

The documented positions in that last column are assertions in
`tests/test_reference_examples.py`, so a change that breaks the coordinate or
strand convention breaks the suite.

## Development

```bash
python -m pip install -e ".[dev]"
python -m pytest
python -m pyflakes src/motif_scan tests scripts examples
python scripts/check_wheel.py
```

`check_wheel.py` builds the distributions, installs the wheel into a throwaway
virtual environment with `PYTHONPATH` stripped, and runs the whole suite plus a
command-line smoke test from outside the repository, which is the only way to
catch a module that never made it into the package.

## Licence

PolyForm Noncommercial 1.0.0. See [LICENSE](LICENSE).
