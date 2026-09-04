# Stage 9 — the genomic and construct-safety arm

Written before `stages/construct_safety.py` exists. The criteria in §6 are fixed
here and committed before any result.

This closes four of the six named gaps in deliverable 6: recombination-prone
regions, cryptic splice sites, unwanted open reading frames, and sequence
repeats. It needs no new data source, no model and no external call — every
input is already emitted by Stage 6.

---

## 1. The constraint that decides the whole design

Stage 6 emits a DNA map, and §2 of its specification is explicit about what that
map is:

> The DNA is a **map, not an optimised sequence**. It is produced by reverse
> translation under one fixed codon per amino acid … It is not a codon-optimised
> ordering sequence and must not be read as one.

Every residue of a given kind therefore receives the same codon: every leucine is
`CTG`, every serine `AGC`, every alanine `GCC`. **This is not the sequence anyone
would order, and a nucleotide-level finding on it is not automatically a finding
about a manufacturable construct.**

The consequences are specific and they cut both ways:

- **Repeats are inflated.** Two occurrences of the same short peptide produce
  *identical* nucleotides under a fixed codon table, where a real ordering
  sequence would deliberately diversify them. A nucleotide repeat count on this
  map is an upper bound, not a measurement.
- **Splice motifs are arbitrary.** Whether `GT` falls at a given position is a
  property of the codon assignment. A different table moves every motif. A site
  found here may not exist in the ordered construct, and a site in the ordered
  construct may not appear here.
- **Alternate-frame and reverse-strand reading is arbitrary** for the same
  reason.

An arm that reported these four as flat findings would be **stating properties of
an arbitrary encoding as properties of a therapeutic**. That is the failure this
repository keeps finding, and it would be a particularly bad instance because
every number would look like a real sequence analysis.

## 2. What is actually invariant, and what is not

The design turns on separating the two, and every finding carries its class.

**`CODON_INVARIANT`** — true under any codon assignment, because it is a property
of the protein or of the domain layout:

- **Repeated parts in one construct.** A dual design contains the CD8A leader,
  hinge and transmembrane segment *twice*, by construction. That is read from the
  domain map, exactly, and no codon choice changes it. It is also the dominant
  real recombination hazard in these constructs: long, perfectly identical
  stretches that a vector can loop out.
- **Repeated peptides.** A peptide occurring twice is repeated under every
  encoding. The nucleotides may differ in an ordered sequence; the homology
  pressure does not vanish.
- **Internal methionines**, which are candidate internal initiation sites in the
  coding frame regardless of encoding.

**`MAP_SPECIFIC`** — a property of this particular reverse translation:

- Nucleotide-level direct and inverted repeats.
- Splice donor and acceptor motifs.
- Open reading frames in frames 2 and 3 and on the reverse strand.

**Both are reported. Neither is dropped, and neither is presented as the other.**
A `MAP_SPECIFIC` finding is a bound on what an ordering sequence might contain
and is labelled as one, which is the same rule the platform already applies to
risk figures that rest on an unmeasured organ.

## 3. What the arm computes

Per construct, from `construct.dna`, `construct.amino_acid_sequence` and
`construct.segments` — nothing else.

| finding | class | rule |
| --- | --- | --- |
| `repeated_part` | invariant | two segments sharing accession and residue range |
| `repeated_peptide` | invariant | a peptide of ≥ `PEPTIDE_REPEAT_MIN` residues occurring more than once |
| `internal_methionine` | invariant | an `M` after position 1, with its residue index |
| `direct_repeat` | map-specific | a nucleotide substring of ≥ `NT_REPEAT_MIN` occurring more than once |
| `inverted_repeat` | map-specific | a substring of ≥ `NT_REPEAT_MIN` whose reverse complement also occurs |
| `splice_donor` | map-specific | `GT` at an exon-boundary-like context, scored by the canonical donor consensus |
| `splice_acceptor` | map-specific | `AG` preceded by a pyrimidine tract of ≥ `TRACT_MIN` |
| `alternate_orf` | map-specific | an ATG-to-stop run of ≥ `ORF_MIN` codons in frame 2 or 3 |
| `reverse_orf` | map-specific | the same on the reverse complement |
| `homopolymer` | map-specific | a run of ≥ `HOMOPOLYMER_MIN` identical bases |

Every threshold is a named constant, fixed in this document before any run, and
none is tuned afterwards:

| constant | value | why this value |
| --- | --- | --- |
| `PEPTIDE_REPEAT_MIN` | 8 residues | below this, repeats are ubiquitous in any protein and carry no signal |
| `NT_REPEAT_MIN` | 24 bp | the length below which homologous recombination is not efficiently promoted in the literature this platform can cite; also 8 codons, so it is the nucleotide image of the peptide floor |
| `TRACT_MIN` | 8 pyrimidines | the short end of a functional polypyrimidine tract |
| `ORF_MIN` | 30 codons | shorter reading frames occur by chance in any sequence of this length |
| `HOMOPOLYMER_MIN` | 8 bases | the run length at which polymerase slippage becomes a documented synthesis problem |

**These are stated as conventional working values and not as measurements this
platform made.** Nothing here calibrates them against an outcome, because the
platform has no outcome data to calibrate against. Where that matters, §6's
criteria test the *detector*, never the threshold.

## 4. What it does not do

- **It does not gate.** No construct is blocked, no verdict changes, no ranking
  moves. The arm reports; Stage 9's existing risk gate is untouched. A
  construct-safety finding is information for a reader, and turning it into a
  block would be setting a tolerance this specification has no basis to set.
- **It does not modify the construct.** No sequence is rewritten, no site
  removed, no codon changed. Fixing a cryptic splice site means changing a
  codon, which would make Stage 6 emit a sequence that is no longer the map it
  says it is.
- **It does not claim to have analysed a manufacturable sequence.** §1.
- **It does not add a data source.** Every input is already in memory.

## 5. Where it lives

`car_pipeline/stages/construct_safety.py`, pure functions over a construct.
Wired into `stage9.gate`, which gains `construct_safety` on each `SafetyRecord`
for constructs that assembled. The candidate package's `safety` section carries
it, which closes the gap probe `(key, "safety", "construct_safety")`.

**That probe will trip Q6 on the first run after this lands, and that is the
mechanism working.** The gap table declares the element missing; once it is
present the declaration is false and the criterion says so. The table entry is
removed in the same change, and the run that follows is the evidence.

## 6. Rejection criteria — fixed before the run

Each detector is tested against a **known answer** — a sequence built to contain
exactly the thing being looked for, and one built to contain none of it. That is
the pattern the developability stage already uses for its charge model, and it
is the only way to test a detector on a platform with no labelled data.

| id | trips when |
| --- | --- |
| **S8** | a 60 bp direct repeat planted at known positions is not reported at those positions, or a control sequence with no repeat above the floor reports one. Both directions, because a detector that finds everything and one that finds nothing both pass a single-sided check |
| **S9** | a canonical donor `GTAAGT` and an acceptor with a 10-pyrimidine tract, both planted, are not reported; or a control containing neither reports one |
| **S10** | a 40-codon ORF planted in frame 2 is not reported, or the construct's own coding frame is reported as an alternate ORF — the annotated frame is never an unwanted one |
| **S11** | any finding lacks a class, or re-running the whole arm under a **different codon assignment** changes any `CODON_INVARIANT` finding. This is the criterion that enforces §2: the labels are not commentary, they are a claim that can be falsified by re-encoding the same protein |
| **S12** | a construct whose domain map repeats a part does not report `repeated_part`, or the adaptor designs — which repeat none — report one |
| **S13** | any construct's sequence, DNA, segments or verdict differs before and after the arm runs |
| **S14** | the arm reports a finding for a construct that assembled nothing, or omits the section for one that did |

**S11 is the one that matters.** S8 through S10 test that the detectors work.
S11 tests that the honesty label is true, by re-encoding the same protein under a
shuffled codon table and requiring the invariant findings to be identical and the
map-specific ones to be free to move. A label nobody can falsify is decoration.

### Explicitly not grounds for rejection

- A construct carrying many map-specific findings. That is expected: a fixed
  codon table maximises nucleotide self-similarity, and §1 says so in advance.
- Zero `repeated_part` findings across the shipping designs. All five are
  adaptors, which repeat no part. S12's positive half is a known-answer control
  for exactly that reason.
- The thresholds in §3 being conventional rather than derived. §3 says so.
