# docking-ai4 — Protein-ligand binding interaction descriptions

A pipeline that takes protein-ligand crystal structures, extracts geometric
contacts (hydrogen bonds, salt bridges, hydrophobic shells, pi-stacking,
covalent linkages, oxyanion holes, water bridges), and generates
crystallography-paper-style narrative descriptions via an LLM grounded in
per-target biological context.

The project started as a single example — three TEM-1 β-lactamase / boronic
acid inhibitor complexes from [Ness et al., *Biochemistry* 2000] — and now
covers **255 CASP protein-ligand complexes** from CASP15 and CASP16 with
generated narratives for each.

## What's here

```
docking-ai4/
├── describe_binding.py          # the whole pipeline, ~900 lines
├── CLAUDE.md                    # architecture notes for future Claude sessions
├── README.md                    # this file
│
├── 1erm.pdb / 1ero.pdb / 1erq.pdb   # TEM-1 beta-lactamase co-crystals (Ness 2000)
├── 1erm.md  / 1ero.md  / 1erq.md    # generated narratives for those
├── 1erm.json / 1ero.json / 1erq.json
├── comparison.md                     # hand-written cross-structure comparison
├── paper_vs_tool.md                  # validation vs. the published paper
│
├── casp15/
│   ├── targets.tsv              # 27-row target -> ligand mapping
│   ├── run_casp15.py            # driver
│   └── out/                     # 22 × {json,txt,md}
│
└── casp16/
    ├── targets.tsv              # 233-row target -> (RCSB, series, affinity, covalent) map
    ├── build_targets.py         # joiner: CASP archive + affinity CSV + Zenodo SI
    ├── run_casp16.py            # driver (split protein+ligand input mode)
    ├── context/
    │   ├── chymase.txt          # per-family biology context
    │   ├── cathepsin_g.txt      # injected into LLM prompts
    │   ├── autotaxin.txt
    │   └── mpro.txt
    └── out/                     # 233 × {json,txt,md}
```

Data too large for the repo lives in `~/data/vaults/casp/`:

- `casp15_ligands/targets_ligand/Targets_ligand/` — 27 CASP15 `*_lig.pdb`
  files (5.3 MB, the full crystal complexes)
- `casp16/pharma_ligands/exper_struct/L{1,2,3,4}000/` — 233 per-target
  directories, each containing `protein_aligned.pdb` + one or more
  `ligand_*.pdb` files (57 MB total)
- `casp16/zenodo/` — 5 Zenodo SI files from Skrinjar et al. 2025: the
  authoritative compound ID → RCSB PDB → chemical series → IC50 tables
- `casp16/papers/` — the Gilson assessment paper PDF

## Quick start

One-shot on the TEM-1 example:

```bash
python3 describe_binding.py 1erm.pdb --ligand BJI           # full LLM narrative
python3 describe_binding.py 1erm.pdb --ligand BJI --no-llm  # plain text
python3 describe_binding.py 1erm.pdb --ligand BJI --json    # raw interaction data
```

Regenerate the CASP15 or CASP16 outputs:

```bash
# CASP15: geometry only (no API spend)
python3 casp15/run_casp15.py

# CASP15: geometry + LLM narratives via o4-mini
python3 casp15/run_casp15.py --with-llm -m o4-mini

# CASP16 per-family
python3 casp16/run_casp16.py --only-family chymase --with-llm
python3 casp16/run_casp16.py --skip-family autotaxin --with-llm  # 44 targets
python3 casp16/run_casp16.py --only-family autotaxin --with-llm  # 189 targets
python3 casp16/run_casp16.py --only L1003,L4004 --with-llm       # targeted

# Rebuild the CASP16 metadata join from the Zenodo SI
python3 casp16/build_targets.py
```

Total API cost for a full pass is around $12 at o4-mini (255 narratives ×
~$0.05 each). Dropping to `-m gpt-4.1-mini` would cut that by a further ~3×
at some cost in accuracy.

## Pipeline architecture

Four stages, all inside `describe_binding.py`:

1. **Parse & classify atoms.** Walk the Bio.PDB structure, bucket atoms into
   ligand / protein / water. Non-ligand HETATMs that aren't ions or
   crystallization artifacts are treated as part of the protein (handles
   modified residues like covalently-modified catalytic serines).

2. **Geometric interaction detection.** A single `NeighborSearch` over
   protein+water atoms powers seven detectors with hardcoded cutoffs at the
   top of the module:

   | Detector | Cutoff | What it finds |
   |---|---|---|
   | covalent | 2.0 Å | heavy-atom close contacts, skipping hydrogens |
   | hydrogen bond | 3.5 Å | donor/acceptor element check on N/O/S |
   | salt bridge | 4.0 Å | charged ligand N/O ↔ Arg/Lys/His or Asp/Glu |
   | hydrophobic | 4.5 Å | carbon–carbon contacts |
   | pi-stacking | 5.5 Å | centroid–centroid between ligand rings and protein aromatic side chains, with parallel / T-shaped / offset classification from plane angle |
   | water bridge | 3.5 Å per leg | waters within H-bond range of both ligand and protein |
   | oxyanion hole | 3.5 Å | ≥2 backbone NH within H-bond range of a ligand oxygen |

   Ligand-side aromatic rings are detected by building a 1.85-Å connectivity
   graph over ligand heavy atoms, enumerating 5- and 6-cycles, and keeping
   only those whose atoms lie in a plane (smallest SVD singular value
   < 0.25 Å).

3. **Summarization.** Collapses the raw per-atom contact lists into
   per-residue groupings suitable for prose generation — H-bonds bucketed by
   residue, hydrophobic contacts reduced to (residue, closest distance,
   count), pi interactions deduplicated.

4. **Narrative generation.** A long crystallography-paper-style prompt
   instructs the LLM to write 2–4 paragraphs in past tense, passive voice,
   with distances in parenthesized Ångstroms and banned editorializing
   vocabulary. Optional `--context` flag injects authoritative biology
   context (catalytic residues, ligand chemical class, covalent/noncovalent
   flag) that gets marked as "ground truth, do not contradict" in the
   prompt. OpenAI is tried first with Anthropic as a fallback.

## The original example: TEM-1 β-lactamase + boronic acid inhibitors

The three PDB structures `1erm` / `1ero` / `1erq` are the transition-state
analogue inhibitors from [Ness et al., *Biochemistry* 39, 5312–5321 (2000)].
All three are covalent tetrahedral boronate adducts to Ser70.

- `1erm.md` / `1ero.md` / `1erq.md` — individual generated narratives
- `comparison.md` — hand-written side-by-side comparing the three compounds
  across covalent bond length, oxyanion hole fit, Glu166 / Arg243 /
  Tyr105 / Ser130 contacts, and water-mediated bridges
- `paper_vs_tool.md` — validation pass checking our computational output
  against the published Table 2 and prose description

The `paper_vs_tool.md` document flagged a labeling inconsistency between our
`comparison.md` and the Ness paper's compound numbering (we called 1ERM
"compound 2" when the paper calls it "compound 3"); worth knowing if reading
both. That file also catalogues which interactions our tool recovered and
which it missed — the main miss is that our pi-stacking cutoff is loose and
our geometric analysis can't see the oxyanion-hole resonance behavior.

## CASP15 extension

[CASP15] included 20 ligand-binding prediction targets. Ground-truth
structures are on the CASP prediction center (not RCSB):

```
predictioncenter.org/download_area/CASP15/targets/casp15.targets.ligands.ALL_09.18.2025.tar.gz
```

27 files in that archive; we keep 22 as runnable protein-ligand complexes.
Skipped:

- **4 RNA targets** (R1117, R1117v2, R1126, R1136v1) — out of scope for a
  protein-focused interaction analyzer
- **H1135 (SUN1-KASH6)** — only K⁺/Cl⁻ ions, no small-molecule ligand

### Targets and families

| Family | Targets | Ligand | Notes |
|---|---|---|---|
| Huc cofactor | H1114 | FCO (+3NI, F3S, MQ7) | Multi-metal cluster, treat FCO as representative |
| RuvB nucleotide | H1171v1/v2, H1172v1–v4, T1170 | AGS / ADP | 6 variants across ATP and ADP states |
| FoxA ferrioxamine | T1118v1 | LIG | Custom resname |
| MfnG methyltransferase | T1124 | SAH | Plus free Tyr substrate |
| NATA1 acyltransferase | T1127, T1127v2 | COA | |
| Peptidoglycan hydrolase | T1146 | NAG | No RCSB, prediction-center only |
| LysM domain | T1152 | NAG | |
| MRP4 transporter | T1158v1/v2/v3/v4 | XPG / P2E / XH0 / ATP | 4 different substrates |
| Phage tail protein | T1181 | OAA,OAb | Two OAA conformers merged |
| CTX-M-14 β-lactamase | T1186 | LIG (dicloxacillin) | Covalent to Ser43 |
| Nictaba lectin | T1187 | NAG | Tobacco lectin carbohydrate complex |
| *C. perfringens* chitinase | T1188 | DW0 | Four-Trp box binder |

### CASP15 highlights from the run

- **T1186 dicloxacillin (β-lactamase adduct)** — pipeline correctly found a
  1.40 Å covalent bond from the ligand's β-lactam C to Ser43 OG, exactly as
  expected for a serine β-lactamase acylation.
- **T1170 / H1172 RuvB hexamers** — chain-aware residue labels
  (`Thr52/A`, `Arg203/D`) disambiguate the six nucleotide binding sites.
  Narratives correctly identify the P-loop Lys50/Arg156/Arg203 guanidinium
  finger around γ-phosphate.
- **T1124 MfnG + SAH** — proper ligand-side aromatic ring detection picks up
  the SAH adenine base and finds its parallel π-stacking with Phe253 in
  both protomers of the dimer (4.03 Å, 4.32 Å).
- **T1188 chitinase + DW0** — four-tryptophan aromatic box (Trp35, Trp120,
  Trp252, Trp420) correctly recovered; this is the classical
  carbohydrate-binder signature.
- **T1187 Nictaba lectin + NAG trisaccharide** — all three sugars in the
  trimer get individual contact descriptions, Trp16/Trp42/Trp152/Tyr22
  stacking is captured.

## CASP16 extension

CASP16 introduced the first pharmaceutical ligand-binding challenge with
**263 protein-ligand complexes** across four drug targets from Roche and
Idorsia. We cover **233 of them** — the pose-prediction structures with
public ground-truth PDBs. (The remaining 30 are affinity-only with no
structure.)

Ground-truth location (the local vault had incorrectly concluded these were
proprietary; they are public but not on RCSB):

```
predictioncenter.org/download_area/CASP16/targets/pharma_ligands/
├── L1000_exper_struct.tar.gz (chymase, 17 targets)
├── L2000_exper_struct.tar.gz (cathepsin G, 2 targets)
├── L3000_exper_struct.tar.gz (autotaxin, 189 targets)
└── L4000_exper_struct.tar.gz (Mpro, 25 targets)
```

Each per-target archive expands to `L{1,2,3,4}000_prepared/LNNNN/` with
`protein_aligned.pdb` + one or more `ligand_<compound>_<chain>_<copy>.pdb`
files. Multi-ligand entries contain the drug-like compound plus **incidental
ligands** (sulfate, DMSO, PEG, ethylene glycol, bromide, acetate) — small
cosolutes that happen to sit within 4.5 Å of the drug-like ligand in the
crystal. Our driver filters incidentals via `SKIP_RESNAMES` and picks the
alphabetically first non-incidental file per directory.

### Target metadata — `casp16/targets.tsv`

`build_targets.py` joins three data sources:

1. The CASP archive directory tree (protein + ligand file paths)
2. The local `L{1,3}000_exper_affinity.csv` (CASP target → compound ID)
3. The Zenodo SI files (compound ID → RCSB PDB → series → IC50 → pKa)

Coverage after the join:

| Family | Count | RCSB mapped | Series labeled |
|---|---|---|---|
| chymase (L1000) | 17 | 17 | 17 |
| cathepsin_g (L2000) | 2 | 2 | 2 |
| autotaxin (L3000) | 189 | 178 | 10 |
| mpro (L4000) | 25 | 25 | 0 |

ATX series labels are sparse because the Zenodo SI only annotates 40 ATX
compounds (20 pyrrolo-pyrrole + 20 phenoxymethyl-imidazole). The remaining
149 pose-only targets have SMILES but no series label. A SMILES-pattern
classifier could fill this gap.

Flagged targets:

- **Disclosed pre-CASP** (excluded from official affinity scoring but still
  in our dataset): L1003, L1010, L1013 (chymase, Roche patents); 14 ATX
  targets; L4006–L4010 Mpro (inadvertently published as 7grh/7gri/7grr/
  7grw/7gs2 before the prediction deadline).
- **Covalent Mpro** (geometric analysis): L4003, L4004, L4013, L4019, L4023
  — all five show a sub-1.77 Å covalent bond to `CYS145/A SG`. Skrinjar
  2025 only flags L4004 as "covalent" because their terminology is
  *cocrystallized vs. soaked*, not *covalent vs. noncovalent*. The other
  four were soaked but still form the thioether.

### Per-family context blocks

`casp16/context/{chymase,cathepsin_g,autotaxin,mpro}.txt` are 4 per-family
biology context blocks written directly from the Tosstorff/Skrinjar dataset
paper. They're injected into the LLM prompt via `--context` to ground the
narrative — catalytic triads, electrostatic pocket character, the
phenylsulfonamide induced-fit Leu107 note, the pyrrolo-pyrrole vs
phenoxymethyl-imidazole Zn-coordination distinction, the covalent L4004
flag.

### CASP16 highlights from the run

- **L4004 covalent Mpro inhibitor.** The pipeline finds a 1.76 Å covalent
  bond from ligand atom `C02` to `CYS145/A SG` and the narrative correctly
  describes it as "formation of a tetrahedral adduct at the catalytic
  dyad." The LLM was not told in the context block that L4004 was covalent
  — it read the distance from the JSON and drew the right conclusion.
- **L1003 chymase phenylsulfonamide.** This is the diagnostic target.
  Skrinjar describes phenylsulfonamides as binding via induced fit to a
  cryptic sub-pocket that opens when Leu107's side chain moves, and
  — unlike the other three chymase series — they do *not* contact Ser203.
  The geometric pipeline found **zero Ser residues in the hydrogen bond
  list** and **five hydrophobic contacts to Leu107**. The narrative
  correctly states "no direct interaction with Ser203 ... in agreement
  with the mechanistic distinction of this inhibitor series" and lists
  Leu107 in the hydrophobic shell. This validates the context-injection
  strategy end-to-end.
- **L1001 chymase indole-carboxylic acid.** Engages the electropositive S1
  pocket as expected: Gly201 backbone, Ser203 OG side chain, Lys49/His66
  salt-bridge character, His66 imidazole π-stacking at 3.65 Å. The
  acidic-warhead + electropositive-pocket pattern Skrinjar describes is
  recovered cleanly.
- **L4006 pre-disclosed Mpro fragment (PDB 7gs2).** JSON is nearly empty
  (zero H-bonds, zero salt bridges, zero π-stacking). The narrative
  correctly reports the zeros, describes the 6 hydrophobic contacts to the
  S2–S4 region, and is appropriately short. A good demonstration that the
  pipeline doesn't fabricate content when data is sparse.

## Pipeline fixes made in this iteration

Eleven fixes to `describe_binding.py` landed while extending to the
CASP15/16 dataset. Each came from a concrete failure mode:

| Fix | Problem it solved |
|---|---|
| `SKIP_RESNAMES` = ions ∪ crystallization additives | T1127 had EPE/MPD contaminants; CASP16 has DMSO/PEG/sulfate incidentals flanking every multi-ligand target |
| Skip hydrogens in covalent/salt detection | CASP files include explicit H atoms from refinement; produced false covalent bonds at ~1.96 Å |
| `--ligand OAA,OAb` comma-separated lists | T1181 has two conformers of oxaloacetate; without merging, the two bond to each other and look like 5 false covalent bonds |
| Chain-aware residue labels auto-enabled | RuvB hexamers and MRP4 dimers were ambiguously labeled; `/A` `/B` `/D` suffixes added when >1 protein chain present |
| Ligand-side aromatic ring detection (planarity SVD) | SAH adenine, ATP adenine, purine cofactors, imidazole rings — none were detected by the old "any ligand C near protein aromatic centroid" heuristic |
| `--context` flag for per-target biology | Per-family and per-target biology needs to reach the LLM without bloating `run_*.py` with f-string prompt construction |
| Missing `os` import | Broke `--context` file loading for every CASP16 LLM call; empty outputs across all 44 chymase/catG/mpro targets in the first background run |
| Stronger anti-hallucination prompt rules | Early narratives ended with fabricated protein-protein H-bond paragraphs. Now partially suppressed (see Known issues below) |
| ATX altloc filename regex | ~20 ATX targets use `ligand_L0R_C_1_A.pdb` / `_B.pdb` alternate-conformation split files; old regex rejected the trailing suffix |
| `UNL`/`UNK`/`UNX` removed from `SKIP_RESNAMES` | 6 ATX targets use the generic "unknown" PDB code as the drug-like ligand resname |
| `max_completion_tokens` 4000 → 16000 for reasoning models | At 4000 on o4-mini, ~20% of larger CASP16 targets came back with empty output because the hidden reasoning budget consumed the whole allowance |

## Quality assessment

Spot-checked 12 narratives spanning both CASP rounds and all four CASP16
families. Scoring rubric:

- **A** — accurate, biology-grounded, no hallucination, good prose
- **B** — mostly accurate, minor issues
- **C** — has hallucinations or factual errors

| Target | Family / class | Score | Key finding |
|---|---|---|---|
| T1124 | MfnG + SAH | A− | Correct adenine H-bonds + Asp225 ribose motif + Phe253 π-stacking |
| T1158v1 | MRP4 + PGE1 | A | Arg946/Arg998 carboxylate pincer + Asp842/Gln845 polar head + correctly zero π |
| T1187 | Nictaba lectin | A− | Trp-box recovery, minor "oxyanion hole" mislabel (lectin misnomer) |
| L1001 | chymase indole-CA | B | Correct S1 + Ser203 + π-stack, hallucinated final paragraph |
| L1002 | chymase 5-ring vinyl | C | Correct + partial occupancy note, 5 fabricated protein-protein bonds |
| **L1003** | chymase phenylsulfonamide | B+ | *Headline success* — no Ser contacts, Leu107 pocket recovered |
| L1014 | chymase indole-CA | B | Weak binder correctly described, hallucinated final paragraph |
| L1017 | chymase indole-CA | B | Strong binder correctly described, hallucinated final paragraph |
| **L4004** | Mpro covalent (cocrystal) | B | 1.76 Å Cys145 bond correctly described as tetrahedral adduct |
| L4019 | Mpro covalent (soaked) | B+ | 1.68 Å + 0.80 partial occupancy, shorter hallucination |
| L4001 | Mpro noncovalent | B− | Correct bidentate His163, worst hallucination (14 fake distances) |
| **L4006** | Mpro fragment (7gs2) | A | Sparse JSON → sparse narrative, zero hallucination |

### The consistent failure mode

**7 of 12 narratives end with a paragraph hallucinating protein-protein
hydrogen bonds** — backbone-backbone contacts and catalytic triad geometry
that are *not* in the JSON. The pattern is identical across targets and
always takes the form:

> "A network of protein–protein hydrogen bonds was maintained in the active
> site near the ligand: Ser203 OG to His66 NE2 (2.74 Å), Lys49 NZ to His66
> O (2.95 Å), …"

This is o4-mini's reasoning pulling structural-paper-genre patterns from
its training data when the ligand-interaction section runs short. I added
a prompt rule forbidding protein-protein descriptions, which partially
suppressed the behavior but didn't eliminate it. Two remedies considered:

1. **Post-process filter** — regex-detect and strip the stock hallucinated
   paragraph. Zero API cost, deterministic.
2. **Stronger prompt rewrite + full regeneration** — ~$10 more at o4-mini,
   ~40 minutes wall time.

Not addressed yet — leaving for later.

### The headline successes

- **L1003** end-to-end proves the context-injection strategy works:
  Skrinjar's "phenylsulfonamides bind Leu107 induced-fit pocket, do not
  contact Ser203" narrative fact propagated from the paper → per-family
  context block → LLM prompt → final narrative, and the geometric data
  independently confirms it.
- **L4004** covalent linkage was correctly described as a tetrahedral
  adduct without the context block being told L4004 is covalent — the
  LLM read the 1.76 Å distance from the JSON and inferred the chemistry
  correctly.
- **L4006** demonstrates that sparse JSON yields sparse prose: when there
  is nothing to say, the model says little, rather than filling with
  fabricated content. This is the desired degenerate case.

## Known issues and future work

- **Hallucinated final paragraph** in ~60% of narratives (see Quality
  assessment above). Post-process filter is the planned remedy.
- **"Oxyanion hole" label for lectins.** The detector fires on any ligand
  oxygen with ≥2 backbone NH within 3.5 Å. For serine proteases and
  β-lactamases this is the real oxyanion hole; for a carbohydrate-binding
  lectin like Nictaba it's a sugar-recognition motif using backbone NH
  donors. The label propagates into the narrative and is wrong for
  non-protease ligand classes. Consider renaming the output field or
  conditioning the label on the protein family.
- **ATX series labels (10/189).** Most ATX targets have no chemical-class
  annotation in the Zenodo SI. A simple SMILES-pattern classifier could
  distinguish pyrrolo-pyrrole from phenoxymethyl-imidazole from
  miscellaneous, and would improve per-target narratives.
- **L4004 `COVALENT_MPRO` set in `build_targets.py` is incomplete.** It
  includes only L4004 (following Skrinjar's "cocrystallized" terminology)
  but the geometric analysis shows L4003/L4013/L4019/L4023 are also
  covalent to Cys145. Should extend the set to all five.
- **Chymase UniProt vs PDB numbering shift.** Skrinjar uses UniProt
  numbering (Ser203, His66, Asp113); the actual PDB files may use a
  different offset. Narratives use whatever numbers appear in the
  coordinate file — worth verifying consistency against Skrinjar's text
  before publishing any narrative.
- **Per-target comparisons** like the original `comparison.md` are a
  natural next deliverable for the three big "sibling" groups:
  - **MRP4 family** (T1158v1/v2/v3/v4): same transporter, four different
    substrates (PGE1, PGE2, DHEAS, ATP). Good for studying substrate
    promiscuity.
  - **RuvB hexamer** (H1171v1/v2, H1172v1–v4, T1170): 6 nucleotide
    binding sites across AGS and ADP states.
  - **Mpro covalent series** (L4003/L4004/L4013/L4019/L4023): 5 different
    covalent warheads attaching to the same Cys145.
- **Water-mediated interactions** are zero across all CASP15/16 targets
  because the CASP `*_lig.pdb` files have no crystal waters — the
  organizers stripped them. For datasets that retain waters, the detector
  already works.
- **Pi-stacking classification** is new and useful but not validated
  quantitatively. The plane-angle thresholds (30° / 60° for
  parallel/T-shaped/offset) are defaults, not tuned.

## References

Papers consumed while building this:

- **Ness et al. (2000)** "Structure-based design guides the improved
  efficacy of deacylation transition state analogue inhibitors of TEM-1
  β-lactamase." *Biochemistry* 39, 5312–5321. The source of the 1erm/1ero/
  1erq example structures.
- **Tosstorff, Skrinjar et al. (2026)** "The CASP 16 Experimental
  Protein-Ligand Datasets." *Proteins* 94, 79–85.
  [10.1002/prot.70053](https://doi.org/10.1002/prot.70053). The
  authoritative description of the Roche chymase/catG/ATX datasets and
  the Idorsia Mpro dataset; the source of our per-family context blocks.
  Zenodo data at [zenodo.org/records/16762332](https://zenodo.org/records/16762332).
- **Gilson et al. (2026)** "Assessment of Pharmaceutical Protein–Ligand
  Pose and Affinity Predictions in CASP16." *Proteins* 94, 249–266.
  [10.1002/prot.70061](https://doi.org/10.1002/prot.70061). Describes the
  assessment methodology and per-target-family performance; useful
  background for interpreting CASP16 targets.
- **Alexander, Kryshtafovych et al. (2025)** "Protein Target Highlights
  in CASP16: Insights From the Structure Providers." *Proteins*.
  [10.1002/prot.70025](https://doi.org/10.1002/prot.70025). Covers the
  non-pharma CASP16 targets from the structure-provider perspective.

Related local documentation:

- `~/data/vaults/docs/CASP15-LIGANDS-REPORT.md`
- `~/data/vaults/docs/CASP16-LIGAND-DEEP-DIVE.md`
- `~/data/vaults/docs/CASP16-PAPERS-REPORT.md`
- `~/data/vaults/docs/SYSTEM-casp-pdb-papers.md`
- `~/data/vaults/log/LOG-2026-04-14.md` — session log for this work
