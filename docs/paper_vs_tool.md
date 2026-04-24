# Comparison: Computational Analysis vs. Published Paper

Ness et al. (2000) *Biochemistry* 39, 5312-5321 (PDB: 1ERM, 1ERO, 1ERQ)

This document compares the output of `describe_binding.py` (our computational tool) against the descriptions and Table 2 in the original published paper. The goal is to identify what our tool gets right, what it misses, and where there are discrepancies.

---

## 1. Compound Numbering Mismatch

The paper and our `comparison.md` use different compound labels:

| Paper name | PDB | Ligand ID | Our label (comparison.md) |
|------------|-----|-----------|---------------------------|
| Compound 1 | 1ERO | BJP | Compound 1 (correct) |
| Compound 2 | 1ERQ | BJH | Compound 2f (should be 2) |
| Compound 3 | 1ERM | BJI | Compound 2 (should be 3) |

Our comparison.md incorrectly labeled 1ERM as "compound 2" when the paper calls it "compound 3" (the original N-acetyl inhibitor from the 1996 Strynadka *Nat. Struct. Biol.* paper). And we called 1ERQ "compound 2f" when the paper calls it "compound 2". This is a labeling error in our output that could confuse readers familiar with the paper.

---

## 2. Residue Numbering: Arg244 vs. Arg243

The paper consistently refers to **Arg244**, but the PDB file numbers this residue as **Arg243**. We verified that PDB residue 243 is ARG and residue 244 is GLY. The PDB DBREF record shows residues 26-288 mapped to UniProt 24-286, indicating a possible offset in Ambler numbering. The paper appears to use a numbering convention with a +1 offset relative to the deposited PDB coordinates for this residue.

Our tool correctly identifies the interacting arginine from the coordinates but calls it Arg243 (from the PDB), while the paper calls it Arg244. This could be confusing when cross-referencing.

**Impact**: A student comparing our output to the paper would see different residue numbers for the same interaction.

---

## 3. Atom Naming Conventions

The paper uses different names for the boronate oxygen atoms than the PDB:

| Paper name | PDB name | Role |
|------------|----------|------|
| OH1 | OB1 | Oxyanion hole oxygen |
| OH2 | OB2 | Deacylating water mimic |
| OH3 | OB3 / O5 | Hydroxyl on phenyl (compound 2 only) |

The paper also uses some non-standard enzyme atom names in Table 2 (e.g., "Glu166 OH1" appears to mean Glu166 OE1). Our tool uses the PDB atom names directly, which are the standard and correct names, but this creates a translation barrier when comparing to the paper.

---

## 4. Distance Comparison: Table 2 vs. Computed Values

Here we compare the paper's Table 2 distances with our computed distances from the PDB coordinates. The paper rounds to one decimal place.

### Carboxylate interactions

| Interaction | Paper compd 1 (1ERO) | Computed 1ERO | Paper compd 2 (1ERQ) | Computed 1ERQ | Paper compd 3 (1ERM) | Computed 1ERM |
|-------------|---------------------|---------------|---------------------|---------------|---------------------|---------------|
| O3 -- Ser130 OG | 3.5 | 3.52 | 3.4 | 3.43 | 2.5 | 2.49 |
| O3 -- Ser235 OG | * | 2.29 | * | 2.76 | * | 2.85 |
| O4 -- Ser235 OG | * | 3.03 | * | 3.32 | * | 2.83 |
| O4 -- Arg243 NH1 | * | 3.12 | 2.8 | 2.76 | 2.8 | 2.81 |

\* = value in Table 2 but hard to read precisely from PDF

**Verdict**: The clearly readable values (O3--Ser130 OG, O4--Arg243 NH1) match well between our computed distances and the paper. Differences are within rounding.

### Oxyanion hole

| Interaction | Paper compd 1 | Computed 1ERO | Paper compd 2 | Computed 1ERQ | Paper compd 3 | Computed 1ERM |
|-------------|--------------|---------------|--------------|---------------|--------------|---------------|
| OH1/OB1 -- Ser70 N | 2.9 | 2.85 | 2.8 | 2.80 | * | 2.73 |
| OH1/OB1 -- Ala237 N | * | 3.42 | * | 3.07 | * | 3.01 |

**Verdict**: Good agreement. The Ser70 N distances match the paper's rounded values.

### Glu166 / deacylating water position

| Interaction | Paper compd 1 | Computed 1ERO | Paper compd 2 | Computed 1ERQ | Paper compd 3 | Computed 1ERM |
|-------------|--------------|---------------|--------------|---------------|--------------|---------------|
| OH2/OB2 -- Glu166 OE1 | 2.9 | 2.57 | 2.5 | 2.87 | * | 2.46 |

**Discrepancy**: The paper reports OH2--Glu166 as **2.9 A for compound 1** and **2.5 A for compound 2**. Our computed values are **2.57 for 1ERO (compd 1)** and **2.87 for 1ERQ (compd 2)**. These don't match the paper's values well. Possible explanations:
- The paper may be referring to a different Glu166 atom (OE2 rather than OE1)
- The PDB coordinates may have been re-refined or updated since publication
- The paper's "Glu166 OH1" label in Table 2 may not correspond to OE1

Our computed OB2--Glu166 OE2 distances are: 1ERO = 3.37, 1ERQ = 3.73, 1ERM = 3.25. Neither OE1 nor OE2 consistently matches the Table 2 values, suggesting possible coordinate differences between the published and deposited structures.

### Hydroxyl-specific interactions (compound 2 / 1ERQ only)

| Interaction | Paper | Computed |
|-------------|-------|----------|
| OH3/O5 -- Ser130 OG | 2.9 | 2.73 |
| OH3/O5 -- Ser70 OG | 2.5 (or 2.3) | 2.30 |

**Verdict**: Reasonable agreement. The Ser70 OG distance matches the paper's value of 2.3 exactly. The Ser130 OG distance (2.73 vs 2.9) shows a 0.17 A difference.

---

## 5. Interactions Described in the Paper That Our Tool Misses

### 5a. Lys73--Glu166 hydrogen bond (protein-protein interaction)

This is arguably the most important omission. The paper highlights a hydrogen bond between **Lys73 NZ and Glu166** that forms upon inhibitor binding:
- Compound 3 (1ERM): 2.8 A (paper)
- Compound 2 (1ERQ): 2.7 A (paper)

The paper interprets this as evidence for the role of Lys73 in regenerating the Ser70 nucleophile via proton transfer through the general base Glu166.

Our computed distances from the current PDB coordinates:
- 1ERM: Lys73 NZ -- Glu166 OE1 = **3.18 A**
- 1ERQ: Lys73 NZ -- Glu166 OE1 = **3.57 A**
- 1ERO: Lys73 NZ -- Glu166 OE1 = **3.39 A**

These are substantially longer than the paper's reported values (3.18 vs 2.8, 3.57 vs 2.7). This is a significant discrepancy that likely reflects coordinate differences between the originally published and the currently deposited (possibly re-refined) PDB structures.

**Why our tool misses this**: Our tool only detects protein-ligand interactions, not protein-protein interactions that change upon ligand binding. Even if it did, the current coordinates don't show the short distances the paper reports.

### 5b. Inhibitor occupancy and alternate conformations (compound 1 / 1ERO)

The paper states that compound 1 was refined to **0.50 occupancy** because the phenyl group makes it more hydrophobic and less soluble in crystallization conditions. PDB confirms: the core atoms of BJP are at 0.50 occupancy, and the phenylacetamido group (C19-C26, O6, N1) is split into two alternate conformations at **0.25 occupancy each**.

Our tool ignores occupancy entirely. All distances for 1ERO are computed using half-occupied atoms without any caveat about reliability. The paper explicitly warns that distances from compound 1 should be interpreted with caution.

### 5c. Water displacement

The paper devotes an entire subsection ("Displacement of Ordered Water To Facilitate Increased Binding") to the entropic contribution of water displacement. Key findings:
- Compounds 2 (1ERQ) and 3 (1ERM) retain several conserved ordered water molecules
- Compound 1 (1ERO) displaces at least 3 more ordered water molecules due to its bulkier phenylacetamido group

Our computed data confirms this trend:
- 1ERM: 5 waters within 4 A of ligand, 137 total
- 1ERQ: 4 waters within 4 A of ligand, 121 total
- 1ERO: **0 waters** within 4 A of ligand, **63 total**

Our tool detects water-mediated bridges but does not analyze or report water displacement. For structure-based design, understanding which waters are displaced and which are retained is essential.

### 5d. 25-degree ring rotation

The paper describes a 25-degree rotation of the hydroxyphenyl ring in compound 2 (1ERQ) relative to the carboxyphenyl ring in compound 3 (1ERM). This unanticipated rotation brings the hydroxyl oxygen (OH3/O5) into hydrogen bonding position with Ser70 OG and creates additional contacts.

Our tool does not measure or report ligand conformation differences between structures.

### 5e. Comparison to penicillin G substrate binding

The paper extensively compares the inhibitor binding mode to the natural substrate penicillin G (PDB 1FQG), showing that the inhibitors mimic the substrate's carboxylate position and the C7 carbonyl of the beta-lactam ring. Our tool only analyzes the input PDB and cannot make cross-structure comparisons.

### 5f. Conformational changes in the protein

The paper discusses conformational changes in active site residues upon inhibitor binding compared to the native enzyme. These include side chain rearrangements of Tyr105, Glu166, and others. Our tool uses a static single structure and cannot detect ligand-induced conformational changes.

### 5g. Glu166 to Ser70 distance

The paper notes that "The distance of Glu166 to Ser70 in both the native and inhibitor bound structures is more than 4 A, and thus, barring significant conformational changes, a direct transfer between the two residues seems unlikely." This observation about the Glu166--Ser70 distance being too long for direct proton transfer is mechanistically important.

Our computed Glu166 OE1 -- Ser70 OG distances: 1ERM = 3.49, 1ERO = 3.61, 1ERQ = 4.16. These are consistent with the paper's statement for 1ERQ but shorter for 1ERM/1ERO than the ">4 A" the paper reports, again suggesting possible coordinate differences.

---

## 6. Interactions Our Tool Reports That the Paper Does Not Emphasize

### 6a. Asn170 OD1 contacts

Our tool reports OB2--Asn170 OD1 contacts at 2.56-2.76 A in all three structures. The paper mentions Asn170 in Table 2 (as "Asn170 N" = probably ND2 or OD1) but only lists values for compounds 1 and 3 (2.8 A each). The paper focuses less on the Asn170 interaction than our tool does.

### 6b. N1--Asn132 ND2 distance

The paper's Table 2 lists an N1--Asn132 ND2 interaction at **2.9 A** for compounds 2 and 3. However, our computed distances are:
- 1ERM (compd 3): N1--Asn132 ND2 = **4.40 A**
- 1ERQ (compd 2): N1--Asn132 ND2 = **4.44 A**
- 1ERO (compd 1): N1--Asn132 ND2 = **4.36 A**

This is a major discrepancy. At 4.4 A, this is not a hydrogen bond. The paper's reported 2.9 A contact cannot be reproduced from the current PDB coordinates. This strongly suggests the deposited coordinates have been re-refined since publication, or there is a systematic difference in how the interaction was measured.

Alternatively, the paper may be describing the interaction from the amide carbonyl O2 to Asn132 ND2 rather than from N1:
- 1ERM: O2--Asn132 ND2 = **2.85 A** (matches 2.9 well)
- 1ERQ: O2--Asn132 ND2 = **2.86 A** (matches 2.9 well)

This seems likely — the paper's Table 2 labels this as "amide / N1 / Asn132 ND2" but the actual measured distance may be from the carbonyl oxygen of the amide group, not the nitrogen.

---

## 7. Summary of Discrepancies

| Category | Severity | Description |
|----------|----------|-------------|
| Compound numbering | Medium | Our labels don't match the paper's numbering scheme |
| Arg244 vs Arg243 | Low | Numbering offset between paper and PDB |
| Atom naming (OH1/OB1) | Low | Different conventions, same atoms |
| Lys73-Glu166 H-bond | **High** | Key mechanistic interaction completely missed (protein-protein, not protein-ligand) |
| Occupancy/alt confs | **High** | Compound 1 at 0.50 occupancy with alternate conformations — distances unreliable, not flagged |
| Water displacement | **High** | Paper's central argument about entropic contribution not captured |
| Glu166 distances | Medium | Computed distances don't match Table 2 for OB2-Glu166 |
| N1-Asn132 distance | Medium | Paper reports 2.9 A, we compute 4.4 A — likely a mislabel in paper (should be O2, not N1) |
| Ring rotation | Medium | 25-degree conformational difference between compounds not measured |
| Substrate comparison | Low | Cannot compare to penicillin G binding mode |
| Protein conformational changes | Medium | Cannot detect ligand-induced rearrangements |

---

## 8. Recommendations for Improving the Tool

Based on these discrepancies, the following enhancements would bring the computational analysis closer to paper quality:

1. **Add protein-protein interaction detection near the ligand**: Identify hydrogen bonds between protein residues within ~6 A of the ligand. The Lys73-Glu166 interaction is the most important mechanistic observation in the paper, and our tool misses it entirely.

2. **Report and warn about occupancy**: Flag ligand atoms with occupancy < 1.0 and warn that distances may be unreliable.

3. **Analyze water displacement**: Compare ordered waters in the complex to a reference native structure, identifying which conserved waters are displaced by the ligand.

4. **Support multi-structure comparison**: Allow comparing two or more complexes to identify conformational differences (like the 25-degree ring rotation).

5. **Map residue numbering**: Add an option to convert PDB residue numbers to Ambler (or other standard) numbering schemes using the DBREF record.

6. **Translate atom names**: Provide a mapping between PDB boronate atom names (OB1/OB2) and common literature names (OH1/OH2) in the output.
