# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A single-script workflow for generating crystallography-paper-style descriptions of protein-ligand binding interactions from PDB files. The example dataset is three TEM-1 beta-lactamase complexes with boronic acid transition-state analogue inhibitors (1ERM/BJI, 1ERO/BJP, 1ERQ/BJH) from Ness et al., *Biochemistry* 2000 (`Ness_TEM_Biochem_2000.pdf`).

There is no build system, no package manifest, and no test suite — just `describe_binding.py` plus input PDBs and generated `.md`/`.json` outputs alongside it.

## Running the pipeline

```bash
# Download a structure
curl -s "https://files.rcsb.org/download/1ERM.pdb" -o 1erm.pdb

# Full LLM narrative (default: OpenAI o3, falls back to Anthropic)
python3 describe_binding.py 1erm.pdb --ligand BJI

# Structured plain text (no LLM call)
python3 describe_binding.py 1erm.pdb --ligand BJI --no-llm

# Raw interaction JSON (feed to other tools or inspection)
python3 describe_binding.py 1erm.pdb --ligand BJI --json > 1erm.json

# Override model
python3 describe_binding.py 1erm.pdb --ligand BJI -m o4-mini
```

`--ligand` takes a three-letter HETATM residue name (auto-detected if omitted). `--ligand-pdb` merges a separate ligand PDB into the protein structure. `--no-waters` disables water-bridge analysis. Narratives are written to `<pdb>.md` by redirecting stdout; JSON interaction dumps live in `<pdb>.json`.

`comparison.md` is a hand-curated cross-structure comparison built by reading the three `*.md` narratives and `*.json` dumps together — regenerate it after rerunning the individual analyses if distances change.

## Dependencies

Uses the system Python 3 with: `biopython` (PDB parsing + `NeighborSearch`), `numpy`, and one of `openai` or `anthropic` for narrative generation. The script also references `rdkit`, `prolif`, and `MDAnalysis` in allowed-commands (see `.claude/settings.local.json`) but does not currently import them — they're available if extending the interaction analysis.

Set `OPENAI_API_KEY` or `ANTHROPIC_API_KEY` before running without `--no-llm`. If both fail the script falls back to the plain-text formatter with a stderr warning.

## Architecture

Everything lives in `describe_binding.py`. The flow is a four-stage pipeline:

1. **Parse & classify atoms** (`find_ligand_residues`, `InteractionAnalyzer.__init__`): walk the Bio.PDB structure, split atoms into ligand / protein / water buckets. Modified residues (HETATM entries with `H_` flag that aren't HOH and aren't the target ligand) are treated as protein — this matters for structures with covalently modified catalytic residues.

2. **Geometric interaction detection** (`InteractionAnalyzer.find_all`): a single `NeighborSearch` over protein+water atoms powers seven distance-based detectors with hardcoded cutoffs at the top of the file (`COVALENT_DIST=2.0`, `HBOND_DIST=3.5`, `SALT_BRIDGE_DIST=4.0`, `HYDROPHOBIC_DIST=4.5`, `PI_STACK_DIST=5.5`, `WATER_BRIDGE_DIST=3.5`). The detectors are purely geometric — no angle filtering on H-bonds, no SMARTS-based aromaticity. Aromatic rings are recognized only on the protein side via the hardcoded `AROMATIC_ATOMS` table; ligand-side pi detection degenerates to "any ligand carbon within 5.5 Å of a protein aromatic centroid." The oxyanion-hole detector is a special case: a ligand O flanked by ≥2 backbone N atoms within H-bond range.

3. **Summarization** (`summarize_interactions`): collapses raw per-atom lists into per-residue groupings suitable for prose generation — H-bonds bucketed by residue, hydrophobic contacts reduced to (residue, closest distance, count), pi interactions deduplicated to one per residue.

4. **Narrative generation** (`generate_narrative`): a long, strict prompt in `describe_binding.py:476-521` instructs the LLM to write crystallography-paper prose with specific style rules (past tense / passive voice, Å symbol with parenthesized distances, banned editorializing vocabulary, no markdown). The prompt is the primary output-quality knob — edit it there rather than post-processing. OpenAI is tried first; Anthropic is the fallback. Reasoning models (`o1`/`o3`/`o4-mini`) use `max_completion_tokens` instead of `max_tokens`.

## Working on this code

- Cutoff constants are the main tuning surface. Changing them shifts what gets reported and therefore what the narrative claims — rerun all three example structures when adjusting.
- The pi-stacking heuristic is deliberately loose; genuine ring-plane geometry requires detecting ligand rings (RDKit would be the natural path).
- `residue_label` hides chain IDs for chain "A" or blank — structures with multiple chains will render ambiguously.
- When regenerating `*.md` files, redirect stdout: `python3 describe_binding.py 1erm.pdb --ligand BJI | tee 1erm.md`.
