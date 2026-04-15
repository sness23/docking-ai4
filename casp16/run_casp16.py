#!/usr/bin/env python3
"""Run describe_binding.py across CASP16 pharma ligand targets.

Reads casp16/targets.tsv (built by build_targets.py) and runs the analysis
in split-protein/ligand mode:

    python3 describe_binding.py <protein_aligned.pdb> \
        --ligand-pdb <ligand.pdb> --ligand <resname> \
        --context <per-family context file>

Outputs go to casp16/out/<target>.{json,txt,md} — same layout as CASP15.

Usage:
    python3 casp16/run_casp16.py --only-family chymase           # geometry only
    python3 casp16/run_casp16.py --only-family chymase --with-llm
    python3 casp16/run_casp16.py --only L1001,L4004 --with-llm
    python3 casp16/run_casp16.py --skip-family autotaxin         # everything but ATX
"""

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "describe_binding.py"
TARGETS_TSV = REPO / "casp16" / "targets.tsv"
OUT_DIR = REPO / "casp16" / "out"
CTX_DIR = REPO / "casp16" / "context"


FAMILY_CONTEXT = {
    "chymase": """\
Human chymase (EC 3.4.21.39) is a chymotryptic serine protease (S1 family)
stored in mast cell granules. The catalytic triad is Ser203, His66, Asp113 in
UniProt numbering (note: the PDB sequence numbering used in the coordinate
file is shifted relative to UniProt; favor whatever residue numbers appear in
the coordinate data). The S1 pocket surface is electropositive — lined by
Lys49, His66, Arg151, and Lys200 — which favors acidic ligand motifs and is
responsible for the acidic group present in every Roche chymase series.

Four chemical series are represented in this dataset: indole-carboxylic
acids, 4-ring vinylogous acids, 5-ring vinylogous acids, and phenylsulfonamides.
The phenylsulfonamide class is mechanistically distinct — it binds via
induced fit to a previously undisclosed sub-pocket that opens when the Leu107
side chain moves aside, and — unlike the other three series — does NOT
interact directly with the Ser203 side chain. If the interaction data shows
no Ser203 side-chain contact, do not try to invent one.
""",

    "cathepsin_g": """\
Human cathepsin G (EC 3.4.21.20) is a neutrophil chymotryptic serine protease
of the S1 family, structurally similar to chymase. The catalytic triad is
Ser195, His57, Asp102 in classical chymotrypsin numbering. Cathepsin G has
unusual dual substrate specificity because Glu226 at the bottom of the S1
pocket divides it into two compartments, accommodating either positively
charged (Lys, Arg) or hydrophobic (Phe, Tyr) P1 residues. Only two co-crystal
structures are provided in the CASP16 dataset — both were obtained during a
chymase drug-discovery project as off-target selectivity controls. No affinity
predictions are required for this target.
""",

    "autotaxin": """\
Human autotaxin (ATX, ENPP2, EC 3.1.4.39) is a secreted nucleotide
pyrophosphatase/phosphodiesterase that hydrolyzes lysophosphatidylcholine to
generate lysophosphatidic acid (LPA). The active site contains TWO catalytic
Zn²⁺ ions, plus a hydrophobic pocket and a lipid channel. The enzyme carries
an N-glycosylation on Asn497 that sits far from the active site.

Two inhibitor classes dominate the CASP16 dataset: pyrrolo-pyrroles contain a
zinc-binding sulfonamide group and directly coordinate the two catalytic
zinc ions; phenoxymethyl-imidazoles occupy the hydrophobic pocket and channel
but do NOT coordinate the zinc ions. When describing binding geometry, be
explicit about whether the ligand coordinates the zinc ions if the
interaction data is informative on that point.
""",

    "mpro": """\
SARS-CoV-2 main protease (Mpro, 3CLpro, NSP5, EC 3.4.22.69) is a cysteine
protease essential for coronavirus polyprotein processing. It functions as a
homodimer (306 residues per monomer); the two active sites are symmetric and
one is sufficient to describe. The catalytic dyad is Cys145 and His41. The
S1 sub-pocket prefers Gln at P1; the S2 pocket is hydrophobic and leucine-
shaped; the S4 pocket is solvent-exposed.

The Idorsia CASP16 dataset was obtained in the P2₁2₁2₁ crystal form, which
has a somewhat more open and less constrained active-site conformation than
the typical C2 crystal form used in most published Mpro inhibitor work.
""",
}


# Per-target refinements layered on top of the family context.
TARGET_NOTES = {
    "L4004": "This specific target (L4004) is the only cocrystallized covalent "
             "inhibitor in the Idorsia set — a covalent linkage to the "
             "catalytic Cys145 is expected and should be described if the "
             "geometric data shows a sub-2 Å contact there.",
}


def write_context_files():
    CTX_DIR.mkdir(parents=True, exist_ok=True)
    for fam, text in FAMILY_CONTEXT.items():
        (CTX_DIR / f"{fam}.txt").write_text(text)


def load_targets():
    with open(TARGETS_TSV) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def build_context(row):
    """Family context + per-target notes + per-target metadata (series,
    RCSB code, disclosed flag) as a single text blob."""
    fam = row["family"]
    parts = [FAMILY_CONTEXT.get(fam, "")]
    meta = []
    if row.get("rcsb"):
        meta.append(f"Reference PDB: {row['rcsb']}")
    if row.get("series"):
        meta.append(f"Chemical series: {row['series']}")
    if row.get("ic50_um"):
        meta.append(f"IC50: {row['ic50_um']} μM")
    if row.get("covalent") == "yes":
        meta.append("Covalent inhibitor (bound via covalent linkage).")
    if row.get("disclosed"):
        meta.append(f"Disclosure status: {row['disclosed']}")
    if meta:
        parts.append("Target-specific metadata: " + "; ".join(meta) + ".")
    note = TARGET_NOTES.get(row["target"])
    if note:
        parts.append(note)
    return "\n\n".join(p for p in parts if p)


def run_one(row, with_llm, model):
    target = row["target"]
    protein = row["protein_pdb"]
    lig_pdb = row["ligand_pdb"]
    resname = row["ligand_resname"]

    if not protein or not os.path.exists(protein):
        return target, "MISSING_PROTEIN", protein
    if not lig_pdb or not os.path.exists(lig_pdb):
        return target, "MISSING_LIGAND", lig_pdb

    out_json = OUT_DIR / f"{target}.json"
    out_txt = OUT_DIR / f"{target}.txt"
    out_md = OUT_DIR / f"{target}.md"

    base = ["python3", str(SCRIPT), protein,
            "--ligand-pdb", lig_pdb, "--ligand", resname]

    # JSON pass
    with open(out_json, "w") as f:
        r = subprocess.run(base + ["--json"], stdout=f, stderr=subprocess.PIPE)
    if r.returncode != 0:
        return target, "JSON_FAIL", r.stderr.decode()[:300]

    # Plain-text pass
    with open(out_txt, "w") as f:
        r = subprocess.run(base + ["--no-llm"], stdout=f, stderr=subprocess.PIPE)
    if r.returncode != 0:
        return target, "TEXT_FAIL", r.stderr.decode()[:300]

    if with_llm:
        ctx_text = build_context(row)
        ctx_file = OUT_DIR / f"{target}.ctx.txt"
        ctx_file.write_text(ctx_text)
        with open(out_md, "w") as f:
            cmd = base + ["-m", model, "--context", str(ctx_file)]
            r = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)
        if r.returncode != 0:
            return target, "LLM_FAIL", r.stderr.decode()[:300]
        return target, "OK+LLM", ""

    return target, "OK", ""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--with-llm", action="store_true")
    parser.add_argument("-m", "--model", default="o4-mini")
    parser.add_argument("--only-family", default="",
                        help="Restrict to one family (chymase, cathepsin_g, autotaxin, mpro)")
    parser.add_argument("--skip-family", default="",
                        help="Exclude a family")
    parser.add_argument("--only", default="",
                        help="Comma-separated target list")
    parser.add_argument("--limit", type=int, default=0,
                        help="Stop after N successful targets (0 = no limit)")
    args = parser.parse_args()

    write_context_files()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    targets = load_targets()
    if args.only_family:
        targets = [r for r in targets if r["family"] == args.only_family]
    if args.skip_family:
        targets = [r for r in targets if r["family"] != args.skip_family]
    if args.only:
        wanted = set(s.strip() for s in args.only.split(","))
        targets = [r for r in targets if r["target"] in wanted]

    print(f"Running {len(targets)} targets (LLM={'yes' if args.with_llm else 'no'}, "
          f"model={args.model})", file=sys.stderr)

    results = []
    done_ok = 0
    for row in targets:
        target, status, err = run_one(row, args.with_llm, args.model)
        results.append((target, status, err))
        marker = "✓" if status.startswith("OK") else "✗"
        print(f"  {marker} {target:<8} {row['family']:<12} {status}"
              f"{'  ' + err[:60] if err else ''}", file=sys.stderr)
        if status.startswith("OK"):
            done_ok += 1
            if args.limit and done_ok >= args.limit:
                break

    fail = [r for r in results if not r[1].startswith("OK")]
    print(f"\nDone. {len(results) - len(fail)}/{len(results)} succeeded.",
          file=sys.stderr)
    return 0 if not fail else 1


if __name__ == "__main__":
    sys.exit(main())
