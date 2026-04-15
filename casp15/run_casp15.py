#!/usr/bin/env python3
"""Run describe_binding.py across all runnable CASP15 protein-ligand targets.

Reads targets.tsv, skips rows marked "skip", and for each remaining target
emits:
    casp15/out/<target>.json — raw interaction data
    casp15/out/<target>.txt  — --no-llm structured text
    casp15/out/<target>.md   — LLM narrative (only if --with-llm)

Usage:
    python3 casp15/run_casp15.py           # geometry only, no LLM
    python3 casp15/run_casp15.py --with-llm   # also call the LLM
    python3 casp15/run_casp15.py --only T1124,T1188  # subset
"""

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "describe_binding.py"
TARGETS_TSV = REPO / "casp15" / "targets.tsv"
LIG_DIR = Path.home() / "data/vaults/casp/casp15_ligands/targets_ligand/Targets_ligand"
OUT_DIR = REPO / "casp15" / "out"


def load_targets():
    with open(TARGETS_TSV) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return [row for row in reader]


def run_one(row, with_llm, model):
    target = row["target"]
    lig_pdb = LIG_DIR / row["lig_pdb"]
    resname = row["ligand_resname"]
    out_json = OUT_DIR / f"{target}.json"
    out_txt = OUT_DIR / f"{target}.txt"
    out_md = OUT_DIR / f"{target}.md"

    if not lig_pdb.exists():
        return target, "MISSING_FILE", f"{lig_pdb} not found"

    base = ["python3", str(SCRIPT), str(lig_pdb), "--ligand", resname]

    # JSON pass
    with open(out_json, "w") as f:
        r = subprocess.run(base + ["--json"], stdout=f, stderr=subprocess.PIPE)
    if r.returncode != 0:
        return target, "JSON_FAIL", r.stderr.decode()[:300]

    # No-LLM structured text
    with open(out_txt, "w") as f:
        r = subprocess.run(base + ["--no-llm"], stdout=f, stderr=subprocess.PIPE)
    if r.returncode != 0:
        return target, "TEXT_FAIL", r.stderr.decode()[:300]

    if with_llm:
        with open(out_md, "w") as f:
            cmd = base + ["-m", model]
            r = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)
        if r.returncode != 0:
            return target, "LLM_FAIL", r.stderr.decode()[:300]
        return target, "OK+LLM", ""

    return target, "OK", ""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--with-llm", action="store_true",
                        help="Also generate LLM narratives (costs API calls).")
    parser.add_argument("-m", "--model", default="o4-mini",
                        help="LLM model to pass to describe_binding.py (default: o4-mini)")
    parser.add_argument("--only", default="",
                        help="Comma-separated target list to restrict run")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    targets = load_targets()
    runnable = [r for r in targets if r["family"] != "skip"]
    if args.only:
        wanted = set(s.strip() for s in args.only.split(","))
        runnable = [r for r in runnable if r["target"] in wanted]

    print(f"Running {len(runnable)} targets "
          f"(LLM={'yes' if args.with_llm else 'no'}, model={args.model})",
          file=sys.stderr)

    results = []
    for row in runnable:
        target, status, err = run_one(row, args.with_llm, args.model)
        results.append((target, status, err))
        marker = "✓" if status.startswith("OK") else "✗"
        print(f"  {marker} {target:<10} {status}"
              f"{'  ' + err[:80] if err else ''}",
              file=sys.stderr)

    fail = [r for r in results if not r[1].startswith("OK")]
    print(f"\nDone. {len(results) - len(fail)}/{len(results)} succeeded.",
          file=sys.stderr)
    return 0 if not fail else 1


if __name__ == "__main__":
    sys.exit(main())
