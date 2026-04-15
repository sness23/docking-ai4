#!/usr/bin/env python3
"""Build casp16/targets.tsv by joining:
    - The CASP prediction-center archive (protein_aligned.pdb + ligand_*.pdb)
    - The local L#000_exper_affinity.csv files (CASP target -> Compound ID)
    - The Zenodo SI files (Compound ID -> RCSB PDB, series, affinity)

Per row we record: target, family, protein_pdb, ligand_pdb, ligand_resname,
rcsb, series, ic50, covalent, disclosed, note.

Runs idempotently — reads from files, writes targets.tsv.
"""

import csv
import os
import re
from pathlib import Path

VAULT = Path.home() / "data/vaults/casp/casp16"
EXPER = VAULT / "pharma_ligands/exper_struct"
ZENODO = VAULT / "zenodo"
OUT = Path(__file__).resolve().parent / "targets.tsv"

# Pre-disclosed targets (Skrinjar 2026, prot.70053, Section 2.1)
DISCLOSED_CHYMASE = {"L1003", "L1010", "L1013"}
DISCLOSED_ATX = {
    "L3009", "L3047", "L3196", "L3197", "L3211",
    "L3213", "L3214", "L3216", "L3217", "L3219",
    "L3222", "L3224", "L3225", "L3226", "L3229",
}
DISCLOSED_MPRO = {"L4006", "L4007", "L4008", "L4009", "L4010"}

# L4004 is the only cocrystallized covalent inhibitor per Skrinjar Section 3.2.
# (Gilson's deep-dive mentions more but Skrinjar is the authoritative source.)
COVALENT_MPRO = {"L4004"}

FAMILY = {
    "L1000": "chymase",
    "L2000": "cathepsin_g",
    "L3000": "autotaxin",
    "L4000": "mpro",
}

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from describe_binding import SKIP_RESNAMES

# Manual L2000 overrides. Skrinjar's paper (prot.70053, Section 2.2.1)
# lists 7h6h and 7h6g as the two cathepsin G co-crystals; these correspond
# to compound IDs 17 and 16 in chymasecatgSI.txt respectively. The local
# L2001/L2002 .tsv files confirm the match by SMILES.
L2000_OVERRIDES = {
    "L2001": {"compound_id": "16", "rcsb": "7h6g", "series": "5-ring_vinylogous_acids"},
    "L2002": {"compound_id": "17", "rcsb": "7h6h", "series": "indole-carboxylic_acids"},
}


def canonical_smiles(smi):
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        return Chem.MolToSmiles(m, canonical=True, isomericSmiles=False)
    except Exception:
        return None


def load_atx_smiles_map(zen_atx):
    """Build canonical SMILES -> compound_id, for fuzzy ATX joining."""
    smi_to_cid = {}
    zen_raw = ZENODO / "atxSI-2.txt"
    with open(zen_raw) as fh:
        reader = csv.DictReader(fh, delimiter=";")
        for row in reader:
            cid = row["compound id"].strip()
            smi = (row.get("smiles") or "").strip()
            if not smi:
                continue
            canon = canonical_smiles(smi)
            if canon:
                smi_to_cid[canon] = cid
    return smi_to_cid


def load_l3000_tsv_smiles(tgt_set_dir):
    """Read each CASP LNNNN.tsv file to get target -> SMILES."""
    out = {}
    for f in sorted(os.listdir(tgt_set_dir)):
        if not f.endswith(".tsv"):
            continue
        tid = f.replace(".tsv", "")
        with open(tgt_set_dir / f) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                smi = (row.get("SMILES") or "").strip()
                if smi:
                    out[tid] = smi
                    break
    return out


def load_zenodo_chymase():
    """compound_id -> (pdb, series, chymase_ic50_um, catg_ic50_um)"""
    out = {}
    f = ZENODO / "chymasecatgSI.txt"
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter=";")
        for row in reader:
            cid = row["compound id"].strip()
            out[cid] = {
                "pdb": (row.get("pdb") or "").strip(),
                "series": (row.get("series") or "").strip(),
                "chymase_ic50_um": (row.get("chymase IC50 (uM)") or "").strip(),
                "catg_ic50_um": (row.get("cathepsin G IC50 (uM)") or "").strip(),
            }
    return out


def load_zenodo_atx():
    """compound_id -> (pdb, series, atx_ic50_um)"""
    out = {}
    f = ZENODO / "atxSI-2.txt"
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter=";")
        for row in reader:
            cid = row["compound id"].strip()
            out[cid] = {
                "pdb": (row.get("pdb") or "").strip(),
                "series": (row.get("series") or "").strip(),
                "atx_ic50_um": (row.get("ATX  IC50 (uM)") or row.get("ATX IC50 (uM)") or "").strip(),
            }
    return out


def load_zenodo_mpro():
    """compound_id -> pdb. Note: compound IDs are already 4xxx in the SI."""
    out = {}
    f = ZENODO / "3CLMproSI.txt"
    with open(f) as fh:
        reader = csv.DictReader(fh, delimiter=";")
        for row in reader:
            cid = row["compound id"].strip()
            out[cid] = {"pdb": (row.get("pdb") or "").strip()}
    return out


def load_affinity_csv(path):
    """Read an L####_exper_affinity.csv into target_id -> row dict."""
    out = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            # The file has a BOM on the first header
            tid = row.get("Target ID") or row.get("\ufeffTarget ID")
            if tid:
                out[tid.strip()] = row
    return out


def choose_ligand_file(tgt_dir):
    """Pick the drug-like ligand file: first in alphabetical order whose
    compound token is not an ion or crystallization artifact.

    Handles two filename shapes:
      ligand_<compound>_<chain>_<copy>.pdb                (most targets)
      ligand_<compound>_<chain>_<copy>_<altloc>.pdb       (ATX with altlocs)
    """
    cands = []
    for f in sorted(os.listdir(tgt_dir)):
        if not f.startswith("ligand_"):
            continue
        m = re.match(
            r"ligand_(.+?)_([A-Za-z0-9]+)_(\d+)(?:_([A-Z]))?\.pdb",
            f,
        )
        if not m:
            continue
        token = m.group(1)
        if token in SKIP_RESNAMES:
            continue
        cands.append((f, token))
    return cands[0] if cands else (None, None)


def main():
    zen_chy = load_zenodo_chymase()
    zen_atx = load_zenodo_atx()
    zen_mpro = load_zenodo_mpro()

    # ATX SMILES -> compound_id for fuzzy joining
    atx_smi_to_cid = load_atx_smiles_map(zen_atx)
    # CASP L3000 target -> SMILES from per-target TSVs
    l3000_tsvs = VAULT / "pharma_ligands/L3000"
    atx_tgt_to_smi = load_l3000_tsv_smiles(l3000_tsvs) if l3000_tsvs.exists() else {}

    # Local affinity CSVs provide the CASP target -> Compound ID mapping
    aff_l1000 = load_affinity_csv(VAULT / "pharma_ligands/L1000_exper_affinity.csv")
    aff_l3000 = load_affinity_csv(VAULT / "pharma_ligands/L3000_exper_affinity.csv")
    # L2000 has no affinity CSV with the same structure; its 2 targets we map manually.
    # L4000 has no affinity CSV; we use target number directly (L4001 -> compound 4001).

    rows = []

    for tgt_set, family in FAMILY.items():
        root = EXPER / tgt_set / f"{tgt_set}_prepared"
        if not root.is_dir():
            continue
        for tgt in sorted(os.listdir(root)):
            tgt_dir = root / tgt
            if not tgt_dir.is_dir():
                continue
            protein = tgt_dir / "protein_aligned.pdb"
            lig_file, lig_token = choose_ligand_file(tgt_dir)
            if not lig_file:
                rows.append({
                    "target": tgt, "family": family, "status": "NO_LIGAND",
                    "protein_pdb": str(protein), "ligand_pdb": "",
                    "ligand_resname": "", "rcsb": "", "series": "",
                    "ic50_um": "", "covalent": "", "disclosed": "",
                    "note": "no non-artifact ligand file found",
                })
                continue

            row = {
                "target": tgt,
                "family": family,
                "status": "OK",
                "protein_pdb": str(protein),
                "ligand_pdb": str(tgt_dir / lig_file),
                "ligand_resname": lig_token,
                "rcsb": "",
                "series": "",
                "ic50_um": "",
                "covalent": "",
                "disclosed": "",
                "note": "",
            }

            if family == "chymase":
                af = aff_l1000.get(tgt, {})
                cid = (af.get("Compound ID") or "").strip()
                zen = zen_chy.get(cid)
                if zen:
                    row["rcsb"] = zen["pdb"]
                    row["series"] = zen["series"]
                    row["ic50_um"] = zen["chymase_ic50_um"]
                row["disclosed"] = "yes" if tgt in DISCLOSED_CHYMASE else ""
            elif family == "cathepsin_g":
                ov = L2000_OVERRIDES.get(tgt)
                if ov:
                    zen = zen_chy.get(ov["compound_id"])
                    row["rcsb"] = ov["rcsb"]
                    row["series"] = ov["series"]
                    if zen:
                        row["ic50_um"] = zen["catg_ic50_um"]
            elif family == "autotaxin":
                # Primary path: local CSV -> compound id -> Zenodo row.
                af = aff_l3000.get(tgt, {})
                cid = (af.get("Compound ID") or "").strip()
                zen = zen_atx.get(cid) if cid else None
                # Fallback path: canonicalize the CASP target's SMILES and
                # look it up in the Zenodo atxSI. Catches the pose-only
                # targets that don't appear in the affinity CSV.
                if not zen:
                    smi = atx_tgt_to_smi.get(tgt, "")
                    canon = canonical_smiles(smi)
                    if canon and canon in atx_smi_to_cid:
                        cid = atx_smi_to_cid[canon]
                        zen = zen_atx.get(cid)
                if zen:
                    row["rcsb"] = zen["pdb"]
                    row["series"] = zen["series"]
                    row["ic50_um"] = zen["atx_ic50_um"]
                row["disclosed"] = "yes" if tgt in DISCLOSED_ATX else ""
            elif family == "mpro":
                # L4001 -> compound 4001, etc. No gap in Zenodo SI but
                # CASP target numbering skips 4005 and 4012.
                n = tgt.lstrip("L")
                zen = zen_mpro.get(n)
                if zen:
                    row["rcsb"] = zen["pdb"]
                row["covalent"] = "yes" if tgt in COVALENT_MPRO else ""
                row["disclosed"] = "yes (pre-CASP)" if tgt in DISCLOSED_MPRO else ""

            rows.append(row)

    # Write TSV
    cols = ["target", "family", "status", "rcsb", "series", "ic50_um",
            "covalent", "disclosed", "ligand_resname", "protein_pdb",
            "ligand_pdb", "note"]
    with open(OUT, "w") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in cols})

    # Summary
    from collections import Counter
    by_fam = Counter(r["family"] for r in rows)
    with_rcsb = Counter(r["family"] for r in rows if r["rcsb"])
    with_series = Counter(r["family"] for r in rows if r["series"])
    disclosed = Counter(r["family"] for r in rows if r["disclosed"])
    covalent = Counter(r["family"] for r in rows if r["covalent"])
    print(f"Wrote {len(rows)} rows to {OUT}")
    print(f"  by family: {dict(by_fam)}")
    print(f"  with RCSB code: {dict(with_rcsb)}")
    print(f"  with series label: {dict(with_series)}")
    print(f"  disclosed: {dict(disclosed)}")
    print(f"  covalent: {dict(covalent)}")


if __name__ == "__main__":
    main()
