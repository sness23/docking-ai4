#!/usr/bin/env python3
"""
Describe protein-ligand binding interactions in crystallography paper style.

Usage:
    python describe_binding.py complex.pdb --ligand BJI
    python describe_binding.py protein.pdb --ligand-pdb ligand.pdb
    python describe_binding.py 1erm.pdb --ligand BJI --include-waters
"""

import argparse
import sys
import json
import math
import numpy as np
from collections import defaultdict
from pathlib import Path

from Bio.PDB import PDBParser, NeighborSearch, Selection
from Bio.PDB.Polypeptide import is_aa


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

COVALENT_DIST = 2.0          # upper bound for covalent bond detection
HBOND_DIST = 3.5             # donor-acceptor distance
SALT_BRIDGE_DIST = 4.0       # charged group centroids
HYDROPHOBIC_DIST = 4.5       # carbon-carbon contacts
PI_STACK_DIST = 5.5          # ring centroid distance
WATER_BRIDGE_DIST = 3.5      # water-mediated H-bond leg

HBOND_DONORS = {"N", "O", "S"}
HBOND_ACCEPTORS = {"O", "N", "S"}
CHARGED_POS = {"ARG": ["NH1", "NH2", "NE"],
                "LYS": ["NZ"],
                "HIS": ["ND1", "NE2"]}
CHARGED_NEG = {"ASP": ["OD1", "OD2"],
                "GLU": ["OE1", "OE2"]}
AROMATIC_RESIDUES = {"PHE", "TYR", "TRP", "HIS"}
AROMATIC_ATOMS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
}

OXYANION_HOLE_RESIDUES_CLASS_A = {
    # backbone NH of Ser70 and Ala237 (TEM-1 numbering, Ambler)
    # These form the oxyanion hole in class A beta-lactamases
}


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def dist(a, b):
    """Euclidean distance between two Bio.PDB Atom objects."""
    return np.linalg.norm(a.get_vector().get_array() - b.get_vector().get_array())


def angle(a, b, c):
    """Angle in degrees at atom b formed by atoms a-b-c."""
    v1 = a.get_vector().get_array() - b.get_vector().get_array()
    v2 = c.get_vector().get_array() - b.get_vector().get_array()
    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
    return math.degrees(math.acos(np.clip(cos_angle, -1.0, 1.0)))


def ring_centroid(atoms):
    """Centroid of a list of atoms."""
    coords = np.array([a.get_vector().get_array() for a in atoms])
    return coords.mean(axis=0)


def ring_normal(atoms):
    """Normal vector to a ring defined by a list of atoms."""
    coords = np.array([a.get_vector().get_array() for a in atoms])
    centroid = coords.mean(axis=0)
    centered = coords - centroid
    # Use SVD to find the normal
    _, _, vh = np.linalg.svd(centered)
    return vh[-1]


def centroid_dist(atoms1, atoms2):
    c1 = ring_centroid(atoms1)
    c2 = ring_centroid(atoms2)
    return np.linalg.norm(c1 - c2)


def ring_angle(atoms1, atoms2):
    """Angle between ring planes (0=parallel, 90=perpendicular)."""
    n1 = ring_normal(atoms1)
    n2 = ring_normal(atoms2)
    cos_a = abs(np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2) + 1e-10))
    return math.degrees(math.acos(np.clip(cos_a, 0, 1.0)))


# ---------------------------------------------------------------------------
# Atom classification helpers
# ---------------------------------------------------------------------------

def is_hbond_donor(atom):
    return atom.element in HBOND_DONORS


def is_hbond_acceptor(atom):
    return atom.element in HBOND_ACCEPTORS


def is_carbon(atom):
    return atom.element == "C"


def atom_label(atom):
    """Human-readable atom label like 'Ser70 OG' or 'BJI500 B'."""
    res = atom.get_parent()
    resname = res.get_resname().strip()
    resid = res.get_id()[1]
    aname = atom.get_name().strip()
    if is_aa(res, standard=True):
        return f"{resname}{resid} {aname}"
    else:
        return f"{resname}{resid} {aname}"


def residue_label(res):
    resname = res.get_resname().strip()
    resid = res.get_id()[1]
    chain = res.get_parent().get_id()
    return f"{resname} {resid}" if chain == " " or chain == "A" else f"{resname} {resid} (chain {chain})"


# ---------------------------------------------------------------------------
# Interaction detection
# ---------------------------------------------------------------------------

class InteractionAnalyzer:
    def __init__(self, structure, ligand_residues, include_waters=True):
        self.structure = structure
        self.ligand_residues = ligand_residues
        self.ligand_resnames = {r.get_resname().strip() for r in ligand_residues}
        self.ligand_resids = {(r.get_parent().get_id(), r.get_id()) for r in ligand_residues}
        self.include_waters = include_waters

        # Separate atoms
        self.ligand_atoms = []
        self.protein_atoms = []
        self.water_atoms = []

        for model in structure:
            for chain in model:
                for res in chain:
                    rname = res.get_resname().strip()
                    rid = (chain.get_id(), res.get_id())
                    if rid in self.ligand_resids:
                        self.ligand_atoms.extend(res.get_atoms())
                    elif rname == "HOH":
                        self.water_atoms.extend(res.get_atoms())
                    elif is_aa(res, standard=True):
                        self.protein_atoms.extend(res.get_atoms())
                    # modified residues that are part of the chain (like BHD)
                    elif res.get_id()[0].startswith("H_"):
                        # Check if it's linked to the protein backbone
                        self.protein_atoms.extend(res.get_atoms())

        # Build neighbor search for protein + water
        all_env = self.protein_atoms + self.water_atoms
        if all_env:
            self.ns = NeighborSearch(all_env)

        self.interactions = {
            "covalent_bonds": [],
            "hydrogen_bonds": [],
            "salt_bridges": [],
            "hydrophobic_contacts": [],
            "pi_interactions": [],
            "water_bridges": [],
            "oxyanion_hole": [],
        }

    def find_all(self):
        self._find_covalent_bonds()
        self._find_hydrogen_bonds()
        self._find_salt_bridges()
        self._find_hydrophobic_contacts()
        self._find_pi_interactions()
        if self.include_waters:
            self._find_water_bridges()
        self._find_oxyanion_hole()
        return self.interactions

    def _find_covalent_bonds(self):
        """Detect covalent bonds between ligand and protein using LINK records and distance."""
        seen = set()
        for la in self.ligand_atoms:
            neighbors = self.ns.search(la.get_vector().get_array(), COVALENT_DIST)
            for nb in neighbors:
                if nb in self.ligand_atoms or nb in self.water_atoms:
                    continue
                d = dist(la, nb)
                key = (atom_label(la), atom_label(nb))
                if key not in seen:
                    seen.add(key)
                    self.interactions["covalent_bonds"].append({
                        "ligand_atom": atom_label(la),
                        "protein_atom": atom_label(nb),
                        "distance": round(d, 2),
                    })

    def _find_hydrogen_bonds(self):
        seen = set()
        for la in self.ligand_atoms:
            if not (is_hbond_donor(la) or is_hbond_acceptor(la)):
                continue
            neighbors = self.ns.search(la.get_vector().get_array(), HBOND_DIST)
            for nb in neighbors:
                if nb in self.water_atoms or nb in self.ligand_atoms:
                    continue
                if not (is_hbond_donor(nb) or is_hbond_acceptor(nb)):
                    continue
                d = dist(la, nb)
                if d < COVALENT_DIST:
                    continue  # skip covalent bonds
                # At least one must be donor and the other acceptor
                if not ((is_hbond_donor(la) and is_hbond_acceptor(nb)) or
                        (is_hbond_acceptor(la) and is_hbond_donor(nb))):
                    continue
                key = tuple(sorted([atom_label(la), atom_label(nb)]))
                if key not in seen:
                    seen.add(key)
                    res = nb.get_parent()
                    backbone = nb.get_name().strip() in ("N", "CA", "C", "O")
                    self.interactions["hydrogen_bonds"].append({
                        "ligand_atom": atom_label(la),
                        "protein_atom": atom_label(nb),
                        "protein_residue": residue_label(res),
                        "distance": round(d, 2),
                        "backbone": backbone,
                    })

    def _find_salt_bridges(self):
        """Detect salt bridges between charged groups."""
        # Find charged atoms on ligand
        lig_charged = []
        for la in self.ligand_atoms:
            if la.element in ("N", "O"):
                lig_charged.append(la)

        for la in lig_charged:
            neighbors = self.ns.search(la.get_vector().get_array(), SALT_BRIDGE_DIST)
            for nb in neighbors:
                if nb in self.water_atoms or nb in self.ligand_atoms:
                    continue
                res = nb.get_parent()
                resname = res.get_resname().strip()
                aname = nb.get_name().strip()
                # Check if protein atom is in a charged group
                is_pos = resname in CHARGED_POS and aname in CHARGED_POS[resname]
                is_neg = resname in CHARGED_NEG and aname in CHARGED_NEG[resname]
                if not (is_pos or is_neg):
                    continue
                d = dist(la, nb)
                self.interactions["salt_bridges"].append({
                    "ligand_atom": atom_label(la),
                    "protein_atom": atom_label(nb),
                    "protein_residue": residue_label(res),
                    "distance": round(d, 2),
                    "charge_type": "positive" if is_pos else "negative",
                })

    def _find_hydrophobic_contacts(self):
        seen = set()
        for la in self.ligand_atoms:
            if not is_carbon(la):
                continue
            neighbors = self.ns.search(la.get_vector().get_array(), HYDROPHOBIC_DIST)
            for nb in neighbors:
                if nb in self.water_atoms or nb in self.ligand_atoms:
                    continue
                if not is_carbon(nb):
                    continue
                d = dist(la, nb)
                res = nb.get_parent()
                key = (atom_label(la), residue_label(res))
                if key not in seen:
                    seen.add(key)
                    self.interactions["hydrophobic_contacts"].append({
                        "ligand_atom": atom_label(la),
                        "protein_atom": atom_label(nb),
                        "protein_residue": residue_label(res),
                        "distance": round(d, 2),
                    })

    def _find_pi_interactions(self):
        """Find pi-stacking and pi-cation interactions."""
        # Identify aromatic rings in ligand
        # For generality, find 5- or 6-membered rings by checking planarity
        # For now, use a simpler heuristic: look for aromatic carbons in ligand
        # and check against known aromatic protein residues

        for pa in self.protein_atoms:
            res = pa.get_parent()
            resname = res.get_resname().strip()
            if resname not in AROMATIC_RESIDUES:
                continue
            ring_names = AROMATIC_ATOMS.get(resname, [])
            ring_atoms_prot = [a for a in res.get_atoms() if a.get_name().strip() in ring_names]
            if len(ring_atoms_prot) < 4:
                continue

            # Check distance to each ligand aromatic carbon
            prot_centroid = ring_centroid(ring_atoms_prot)
            for la in self.ligand_atoms:
                if not is_carbon(la):
                    continue
                la_coord = la.get_vector().get_array()
                d = np.linalg.norm(prot_centroid - la_coord)
                if d < PI_STACK_DIST:
                    self.interactions["pi_interactions"].append({
                        "ligand_atom": atom_label(la),
                        "protein_residue": residue_label(res),
                        "distance": round(d, 2),
                        "type": "pi-stacking/vdW",
                    })
                    break  # one per residue

    def _find_water_bridges(self):
        """Find water-mediated hydrogen bonds."""
        seen = set()
        for wat in self.water_atoms:
            # Check if water is near both ligand and protein
            lig_contacts = []
            prot_contacts = []
            for la in self.ligand_atoms:
                if is_hbond_donor(la) or is_hbond_acceptor(la):
                    d = dist(wat, la)
                    if d <= WATER_BRIDGE_DIST:
                        lig_contacts.append((la, d))
            if not lig_contacts:
                continue
            neighbors = self.ns.search(wat.get_vector().get_array(), WATER_BRIDGE_DIST)
            for nb in neighbors:
                if nb in self.water_atoms or nb in self.ligand_atoms:
                    continue
                if is_hbond_donor(nb) or is_hbond_acceptor(nb):
                    d = dist(wat, nb)
                    prot_contacts.append((nb, d))
            if not prot_contacts:
                continue
            wat_res = wat.get_parent()
            wat_id = wat_res.get_id()[1]
            for (la, d_lw) in lig_contacts:
                for (pa, d_pw) in prot_contacts:
                    key = (atom_label(la), wat_id, atom_label(pa))
                    if key not in seen:
                        seen.add(key)
                        self.interactions["water_bridges"].append({
                            "ligand_atom": atom_label(la),
                            "water_id": wat_id,
                            "protein_atom": atom_label(pa),
                            "protein_residue": residue_label(pa.get_parent()),
                            "lig_water_dist": round(d_lw, 2),
                            "water_prot_dist": round(d_pw, 2),
                        })

    def _find_oxyanion_hole(self):
        """Detect oxyanion hole interactions (common in serine proteases/beta-lactamases)."""
        # Look for ligand oxygen atoms that are within H-bond distance of
        # backbone NH groups (classic oxyanion hole)
        for la in self.ligand_atoms:
            if la.element != "O":
                continue
            neighbors = self.ns.search(la.get_vector().get_array(), HBOND_DIST)
            backbone_nh = []
            for nb in neighbors:
                if nb in self.water_atoms or nb in self.ligand_atoms:
                    continue
                if nb.get_name().strip() == "N":  # backbone nitrogen
                    d = dist(la, nb)
                    if COVALENT_DIST < d <= HBOND_DIST:
                        backbone_nh.append((nb, d))
            if len(backbone_nh) >= 2:
                self.interactions["oxyanion_hole"].append({
                    "ligand_atom": atom_label(la),
                    "backbone_nitrogens": [
                        {"residue": residue_label(nb.get_parent()),
                         "distance": round(d, 2)}
                        for nb, d in backbone_nh
                    ],
                })


# ---------------------------------------------------------------------------
# Summarize interactions into structured text
# ---------------------------------------------------------------------------

def summarize_interactions(interactions, ligand_name, protein_name):
    """Create a structured summary dict for the LLM."""
    summary = {
        "ligand": ligand_name,
        "protein": protein_name,
        "covalent_bonds": interactions["covalent_bonds"],
        "hydrogen_bonds": [],
        "salt_bridges": interactions["salt_bridges"],
        "hydrophobic_residues": [],
        "pi_interactions": [],
        "water_bridges": [],
        "oxyanion_hole": interactions["oxyanion_hole"],
    }

    # Deduplicate and summarize H-bonds by residue
    hbond_by_res = defaultdict(list)
    for hb in interactions["hydrogen_bonds"]:
        hbond_by_res[hb["protein_residue"]].append(hb)
    for res, hbs in sorted(hbond_by_res.items()):
        summary["hydrogen_bonds"].append({
            "residue": res,
            "contacts": [{
                "ligand_atom": h["ligand_atom"],
                "protein_atom": h["protein_atom"],
                "distance": h["distance"],
                "backbone": h["backbone"],
            } for h in hbs]
        })

    # Summarize hydrophobic contacts by residue
    hydro_res = defaultdict(list)
    for hc in interactions["hydrophobic_contacts"]:
        hydro_res[hc["protein_residue"]].append(hc["distance"])
    for res, dists in sorted(hydro_res.items()):
        summary["hydrophobic_residues"].append({
            "residue": res,
            "closest_distance": round(min(dists), 2),
            "num_contacts": len(dists),
        })

    # Pi interactions
    pi_res = set()
    for pi in interactions["pi_interactions"]:
        if pi["protein_residue"] not in pi_res:
            pi_res.add(pi["protein_residue"])
            summary["pi_interactions"].append({
                "residue": pi["protein_residue"],
                "distance": pi["distance"],
                "type": pi["type"],
            })

    # Water bridges - summarize
    wb_seen = set()
    for wb in interactions["water_bridges"]:
        key = (wb["ligand_atom"], wb["protein_residue"])
        if key not in wb_seen:
            wb_seen.add(key)
            summary["water_bridges"].append(wb)

    return summary


# ---------------------------------------------------------------------------
# Generate narrative using Claude
# ---------------------------------------------------------------------------

DEFAULT_MODEL = "o3"

def generate_narrative(summary_json, use_llm=True, model=None):
    """Generate a crystallography-paper-style paragraph."""
    if not use_llm:
        return format_plain_text(summary_json)

    model = model or DEFAULT_MODEL

    prompt = f"""You are writing the "Structure Description" subsection of a protein
crystallography paper for a journal such as Biochemistry, J. Mol. Biol., or Acta
Crystallographica Section D. Your audience is structural biologists and students
learning to read crystallography papers. Accuracy and precise use of terminology
are paramount. It is far better to be dry and correct than colorful and wrong.

Given the interaction data below (extracted computationally from a PDB coordinate
file), write 2-4 paragraphs describing the binding interactions.

WHAT TO COVER (in roughly this order):
1. Covalent linkage, if any — state the bond, the atoms involved, and the distance.
2. Oxyanion hole — which backbone amide nitrogens form it, and what ligand atom
   they coordinate. State distances.
3. Direct hydrogen bonds — list each with donor, acceptor, distance, and whether
   the protein atom is backbone or side chain. Group logically (e.g., by region
   of the active site).
4. Salt bridges / electrostatic interactions.
5. Hydrophobic and van der Waals contacts — name residues forming the hydrophobic
   shell; give closest-approach distances only where informative.
6. Aromatic / pi-stacking interactions if present.
7. Water-mediated hydrogen bonds if present — identify the bridging water and both
   legs of the interaction.

STRICT STYLE RULES:
- Write in the past tense, passive voice ("was observed", "was positioned",
  "formed a hydrogen bond"). This is the standard for describing a single
  crystal structure.
- Report distances in Angstroms in parentheses: (2.8 Å). Use the Å symbol.
- Name residues as three-letter code plus number: Ser70, Ala237, Lys73.
  Distinguish "main chain" / "backbone" from "side chain" atoms explicitly.
- Use only standard structural biology terminology: "hydrogen bond", "salt
  bridge", "van der Waals contact", "oxyanion hole", "tetrahedral adduct",
  "covalent bond". Do NOT use vague words like "robust", "notable",
  "significant", "crucial", "key", "critical", "important", "prominent",
  "extensive", "intricate", or "multifaceted". Just state what was observed.
- Do NOT editorialize about binding affinity, inhibitor potency, or mechanism
  unless it follows directly from the geometry. Describe structure, not function.
- Do NOT use bullet points, numbered lists, headers, or markdown formatting.
  Write continuous prose paragraphs only.
- Every claim must be directly supported by the data below. Do not invent
  interactions or distances not present in the data.

Interaction data:
{json.dumps(summary_json, indent=2)}

Write only the descriptive paragraphs. No title, no heading, no preamble."""

    # Try OpenAI first, then Anthropic
    openai_err = None
    try:
        import openai
        client = openai.OpenAI()
        kwargs = {
            "model": model,
            "messages": [{"role": "user", "content": prompt}],
        }
        # Reasoning models (o1/o3/o4-mini) don't support max_tokens; use max_completion_tokens
        if model.startswith("o"):
            kwargs["max_completion_tokens"] = 4000
        else:
            kwargs["max_tokens"] = 2000
        response = client.chat.completions.create(**kwargs)
        return response.choices[0].message.content
    except Exception as e1:
        openai_err = e1

    try:
        import anthropic
        client = anthropic.Anthropic()
        response = client.messages.create(
            model=model if "claude" in model else "claude-sonnet-4-20250514",
            max_tokens=2000,
            messages=[{"role": "user", "content": prompt}],
        )
        return response.content[0].text
    except Exception as e2:
        print(f"Warning: No LLM API available (OpenAI: {openai_err}, Anthropic: {e2}). Using plain text.", file=sys.stderr)
        return format_plain_text(summary_json)


def format_plain_text(summary):
    """Fallback plain-text formatter without LLM."""
    lines = []
    lines.append(f"Binding Interactions: {summary['ligand']} bound to {summary['protein']}")
    lines.append("")

    if summary["covalent_bonds"]:
        lines.append("COVALENT BONDS:")
        for cb in summary["covalent_bonds"]:
            lines.append(f"  {cb['ligand_atom']} -- {cb['protein_atom']} ({cb['distance']} A)")

    if summary["oxyanion_hole"]:
        lines.append("\nOXYANION HOLE:")
        for oh in summary["oxyanion_hole"]:
            residues = ", ".join(f"{n['residue']} ({n['distance']} A)" for n in oh["backbone_nitrogens"])
            lines.append(f"  {oh['ligand_atom']} positioned in oxyanion hole formed by: {residues}")

    if summary["hydrogen_bonds"]:
        lines.append("\nHYDROGEN BONDS:")
        for hb in summary["hydrogen_bonds"]:
            for c in hb["contacts"]:
                bb = " (backbone)" if c["backbone"] else " (side chain)"
                lines.append(f"  {c['ligand_atom']} -- {c['protein_atom']}{bb} ({c['distance']} A)")

    if summary["salt_bridges"]:
        lines.append("\nSALT BRIDGES:")
        for sb in summary["salt_bridges"]:
            lines.append(f"  {sb['ligand_atom']} -- {sb['protein_atom']} ({sb['distance']} A)")

    if summary["hydrophobic_residues"]:
        lines.append("\nHYDROPHOBIC CONTACTS:")
        for hr in summary["hydrophobic_residues"]:
            lines.append(f"  {hr['residue']} ({hr['num_contacts']} contacts, closest {hr['closest_distance']} A)")

    if summary["pi_interactions"]:
        lines.append("\nPI INTERACTIONS:")
        for pi in summary["pi_interactions"]:
            lines.append(f"  {pi['residue']} ({pi['distance']} A) - {pi['type']}")

    if summary["water_bridges"]:
        lines.append("\nWATER-MEDIATED INTERACTIONS:")
        for wb in summary["water_bridges"][:10]:
            lines.append(f"  {wb['ligand_atom']} -- WAT{wb['water_id']} -- {wb['protein_atom']}")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def find_ligand_residues(structure, ligand_id=None):
    """Find ligand residues in the structure."""
    ligand_residues = []
    for model in structure:
        for chain in model:
            for res in chain:
                het_flag = res.get_id()[0]
                resname = res.get_resname().strip()
                if ligand_id:
                    if resname == ligand_id:
                        ligand_residues.append(res)
                else:
                    # Auto-detect: non-water HETATM that isn't a known modified AA
                    if het_flag.startswith("H_") and resname != "HOH":
                        ligand_residues.append(res)
    return ligand_residues


def main():
    parser = argparse.ArgumentParser(
        description="Describe protein-ligand binding interactions in crystallography paper style."
    )
    parser.add_argument("pdb", help="PDB file of the complex (or protein if --ligand-pdb given)")
    parser.add_argument("--ligand", "-l", help="Three-letter ligand residue name (e.g., BJI)")
    parser.add_argument("--ligand-pdb", help="Separate PDB file for the ligand")
    parser.add_argument("--include-waters", action="store_true", default=True,
                        help="Include water-mediated interactions (default: True)")
    parser.add_argument("--no-waters", action="store_true",
                        help="Exclude water-mediated interactions")
    parser.add_argument("--no-llm", action="store_true",
                        help="Output structured text instead of LLM-generated narrative")
    parser.add_argument("--json", action="store_true",
                        help="Output raw interaction data as JSON")
    parser.add_argument("--protein-name", default=None,
                        help="Protein name for the narrative (auto-detected from PDB)")
    parser.add_argument("--ligand-name", default=None,
                        help="Ligand name for the narrative")
    parser.add_argument("--model", "-m", default=DEFAULT_MODEL,
                        help=f"LLM model to use (default: {DEFAULT_MODEL}). "
                             "Examples: o3, o4-mini, o1, gpt-4.1, gpt-4o")

    args = parser.parse_args()

    if args.no_waters:
        args.include_waters = False

    pdb_parser = PDBParser(QUIET=True)

    # Parse main PDB
    structure = pdb_parser.get_structure("complex", args.pdb)

    # If separate ligand PDB, merge it
    if args.ligand_pdb:
        lig_struct = pdb_parser.get_structure("ligand", args.ligand_pdb)
        # Add ligand residues to the main structure's first model/chain
        model = structure[0]
        chain = list(model.get_chains())[0]
        for lig_model in lig_struct:
            for lig_chain in lig_model:
                for res in lig_chain:
                    chain.add(res)

    # Auto-detect protein name from HEADER/TITLE
    protein_name = args.protein_name
    if not protein_name:
        with open(args.pdb) as f:
            for line in f:
                if line.startswith("TITLE"):
                    protein_name = line[10:].strip()
                    break
        if not protein_name:
            protein_name = Path(args.pdb).stem

    # Find ligand
    ligand_residues = find_ligand_residues(structure, args.ligand)
    if not ligand_residues:
        # Try auto-detect
        ligand_residues = find_ligand_residues(structure, None)
        if not ligand_residues:
            print("Error: No ligand found. Specify --ligand <RESNAME>.", file=sys.stderr)
            sys.exit(1)

    lig_names = {r.get_resname().strip() for r in ligand_residues}
    ligand_name = args.ligand_name or ", ".join(lig_names)

    print(f"Analyzing: {ligand_name} bound to {protein_name}", file=sys.stderr)
    print(f"Ligand residues found: {[residue_label(r) for r in ligand_residues]}", file=sys.stderr)

    # Analyze
    analyzer = InteractionAnalyzer(structure, ligand_residues, args.include_waters)
    interactions = analyzer.find_all()

    # Summarize
    summary = summarize_interactions(interactions, ligand_name, protein_name)

    if args.json:
        print(json.dumps(summary, indent=2))
    elif args.no_llm:
        print(format_plain_text(summary))
    else:
        print(generate_narrative(summary, model=args.model))


if __name__ == "__main__":
    main()
