#!/usr/bin/env python3
"""
Describe protein-ligand binding interactions in crystallography paper style.

Usage:
    python describe_binding.py complex.pdb --ligand BJI
    python describe_binding.py protein.pdb --ligand-pdb ligand.pdb
    python describe_binding.py 1erm.pdb --ligand BJI --include-waters
"""

import argparse
import os
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
PROTEIN_PROTEIN_SHELL = 6.0  # radius around ligand for protein-protein interactions
WATER_NEAR_DIST = 4.0        # close waters, likely interacting with ligand
WATER_SHELL_DIST = 6.0       # broader shell for displacement analysis
WATER_MATCH_DIST = 1.5       # max shift to consider a water "conserved" vs "displaced"

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

# Monatomic ions — skipped from auto-detect ligand selection unless --include-ions.
# Also skipped from "ligand_atoms" by default so a Mg2+ next to an ATP doesn't
# get treated as the real ligand.
ION_RESNAMES = {
    "NA", "K", "CL", "MG", "CA", "ZN", "FE", "MN", "CU", "NI", "CO", "CD",
    "BR", "IOD", "HG", "SR", "BA", "CS", "RB", "LI", "AL", "CU1", "FE2",
    "3NI", "NI2", "CO3", "CU2", "ZN2", "MO", "W", "V", "SE", "F", "I",
}

# Common crystallization buffers, cryoprotectants, and additives that land
# in HETATM records but aren't the ligand of interest. Skipped from auto-detect.
CRYSTALLIZATION_RESNAMES = {
    "HOH", "DOD", "EPE", "MPD", "TRS", "TRIS", "MES", "BIS", "BTB", "CHES",
    "HEPES", "PIPES", "BES", "MOPS", "ADA", "BME", "DMS", "DMSO", "GOL",
    "EDO", "PEG", "PG4", "PGE", "P6G", "1PE", "2PE", "PE3", "PE4", "PE5",
    "PE8", "15P", "7PE", "M2M", "PG0", "XPE", "MRD",
    "FMT", "ACT", "ACY", "IMD", "SO4", "PO4", "CIT", "TLA", "MLI", "MLA",
    "FLC", "NH4", "NO3", "BCT", "CO3", "CAC", "OXL", "MAE", "TPO",
    "BOG", "LDA", "OGA", "OLC", "OLA", "OLB", "HEX", "UND", "BAM", "DD9",
    "LMT", "DMU", "MC3", "LMN", "C8E", "F6H", "C10", "DAO",
    "EPH", "PEP", "6JZ", "PTL",
    "DTT", "DTD", "MRZ", "TCE", "CME", "NDG",
    # NOTE: UNK/UNL/UNX are "unknown ligand" resnames. In CASP16 ATX targets
    # these are the actual drug-like ligands with unspecified chemistry, so
    # they are NOT skipped here.
}

SKIP_RESNAMES = ION_RESNAMES | CRYSTALLIZATION_RESNAMES


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


def atom_label(atom, show_chain=False):
    """Human-readable atom label like 'Ser70 OG' or 'Ser70/A OG'."""
    res = atom.get_parent()
    resname = res.get_resname().strip()
    resid = res.get_id()[1]
    aname = atom.get_name().strip()
    chain = res.get_parent().get_id()
    if show_chain and chain.strip():
        return f"{resname}{resid}/{chain} {aname}"
    return f"{resname}{resid} {aname}"


def residue_label(res, show_chain=False):
    resname = res.get_resname().strip()
    resid = res.get_id()[1]
    chain = res.get_parent().get_id()
    if show_chain and chain.strip():
        return f"{resname}{resid}/{chain}"
    return f"{resname} {resid}"


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
        protein_chains = set()

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
                        protein_chains.add(chain.get_id())
                    # modified residues that are part of the chain (like BHD)
                    elif res.get_id()[0].startswith("H_") and rname not in SKIP_RESNAMES:
                        # Check if it's linked to the protein backbone
                        self.protein_atoms.extend(res.get_atoms())
                        protein_chains.add(chain.get_id())

        # Show chain IDs in residue labels when the structure has more than
        # one distinct protein chain (e.g., RuvB hexamers, MRP4 dimers).
        self.show_chain = len(protein_chains) > 1

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
            "protein_protein_interactions": [],
            "occupancy_warnings": {},
            "water_analysis": {},
        }

    def find_all(self):
        self._check_occupancy_and_altloc()
        self._find_covalent_bonds()
        self._find_hydrogen_bonds()
        self._find_salt_bridges()
        self._find_hydrophobic_contacts()
        self._find_pi_interactions()
        if self.include_waters:
            self._find_water_bridges()
            self._analyze_waters()
        self._find_oxyanion_hole()
        self._find_nearby_protein_protein_interactions()
        return self.interactions

    def _find_covalent_bonds(self):
        """Detect covalent bonds between ligand and protein using LINK records and distance."""
        seen = set()
        sc = self.show_chain
        for la in self.ligand_atoms:
            if la.element == "H":
                continue
            neighbors = self.ns.search(la.get_vector().get_array(), COVALENT_DIST)
            for nb in neighbors:
                if nb in self.ligand_atoms or nb in self.water_atoms:
                    continue
                if nb.element == "H":
                    continue
                d = dist(la, nb)
                key = (atom_label(la, sc), atom_label(nb, sc))
                if key not in seen:
                    seen.add(key)
                    self.interactions["covalent_bonds"].append({
                        "ligand_atom": atom_label(la, sc),
                        "protein_atom": atom_label(nb, sc),
                        "distance": round(d, 2),
                    })

    def _find_hydrogen_bonds(self):
        seen = set()
        sc = self.show_chain
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
                key = tuple(sorted([atom_label(la, sc), atom_label(nb, sc)]))
                if key not in seen:
                    seen.add(key)
                    res = nb.get_parent()
                    backbone = nb.get_name().strip() in ("N", "CA", "C", "O")
                    self.interactions["hydrogen_bonds"].append({
                        "ligand_atom": atom_label(la, sc),
                        "protein_atom": atom_label(nb, sc),
                        "protein_residue": residue_label(res, sc),
                        "distance": round(d, 2),
                        "backbone": backbone,
                    })

    def _find_salt_bridges(self):
        """Detect salt bridges between charged groups."""
        sc = self.show_chain
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
                    "ligand_atom": atom_label(la, sc),
                    "protein_atom": atom_label(nb, sc),
                    "protein_residue": residue_label(res, sc),
                    "distance": round(d, 2),
                    "charge_type": "positive" if is_pos else "negative",
                })

    def _find_hydrophobic_contacts(self):
        seen = set()
        sc = self.show_chain
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
                key = (atom_label(la, sc), residue_label(res, sc))
                if key not in seen:
                    seen.add(key)
                    self.interactions["hydrophobic_contacts"].append({
                        "ligand_atom": atom_label(la, sc),
                        "protein_atom": atom_label(nb, sc),
                        "protein_residue": residue_label(res, sc),
                        "distance": round(d, 2),
                    })

    def _find_ligand_aromatic_rings(self, bond_cutoff=1.85, planarity_tol=0.25):
        """Detect planar 5- or 6-membered rings in the ligand.

        Pure-geometry detection: build a connectivity graph from heavy-atom
        distances, enumerate smallest cycles of size 5-6, then keep only
        those whose atoms lie in a plane (smallest singular value of the
        centered coordinates < planarity_tol Å, i.e. out-of-plane thickness
        below ~0.25 Å — typical aromatic rings are ~0.05 Å).
        """
        heavy = [a for a in self.ligand_atoms if a.element != "H"]
        n = len(heavy)
        if n < 5:
            return []
        coords = np.array([a.get_vector().get_array() for a in heavy])
        # Adjacency via distance
        adj = [[] for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if np.linalg.norm(coords[i] - coords[j]) < bond_cutoff:
                    adj[i].append(j)
                    adj[j].append(i)
        # Enumerate 5- and 6-membered cycles, dedup via sorted-tuple key.
        rings = set()
        for start in range(n):
            stack = [(start, [start])]
            while stack:
                node, path = stack.pop()
                if len(path) > 6:
                    continue
                for nb in adj[node]:
                    if nb == start and 5 <= len(path) <= 6:
                        rings.add(tuple(sorted(path)))
                    elif nb > start and nb not in path:
                        stack.append((nb, path + [nb]))
        # Keep only planar rings.
        planar = []
        for ring in rings:
            ring_coords = coords[list(ring)]
            centered = ring_coords - ring_coords.mean(axis=0)
            _, sv, _ = np.linalg.svd(centered)
            if sv[-1] < planarity_tol:
                planar.append([heavy[k] for k in ring])
        return planar

    def _find_pi_interactions(self):
        """Find pi-stacking interactions between ligand rings and protein
        aromatic side chains."""
        sc = self.show_chain
        # Collect protein aromatic rings
        prot_rings = []
        seen_res = set()
        for pa in self.protein_atoms:
            res = pa.get_parent()
            rid = (res.get_parent().get_id(), res.get_id())
            if rid in seen_res:
                continue
            resname = res.get_resname().strip()
            if resname not in AROMATIC_RESIDUES:
                continue
            ring_names = AROMATIC_ATOMS.get(resname, [])
            ring_atoms_prot = [a for a in res.get_atoms() if a.get_name().strip() in ring_names]
            if len(ring_atoms_prot) >= 4:
                prot_rings.append((res, ring_atoms_prot))
                seen_res.add(rid)

        lig_rings = self._find_ligand_aromatic_rings()

        # Ring-ring stacking: both ligand and protein rings detected.
        for lr in lig_rings:
            lc = ring_centroid(lr)
            for res, pr in prot_rings:
                pc = ring_centroid(pr)
                d = float(np.linalg.norm(lc - pc))
                if d < PI_STACK_DIST:
                    plane_angle = ring_angle(lr, pr)
                    # Classify parallel (<30°) vs T-shaped (>60°)
                    if plane_angle < 30:
                        stype = "pi-stacking (parallel)"
                    elif plane_angle > 60:
                        stype = "pi-stacking (T-shaped)"
                    else:
                        stype = "pi-stacking (offset)"
                    self.interactions["pi_interactions"].append({
                        "ligand_atom": atom_label(lr[0], sc),
                        "ligand_ring_size": len(lr),
                        "protein_residue": residue_label(res, sc),
                        "distance": round(d, 2),
                        "plane_angle": round(plane_angle, 1),
                        "type": stype,
                    })

        # Fallback for ligands with no detected aromatic ring (e.g., ions,
        # linear cofactors, or rings missed by the planarity heuristic): use
        # the old carbon-centroid heuristic so we don't lose coverage.
        if not lig_rings:
            for res, pr in prot_rings:
                prot_c = ring_centroid(pr)
                for la in self.ligand_atoms:
                    if not is_carbon(la):
                        continue
                    d = float(np.linalg.norm(prot_c - la.get_vector().get_array()))
                    if d < PI_STACK_DIST:
                        self.interactions["pi_interactions"].append({
                            "ligand_atom": atom_label(la, sc),
                            "protein_residue": residue_label(res, sc),
                            "distance": round(d, 2),
                            "type": "pi-stacking/vdW (heuristic)",
                        })
                        break

    def _find_water_bridges(self):
        """Find water-mediated hydrogen bonds."""
        seen = set()
        sc = self.show_chain
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
                    key = (atom_label(la, sc), wat_id, atom_label(pa, sc))
                    if key not in seen:
                        seen.add(key)
                        self.interactions["water_bridges"].append({
                            "ligand_atom": atom_label(la, sc),
                            "water_id": wat_id,
                            "protein_atom": atom_label(pa, sc),
                            "protein_residue": residue_label(pa.get_parent(), sc),
                            "lig_water_dist": round(d_lw, 2),
                            "water_prot_dist": round(d_pw, 2),
                        })

    def _find_oxyanion_hole(self):
        """Detect oxyanion hole interactions (common in serine proteases/beta-lactamases)."""
        sc = self.show_chain
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
                    "ligand_atom": atom_label(la, sc),
                    "backbone_nitrogens": [
                        {"residue": residue_label(nb.get_parent(), sc),
                         "distance": round(d, 2)}
                        for nb, d in backbone_nh
                    ],
                })

    # ------------------------------------------------------------------
    # Feature 1: Protein-protein interactions near the ligand
    # ------------------------------------------------------------------

    def _find_nearby_protein_protein_interactions(self):
        """Detect H-bonds and salt bridges between protein residues within
        PROTEIN_PROTEIN_SHELL of the ligand.  These may be mechanistically
        important (e.g. Lys73-Glu166 in class A beta-lactamases)."""
        sc = self.show_chain

        # Step 1: identify protein residues near the ligand
        near_residues = set()  # set of residue objects
        near_atoms = set()
        for la in self.ligand_atoms:
            neighbors = self.ns.search(la.get_vector().get_array(), PROTEIN_PROTEIN_SHELL)
            for nb in neighbors:
                if nb in self.water_atoms or nb in self.ligand_atoms:
                    continue
                res = nb.get_parent()
                if is_aa(res, standard=True):
                    near_residues.add(res)

        # Collect atoms from near residues
        near_atom_list = []
        for res in near_residues:
            for a in res.get_atoms():
                near_atoms.add(a)
                near_atom_list.append(a)

        if len(near_atom_list) < 2:
            return

        ns_local = NeighborSearch(near_atom_list)

        # Step 2: find H-bonds and salt bridges between these residues
        seen = set()
        results = []

        for atom_a in near_atom_list:
            res_a = atom_a.get_parent()
            neighbors = ns_local.search(atom_a.get_vector().get_array(), SALT_BRIDGE_DIST)
            for atom_b in neighbors:
                res_b = atom_b.get_parent()
                if res_a is res_b:
                    continue
                # Deduplicate
                key = tuple(sorted([
                    (residue_label(res_a, sc), atom_a.get_name().strip()),
                    (residue_label(res_b, sc), atom_b.get_name().strip()),
                ]))
                if key in seen:
                    continue

                d = dist(atom_a, atom_b)
                interaction_type = None

                # Check salt bridge
                rn_a = res_a.get_resname().strip()
                rn_b = res_b.get_resname().strip()
                an_a = atom_a.get_name().strip()
                an_b = atom_b.get_name().strip()
                charged_a = (rn_a in CHARGED_POS and an_a in CHARGED_POS[rn_a]) or \
                            (rn_a in CHARGED_NEG and an_a in CHARGED_NEG[rn_a])
                charged_b = (rn_b in CHARGED_POS and an_b in CHARGED_POS[rn_b]) or \
                            (rn_b in CHARGED_NEG and an_b in CHARGED_NEG[rn_b])
                if charged_a and charged_b and d <= SALT_BRIDGE_DIST:
                    interaction_type = "salt_bridge"

                # Check H-bond
                if d <= HBOND_DIST and d > COVALENT_DIST:
                    da = is_hbond_donor(atom_a) or is_hbond_acceptor(atom_a)
                    db = is_hbond_donor(atom_b) or is_hbond_acceptor(atom_b)
                    if da and db:
                        if (is_hbond_donor(atom_a) and is_hbond_acceptor(atom_b)) or \
                           (is_hbond_acceptor(atom_a) and is_hbond_donor(atom_b)):
                            if interaction_type is None:
                                interaction_type = "hydrogen_bond"
                            else:
                                interaction_type = "salt_bridge/hydrogen_bond"

                if interaction_type is None:
                    continue

                seen.add(key)

                # Compute min distance of each residue to ligand
                min_d_a = min(dist(la, atom_a) for la in self.ligand_atoms)
                min_d_b = min(dist(la, atom_b) for la in self.ligand_atoms)

                results.append({
                    "type": interaction_type,
                    "residue_a": residue_label(res_a, sc),
                    "atom_a": atom_label(atom_a, sc),
                    "residue_b": residue_label(res_b, sc),
                    "atom_b": atom_label(atom_b, sc),
                    "distance": round(d, 2),
                    "min_dist_to_ligand_a": round(min_d_a, 2),
                    "min_dist_to_ligand_b": round(min_d_b, 2),
                })

        # Sort by sum of distances to ligand (most relevant first)
        results.sort(key=lambda x: x["min_dist_to_ligand_a"] + x["min_dist_to_ligand_b"])
        # Keep top 15 to avoid overwhelming the output
        self.interactions["protein_protein_interactions"] = results[:15]

    # ------------------------------------------------------------------
    # Feature 2: Occupancy and alternate conformation warnings
    # ------------------------------------------------------------------

    def _check_occupancy_and_altloc(self):
        """Flag ligand atoms with occupancy < 1.0 or alternate conformations."""
        occupancies = []
        altloc_atoms = defaultdict(list)
        atoms_below_1 = []

        for atom in self.ligand_atoms:
            occ = atom.get_occupancy()
            occupancies.append(occ)
            altloc = atom.get_altloc()
            if altloc and altloc.strip():
                altloc_atoms[altloc.strip()].append(atom.get_name().strip())
            if occ < 1.0:
                atoms_below_1.append((atom.get_name().strip(), occ))

        if not occupancies:
            return

        min_occ = min(occupancies)
        max_occ = max(occupancies)
        frac_below = sum(1 for o in occupancies if o < 1.0) / len(occupancies)

        has_warnings = min_occ < 1.0 or len(altloc_atoms) > 0

        warning = {}
        if has_warnings:
            warning = {
                "has_warnings": True,
                "min_occupancy": round(min_occ, 2),
                "max_occupancy": round(max_occ, 2),
                "fraction_below_1": round(frac_below, 2),
                "num_atoms_total": len(occupancies),
                "num_atoms_below_1": len(atoms_below_1),
            }
            if altloc_atoms:
                warning["altloc_codes"] = sorted(altloc_atoms.keys())
                warning["altloc_atom_counts"] = {k: len(v) for k, v in sorted(altloc_atoms.items())}
            # Sample of low-occupancy atoms
            if atoms_below_1:
                warning["sample_low_occupancy"] = [
                    f"{name} ({occ:.2f})" for name, occ in atoms_below_1[:8]
                ]
            # Build human-readable warning
            parts = []
            if min_occ < 1.0:
                parts.append(f"Ligand refined at {min_occ:.0%}-{max_occ:.0%} occupancy "
                             f"({len(atoms_below_1)}/{len(occupancies)} atoms below 1.0)")
            if altloc_atoms:
                codes = ", ".join(sorted(altloc_atoms.keys()))
                parts.append(f"Alternate conformations present (altloc codes: {codes})")
            parts.append("Reported distances should be interpreted with caution.")
            warning["reliability_warning"] = ". ".join(parts)
        else:
            warning = {"has_warnings": False}

        self.interactions["occupancy_warnings"] = warning

    # ------------------------------------------------------------------
    # Feature 3: Water analysis
    # ------------------------------------------------------------------

    def _analyze_waters(self):
        """Inventory ordered waters near the ligand."""
        sc = self.show_chain
        nearby_waters = []

        for wat in self.water_atoms:
            wat_coord = wat.get_vector().get_array()
            # Min distance to any ligand atom
            min_d = min(
                (np.linalg.norm(wat_coord - la.get_vector().get_array()) for la in self.ligand_atoms),
                default=999.0,
            )
            if min_d > WATER_SHELL_DIST:
                continue

            wat_res = wat.get_parent()
            wat_id = wat_res.get_id()[1]
            b_factor = round(wat.get_bfactor(), 1)

            # Find protein contacts for this water
            prot_contacts = []
            neighbors = self.ns.search(wat_coord, HBOND_DIST)
            for nb in neighbors:
                if nb in self.water_atoms or nb in self.ligand_atoms:
                    continue
                if is_hbond_donor(nb) or is_hbond_acceptor(nb):
                    d = dist(wat, nb)
                    prot_contacts.append({
                        "atom": atom_label(nb, sc),
                        "residue": residue_label(nb.get_parent(), sc),
                        "distance": round(d, 2),
                    })

            # Find ligand contacts
            lig_contacts = []
            for la in self.ligand_atoms:
                d = np.linalg.norm(wat_coord - la.get_vector().get_array())
                if d <= HBOND_DIST and (is_hbond_donor(la) or is_hbond_acceptor(la)):
                    lig_contacts.append({
                        "atom": atom_label(la, sc),
                        "distance": round(d, 2),
                    })

            nearby_waters.append({
                "water_id": wat_id,
                "min_dist_to_ligand": round(min_d, 2),
                "b_factor": b_factor,
                "within_4A": bool(min_d <= WATER_NEAR_DIST),
                "protein_contacts": prot_contacts[:5],
                "ligand_contacts": lig_contacts,
            })

        nearby_waters.sort(key=lambda w: w["min_dist_to_ligand"])

        self.interactions["water_analysis"] = {
            "total_waters_in_structure": len(self.water_atoms),
            "waters_within_4A": sum(1 for w in nearby_waters if w["within_4A"]),
            "waters_within_6A": len(nearby_waters),
            "nearby_waters": nearby_waters,
        }


def compare_waters(complex_structure, reference_structure, ligand_residues):
    """Compare active-site waters between complex and reference (apo) structure.

    Structures must be pre-superposed (same crystal form / space group).
    Returns a dict with conserved, displaced, and new water lists.
    """
    # Get ligand atom coordinates for defining active site
    lig_coords = []
    for res in ligand_residues:
        for a in res.get_atoms():
            lig_coords.append(a.get_vector().get_array())
    lig_coords = np.array(lig_coords)

    def waters_near_site(structure, shell=WATER_SHELL_DIST):
        """Find waters within shell of ligand position."""
        waters = []
        for model in structure:
            for chain in model:
                for res in chain:
                    if res.get_resname().strip() != "HOH":
                        continue
                    for a in res.get_atoms():
                        coord = a.get_vector().get_array()
                        min_d = np.min(np.linalg.norm(lig_coords - coord, axis=1))
                        if min_d <= shell:
                            waters.append({
                                "id": res.get_id()[1],
                                "coord": coord,
                                "b_factor": round(a.get_bfactor(), 1),
                                "min_dist_to_site": round(min_d, 2),
                            })
        return waters

    complex_waters = waters_near_site(complex_structure)
    ref_waters = waters_near_site(reference_structure)

    # Match waters by spatial proximity
    conserved = []
    displaced = []
    new_waters = []
    matched_complex = set()

    for rw in ref_waters:
        best_d = 999.0
        best_cw = None
        for i, cw in enumerate(complex_waters):
            d = np.linalg.norm(rw["coord"] - cw["coord"])
            if d < best_d:
                best_d = d
                best_cw = i
        if best_d <= WATER_MATCH_DIST and best_cw is not None:
            conserved.append({
                "ref_water_id": rw["id"],
                "complex_water_id": complex_waters[best_cw]["id"],
                "shift": round(best_d, 2),
                "ref_b_factor": rw["b_factor"],
                "complex_b_factor": complex_waters[best_cw]["b_factor"],
            })
            matched_complex.add(best_cw)
        else:
            displaced.append({
                "ref_water_id": rw["id"],
                "ref_b_factor": rw["b_factor"],
                "min_dist_to_site": rw["min_dist_to_site"],
            })

    for i, cw in enumerate(complex_waters):
        if i not in matched_complex:
            new_waters.append({
                "complex_water_id": cw["id"],
                "complex_b_factor": cw["b_factor"],
                "min_dist_to_site": cw["min_dist_to_site"],
            })

    return {
        "reference_active_site_waters": len(ref_waters),
        "complex_active_site_waters": len(complex_waters),
        "conserved": conserved,
        "displaced": displaced,
        "new_in_complex": new_waters,
        "num_conserved": len(conserved),
        "num_displaced": len(displaced),
        "num_new": len(new_waters),
    }


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

    # New features
    if interactions.get("occupancy_warnings", {}).get("has_warnings"):
        summary["occupancy_warnings"] = interactions["occupancy_warnings"]

    if interactions.get("protein_protein_interactions"):
        summary["protein_protein_interactions"] = interactions["protein_protein_interactions"]

    if interactions.get("water_analysis"):
        summary["water_analysis"] = interactions["water_analysis"]

    if interactions.get("water_displacement"):
        summary["water_displacement"] = interactions["water_displacement"]

    return summary


# ---------------------------------------------------------------------------
# Generate narrative using Claude
# ---------------------------------------------------------------------------

DEFAULT_MODEL = "o3"

def generate_narrative(summary_json, use_llm=True, model=None, context=None):
    """Generate a crystallography-paper-style paragraph.

    If `context` is provided, it is inserted before the interaction data block
    and marked as authoritative biology context (catalytic residues, ligand
    chemical class, etc.) that the narrative should respect.
    """
    if not use_llm:
        return format_plain_text(summary_json)

    model = model or DEFAULT_MODEL

    context_block = ""
    if context:
        context_block = f"""

AUTHORITATIVE CONTEXT (treat as ground truth from the structure providers —
do NOT contradict or override with guesses):
{context.strip()}
"""

    prompt = f"""You are writing the "Structure Description" subsection of a protein
crystallography paper for a journal such as Biochemistry, J. Mol. Biol., or Acta
Crystallographica Section D. Your audience is structural biologists and students
learning to read crystallography papers. Accuracy and precise use of terminology
are paramount. It is far better to be dry and correct than colorful and wrong.{context_block}

Given the interaction data below (extracted computationally from a PDB coordinate
file), write 2-4 paragraphs describing the binding interactions.

WHAT TO COVER (in roughly this order):
0. Data quality — if occupancy_warnings are present with has_warnings=true,
   begin with a sentence about ligand occupancy and any alternate conformations.
   State the occupancy value. Warn that distances should be interpreted with
   caution. This is standard practice when reporting partially occupied ligands.
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
8. Ordered water molecules — report the number of ordered waters near the
   ligand (within 4 A and 6 A). If water_displacement data is present,
   describe which conserved waters were retained and which were displaced
   by the ligand, noting implications for binding entropy.
9. Protein-protein interactions near the ligand — if protein_protein_interactions
   data is present, describe hydrogen bonds or salt bridges between protein
   residues within ~6 A of the ligand. These may change upon ligand binding
   and can be mechanistically important (e.g., interactions between catalytic
   residues). State residues, atoms, distances, and proximity to the ligand.

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
- The interaction data lists ONLY protein-ligand contacts. Do NOT describe
  protein-protein hydrogen bonds, backbone-backbone interactions, catalytic
  triad geometry, or any interaction that does not involve the ligand. If the
  data shows zero ligand-protein contacts of a given type, say so plainly; do
  not fabricate catalytic-triad or secondary-structure interactions to fill
  the paragraph. It is acceptable — and in many cases correct — for a
  description to be brief.

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
        # Reasoning models (o1/o3/o4-mini) don't support max_tokens; use
        # max_completion_tokens. This budget includes hidden reasoning tokens,
        # so leave generous room — at 4000 we saw ~20% empty outputs on larger
        # CASP16 targets because the reasoning ate the entire budget before
        # any visible output tokens could be emitted.
        if model.startswith("o"):
            kwargs["max_completion_tokens"] = 16000
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

    if summary.get("occupancy_warnings", {}).get("has_warnings"):
        ow = summary["occupancy_warnings"]
        lines.append("DATA QUALITY WARNINGS:")
        lines.append(f"  {ow['reliability_warning']}")
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

    if summary.get("water_analysis"):
        wa = summary["water_analysis"]
        lines.append(f"\nORDERED WATERS NEAR LIGAND:")
        lines.append(f"  Total waters in structure: {wa['total_waters_in_structure']}")
        lines.append(f"  Waters within 4 A of ligand: {wa['waters_within_4A']}")
        lines.append(f"  Waters within 6 A of ligand: {wa['waters_within_6A']}")
        for w in wa.get("nearby_waters", []):
            tag = "*" if w["within_4A"] else " "
            lines.append(f"  {tag} WAT{w['water_id']} (dist={w['min_dist_to_ligand']} A, "
                         f"B={w['b_factor']})")

    if summary.get("water_displacement"):
        wd = summary["water_displacement"]
        lines.append(f"\nWATER DISPLACEMENT (vs reference):")
        lines.append(f"  Reference active-site waters: {wd['reference_active_site_waters']}")
        lines.append(f"  Complex active-site waters: {wd['complex_active_site_waters']}")
        lines.append(f"  Conserved: {wd['num_conserved']}, Displaced: {wd['num_displaced']}, "
                     f"New: {wd['num_new']}")

    if summary.get("protein_protein_interactions"):
        lines.append("\nPROTEIN-PROTEIN INTERACTIONS NEAR LIGAND:")
        for pp in summary["protein_protein_interactions"]:
            lines.append(f"  {pp['atom_a']} -- {pp['atom_b']} ({pp['distance']} A) "
                         f"[{pp['type']}] "
                         f"(dist to ligand: {pp['min_dist_to_ligand_a']}/{pp['min_dist_to_ligand_b']} A)")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def find_ligand_residues(structure, ligand_id=None, skip_artifacts=True):
    """Find ligand residues in the structure.

    `ligand_id` may be a single 3-letter resname, a comma-separated list
    (e.g. "OAA,OAb"), or None for auto-detect. When given explicitly the
    artifact filter is not applied — the user asked for these by name.

    Auto-detect: any residue that is either a HETATM or a non-standard
    amino acid, skipping water. If skip_artifacts, also drop ions and
    known crystallization additives (see SKIP_RESNAMES).
    """
    wanted_ids = None
    if ligand_id:
        wanted_ids = {s.strip() for s in ligand_id.split(",") if s.strip()}
    ligand_residues = []
    for model in structure:
        for chain in model:
            for res in chain:
                het_flag = res.get_id()[0]
                resname = res.get_resname().strip()
                if wanted_ids is not None:
                    if resname in wanted_ids:
                        ligand_residues.append(res)
                    continue
                # Auto-detect
                if resname == "HOH":
                    continue
                is_het = het_flag.startswith("H_")
                is_nonstandard_aa = (not is_het) and (not is_aa(res, standard=True))
                if not (is_het or is_nonstandard_aa):
                    continue
                if skip_artifacts and resname in SKIP_RESNAMES:
                    continue
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
    parser.add_argument("--include-artifacts", action="store_true",
                        help="Keep ions and crystallization additives in auto-detect "
                             "ligand selection. By default these are skipped — see "
                             "SKIP_RESNAMES at the top of the module.")
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
    parser.add_argument("--context", default=None,
                        help="Path to a text file (or literal string starting with '!') "
                             "with authoritative biology context — catalytic residues, "
                             "ligand chemical class, covalent/noncovalent, etc. Injected "
                             "into the LLM prompt to ground the narrative.")
    parser.add_argument("--reference-pdb", default=None,
                        help="Reference (native/apo) PDB for water displacement analysis. "
                             "Must be pre-superposed onto the complex structure.")

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
    skip_artifacts = not args.include_artifacts
    ligand_residues = find_ligand_residues(structure, args.ligand, skip_artifacts=skip_artifacts)
    if not ligand_residues:
        # Try auto-detect with the same artifact policy
        ligand_residues = find_ligand_residues(structure, None, skip_artifacts=skip_artifacts)
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

    # Water displacement analysis (optional)
    if args.reference_pdb:
        ref_structure = pdb_parser.get_structure("reference", args.reference_pdb)
        displacement = compare_waters(structure, ref_structure, ligand_residues)
        interactions["water_displacement"] = displacement
        print(f"Water displacement: {displacement['num_conserved']} conserved, "
              f"{displacement['num_displaced']} displaced, "
              f"{displacement['num_new']} new", file=sys.stderr)

    # Summarize
    summary = summarize_interactions(interactions, ligand_name, protein_name)

    ctx_text = None
    if args.context:
        if args.context.startswith("!"):
            ctx_text = args.context[1:]
        elif os.path.exists(args.context):
            with open(args.context) as f:
                ctx_text = f.read()
        else:
            ctx_text = args.context

    if args.json:
        print(json.dumps(summary, indent=2))
    elif args.no_llm:
        print(format_plain_text(summary))
    else:
        print(generate_narrative(summary, model=args.model, context=ctx_text))


if __name__ == "__main__":
    main()
