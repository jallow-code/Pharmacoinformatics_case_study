#!/usr/bin/env python
"""Compare a small set of vamifeport protonation/tautomer states by WT re-docking."""

from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdFMCS, rdMolAlign
from rdkit.Chem.MolStandardize import rdMolStandardize


ROOT = Path(__file__).resolve().parents[1]

REFERENCE_LIGAND_PDB = ROOT / "data" / "prepared" / "ligands" / "8C03_SZU_bound_reference.pdb"
REFERENCE_LIGAND_H_SDF = ROOT / "data" / "prepared" / "ligands" / "8C03_SZU_bound_reference_H.sdf"
REFERENCE_LIGAND_SDF = ROOT / "data" / "prepared" / "ligands" / "8C03_SZU_bound_reference.sdf"
RECEPTOR_PDB = ROOT / "data" / "prepared" / "proteins" / "8C03_WT_receptor.pdb"
RECEPTOR_PDBQT = ROOT / "data" / "prepared" / "proteins" / "8C03_WT_validation.pdbqt"
BOX_JSON = ROOT / "config" / "redocking_validation_box.json"
BOX_PDB = ROOT / "data" / "prepared" / "proteins" / "8C03_WT_validation.gpf.pdb"

STATE_DIR = ROOT / "data" / "prepared" / "ligands" / "redocking_state_panel"
DOCKING_DIR = ROOT / "results" / "docking" / "validation" / "state_panel"
TABLE_DIR = ROOT / "results" / "tables"

STATE_MANIFEST = TABLE_DIR / "redocking_state_panel_inputs.csv"
POSE_METRICS = TABLE_DIR / "redocking_state_panel_pose_metrics.csv"
SUMMARY_TABLE = TABLE_DIR / "redocking_state_panel_summary.csv"
PYMOL_TOP = ROOT / "scripts" / "view_redocking_state_panel_top_poses.pml"
PYMOL_BEST = ROOT / "scripts" / "view_redocking_state_panel_best_rmsd_poses.pml"

CONTACT_CUTOFF = 4.0
VINA_SEED = 20260521
CURRENT_TAUTOMER_SMILES = "O=C(NCc1ncccc1F)c1coc(CCNCCc2nc3ccccc3[nH]2)n1"
ALT_TAUTOMER_SMILES = "O=C(NCc1ncccc1F)C1=COC(=CCNCCc2nc3ccccc3[nH]2)N1"
COLORS = ["magenta", "marine", "tv_green", "salmon"]


@dataclass(frozen=True)
class AtomRecord:
    atom_id: int
    atom_name: str
    resname: str
    chain: str
    resid: int
    x: float
    y: float
    z: float
    line: str


@dataclass
class LigandState:
    state_id: str
    display_name: str
    tautomer_family: str
    protonation_state: str
    formal_charge: int
    canonical_smiles: str
    molecule_h: Chem.Mol
    molecule_noh: Chem.Mol
    sdf_path: Path
    pdbqt_path: Path


def run_command(args: list[str]) -> None:
    subprocess.run(args, check=True, cwd=ROOT)


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required input missing: {path.relative_to(ROOT)}")


def parse_atoms(path: Path) -> list[AtomRecord]:
    atoms: list[AtomRecord] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atoms.append(
                AtomRecord(
                    atom_id=int(line[6:11]),
                    atom_name=line[12:16].strip(),
                    resname=line[17:20].strip(),
                    chain=line[21].strip(),
                    resid=int(line[22:26]),
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    line=line,
                )
            )
    return atoms


def distance_xyz(a_xyz: tuple[float, float, float], b_xyz: tuple[float, float, float]) -> float:
    return math.sqrt((a_xyz[0] - b_xyz[0]) ** 2 + (a_xyz[1] - b_xyz[1]) ** 2 + (a_xyz[2] - b_xyz[2]) ** 2)


def load_bound_reference() -> tuple[Chem.Mol, Chem.Mol]:
    mol_h = Chem.SDMolSupplier(str(REFERENCE_LIGAND_H_SDF), removeHs=False)[0]
    if mol_h is None:
        raise ValueError(f"Could not load reference ligand: {REFERENCE_LIGAND_H_SDF.relative_to(ROOT)}")
    mol_noh = Chem.RemoveHs(mol_h)
    return mol_h, mol_noh


def select_tautomer(base_noh: Chem.Mol, target_smiles: str) -> Chem.Mol:
    enumerator = rdMolStandardize.TautomerEnumerator()
    for tautomer in enumerator.Enumerate(base_noh):
        if Chem.MolToSmiles(tautomer) == target_smiles:
            return Chem.Mol(tautomer)
    raise ValueError(f"Target tautomer not found: {target_smiles}")


def add_hydrogens_with_coords(mol_noh: Chem.Mol) -> Chem.Mol:
    mol_h = Chem.AddHs(Chem.Mol(mol_noh), addCoords=True)
    if mol_h.GetNumConformers() == 0 and mol_noh.GetNumConformers():
        mol_h.AddConformer(mol_noh.GetConformer(), assignId=True)
    return mol_h


def protonate_linker_amine(mol_noh: Chem.Mol) -> Chem.Mol:
    rw = Chem.RWMol(Chem.Mol(mol_noh))
    linker_n_idx: int | None = None
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetIsAromatic() or atom.GetDegree() != 2:
            continue
        neighbors = atom.GetNeighbors()
        if not neighbors or any(neighbor.GetAtomicNum() != 6 for neighbor in neighbors):
            continue
        if all(neighbor.GetHybridization() == Chem.rdchem.HybridizationType.SP3 for neighbor in neighbors):
            linker_n_idx = atom.GetIdx()
            break
    if linker_n_idx is None:
        raise ValueError("Could not identify the linker secondary amine for protonation.")

    atom = rw.GetAtomWithIdx(linker_n_idx)
    atom.SetFormalCharge(1)
    atom.SetNumExplicitHs(1)
    atom.SetNoImplicit(True)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def write_mol_sdf(path: Path, mol: Chem.Mol, props: dict[str, str]) -> None:
    out = Chem.Mol(mol)
    for key, value in props.items():
        out.SetProp(key, value)
    writer = Chem.SDWriter(str(path))
    writer.write(out)
    writer.close()


def build_state_panel() -> list[LigandState]:
    bound_h, bound_noh = load_bound_reference()
    current_noh = Chem.Mol(bound_noh)
    alt_noh = select_tautomer(current_noh, ALT_TAUTOMER_SMILES)

    current_h = Chem.Mol(bound_h)
    current_prot_noh = protonate_linker_amine(current_noh)
    current_prot_h = add_hydrogens_with_coords(current_prot_noh)

    alt_h = add_hydrogens_with_coords(alt_noh)
    alt_prot_noh = protonate_linker_amine(alt_noh)
    alt_prot_h = add_hydrogens_with_coords(alt_prot_noh)

    state_defs = [
        (
            "current_neutral",
            "Current neutral",
            "bound_like",
            "neutral",
            current_h,
            current_noh,
        ),
        (
            "current_linker_protonated",
            "Current linker-protonated",
            "bound_like",
            "linker_amine_protonated",
            current_prot_h,
            current_prot_noh,
        ),
        (
            "alt_tautomer_neutral",
            "Alternative tautomer neutral",
            "alternative_keto_tautomer",
            "neutral",
            alt_h,
            alt_noh,
        ),
        (
            "alt_tautomer_linker_protonated",
            "Alternative tautomer linker-protonated",
            "alternative_keto_tautomer",
            "linker_amine_protonated",
            alt_prot_h,
            alt_prot_noh,
        ),
    ]

    states: list[LigandState] = []
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    with STATE_MANIFEST.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "state_id",
                "display_name",
                "tautomer_family",
                "protonation_state",
                "formal_charge",
                "canonical_smiles",
                "sdf_path",
                "pdbqt_path",
            ]
        )
        for state_id, display_name, tautomer_family, protonation_state, mol_h, mol_noh in state_defs:
            sdf_path = STATE_DIR / f"{state_id}.sdf"
            pdbqt_path = STATE_DIR / f"{state_id}.pdbqt"
            canonical_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol_h))
            formal_charge = Chem.GetFormalCharge(mol_h)

            write_mol_sdf(
                sdf_path,
                mol_h,
                {
                    "state_id": state_id,
                    "display_name": display_name,
                    "tautomer_family": tautomer_family,
                    "protonation_state": protonation_state,
                    "formal_charge": str(formal_charge),
                },
            )
            run_command(["mk_prepare_ligand.py", "-i", str(sdf_path), "-o", str(pdbqt_path)])

            writer.writerow(
                [
                    state_id,
                    display_name,
                    tautomer_family,
                    protonation_state,
                    formal_charge,
                    canonical_smiles,
                    str(sdf_path.relative_to(ROOT)),
                    str(pdbqt_path.relative_to(ROOT)),
                ]
            )

            states.append(
                LigandState(
                    state_id=state_id,
                    display_name=display_name,
                    tautomer_family=tautomer_family,
                    protonation_state=protonation_state,
                    formal_charge=formal_charge,
                    canonical_smiles=canonical_smiles,
                    molecule_h=Chem.Mol(mol_h),
                    molecule_noh=Chem.Mol(mol_noh),
                    sdf_path=sdf_path,
                    pdbqt_path=pdbqt_path,
                )
            )

    return states


def load_box() -> dict[str, float]:
    payload = json.loads(BOX_JSON.read_text(encoding="utf-8"))
    return payload["box"]


def parse_vina_scores(out_pdbqt: Path) -> list[float]:
    scores: list[float] = []
    with out_pdbqt.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("REMARK VINA RESULT:"):
                scores.append(float(line.split()[3]))
    return scores


def write_pose_pdb(path: Path, mol: Chem.Mol) -> None:
    Chem.MolToPDBFile(mol, str(path))


def pose_rmsd_against_reference(reference_noh: Chem.Mol, pose_mol: Chem.Mol) -> float:
    reference = Chem.Mol(reference_noh)
    pose = Chem.RemoveHs(Chem.Mol(pose_mol))
    mcs = rdFMCS.FindMCS(
        [reference, pose],
        bondCompare=rdFMCS.BondCompare.CompareAny,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        timeout=10,
    )
    query = Chem.MolFromSmarts(mcs.smartsString)
    reference_match = reference.GetSubstructMatch(query)
    pose_match = pose.GetSubstructMatch(query)
    if len(reference_match) != reference.GetNumAtoms() or len(pose_match) != pose.GetNumAtoms():
        raise ValueError("Could not find a full heavy-atom mapping for RMSD comparison.")
    atom_map = list(zip(pose_match, reference_match))
    return rdMolAlign.AlignMol(pose, reference, atomMap=atom_map)


def residue_contacts_from_pdb_atoms(receptor_atoms: list[AtomRecord], ligand_atoms: list[AtomRecord], cutoff: float) -> set[tuple[str, int, str]]:
    hits: set[tuple[str, int, str]] = set()
    ligand_xyz = [(atom.x, atom.y, atom.z) for atom in ligand_atoms if atom.line[76:78].strip() != "H"]
    for atom in receptor_atoms:
        if atom.line[76:78].strip() == "H":
            continue
        atom_xyz = (atom.x, atom.y, atom.z)
        for lig_xyz in ligand_xyz:
            if distance_xyz(atom_xyz, lig_xyz) <= cutoff:
                hits.add((atom.chain, atom.resid, atom.resname))
                break
    return hits


def residue_contacts_from_mol(receptor_atoms: list[AtomRecord], ligand_mol: Chem.Mol, cutoff: float) -> set[tuple[str, int, str]]:
    hits: set[tuple[str, int, str]] = set()
    ligand_noh = Chem.RemoveHs(Chem.Mol(ligand_mol))
    conformer = ligand_noh.GetConformer()
    ligand_xyz = [
        (
            conformer.GetAtomPosition(atom_idx).x,
            conformer.GetAtomPosition(atom_idx).y,
            conformer.GetAtomPosition(atom_idx).z,
        )
        for atom_idx in range(ligand_noh.GetNumAtoms())
    ]
    for atom in receptor_atoms:
        if atom.line[76:78].strip() == "H":
            continue
        atom_xyz = (atom.x, atom.y, atom.z)
        for lig_xyz in ligand_xyz:
            if distance_xyz(atom_xyz, lig_xyz) <= cutoff:
                hits.add((atom.chain, atom.resid, atom.resname))
                break
    return hits


def format_contacts(contacts: set[tuple[str, int, str]]) -> str:
    return ";".join(f"{chain}:{resid}:{resname}" for chain, resid, resname in sorted(contacts))


def dock_state(state: LigandState, box: dict[str, float], reference_noh: Chem.Mol, receptor_atoms: list[AtomRecord], reference_contacts: set[tuple[str, int, str]]) -> tuple[list[dict[str, object]], dict[str, object]]:
    state_result_dir = DOCKING_DIR / state.state_id
    state_result_dir.mkdir(parents=True, exist_ok=True)

    out_pdbqt = state_result_dir / f"{state.state_id}_dock_out.pdbqt"
    out_sdf = state_result_dir / f"{state.state_id}_dock_out.sdf"
    log_path = state_result_dir / f"{state.state_id}_dock.log"

    cmd = [
        "vina",
        "--receptor",
        str(RECEPTOR_PDBQT),
        "--ligand",
        str(state.pdbqt_path),
        "--center_x",
        str(box["center_x"]),
        "--center_y",
        str(box["center_y"]),
        "--center_z",
        str(box["center_z"]),
        "--size_x",
        str(box["size_x"]),
        "--size_y",
        str(box["size_y"]),
        "--size_z",
        str(box["size_z"]),
        "--exhaustiveness",
        "16",
        "--num_modes",
        "20",
        "--energy_range",
        "4",
        "--seed",
        str(VINA_SEED),
        "--out",
        str(out_pdbqt),
        "--verbosity",
        "1",
    ]

    with log_path.open("w", encoding="utf-8") as handle:
        subprocess.run(cmd, check=True, cwd=ROOT, stdout=handle)

    run_command(["obabel", str(out_pdbqt), "-O", str(out_sdf)])

    poses = [mol for mol in Chem.SDMolSupplier(str(out_sdf), removeHs=False) if mol is not None]
    scores = parse_vina_scores(out_pdbqt)
    if len(poses) != len(scores):
        raise ValueError(f"Pose/score count mismatch for {state.state_id}: {len(poses)} poses vs {len(scores)} scores")

    pose_rows: list[dict[str, object]] = []
    best_row: dict[str, object] | None = None
    for idx, (pose, score) in enumerate(zip(poses, scores), start=1):
        rmsd = pose_rmsd_against_reference(reference_noh, pose)
        contacts = residue_contacts_from_mol(receptor_atoms, pose, CONTACT_CUTOFF)
        shared_contacts = contacts & reference_contacts
        novel_contacts = contacts - reference_contacts

        pose_path = state_result_dir / f"pose_{idx:02d}.pdb"
        write_pose_pdb(pose_path, pose)

        row: dict[str, object] = {
            "state_id": state.state_id,
            "display_name": state.display_name,
            "mode_rank": idx,
            "score_kcal_mol": score,
            "rmsd_to_reference": rmsd,
            "contact_count": len(contacts),
            "shared_contact_count": len(shared_contacts),
            "reference_contact_recovery": len(shared_contacts) / len(reference_contacts) if reference_contacts else 0.0,
            "contacts": format_contacts(contacts),
            "shared_contacts": format_contacts(shared_contacts),
            "novel_contacts": format_contacts(novel_contacts),
            "pose_pdb": str(pose_path.relative_to(ROOT)),
        }
        pose_rows.append(row)

        if best_row is None or rmsd < float(best_row["rmsd_to_reference"]):
            best_row = row

    top_row = pose_rows[0]
    top_pose_src = ROOT / str(top_row["pose_pdb"])
    best_pose_src = ROOT / str(best_row["pose_pdb"])
    top_pose_dst = state_result_dir / "top_pose.pdb"
    best_pose_dst = state_result_dir / "best_rmsd_pose.pdb"
    shutil.copyfile(top_pose_src, top_pose_dst)
    shutil.copyfile(best_pose_src, best_pose_dst)

    summary = {
        "state_id": state.state_id,
        "display_name": state.display_name,
        "tautomer_family": state.tautomer_family,
        "protonation_state": state.protonation_state,
        "formal_charge": state.formal_charge,
        "canonical_smiles": state.canonical_smiles,
        "top_mode_rank": int(top_row["mode_rank"]),
        "top_score_kcal_mol": float(top_row["score_kcal_mol"]),
        "top_rmsd_to_reference": float(top_row["rmsd_to_reference"]),
        "top_shared_contact_count": int(top_row["shared_contact_count"]),
        "top_reference_contact_recovery": float(top_row["reference_contact_recovery"]),
        "top_shared_contacts": top_row["shared_contacts"],
        "top_pose_pdb": str(top_pose_dst.relative_to(ROOT)),
        "best_rmsd_mode_rank": int(best_row["mode_rank"]),
        "best_rmsd_mode_score_kcal_mol": float(best_row["score_kcal_mol"]),
        "best_rmsd_to_reference": float(best_row["rmsd_to_reference"]),
        "best_rmsd_shared_contact_count": int(best_row["shared_contact_count"]),
        "best_rmsd_reference_contact_recovery": float(best_row["reference_contact_recovery"]),
        "best_rmsd_shared_contacts": best_row["shared_contacts"],
        "best_rmsd_pose_pdb": str(best_pose_dst.relative_to(ROOT)),
        "dock_log": str(log_path.relative_to(ROOT)),
        "dock_sdf": str(out_sdf.relative_to(ROOT)),
    }
    return pose_rows, summary


def write_pose_metrics(rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "state_id",
        "display_name",
        "mode_rank",
        "score_kcal_mol",
        "rmsd_to_reference",
        "contact_count",
        "shared_contact_count",
        "reference_contact_recovery",
        "contacts",
        "shared_contacts",
        "novel_contacts",
        "pose_pdb",
    ]
    with POSE_METRICS.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_summary(rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "state_id",
        "display_name",
        "tautomer_family",
        "protonation_state",
        "formal_charge",
        "canonical_smiles",
        "top_mode_rank",
        "top_score_kcal_mol",
        "top_rmsd_to_reference",
        "top_shared_contact_count",
        "top_reference_contact_recovery",
        "top_shared_contacts",
        "top_pose_pdb",
        "best_rmsd_mode_rank",
        "best_rmsd_mode_score_kcal_mol",
        "best_rmsd_to_reference",
        "best_rmsd_shared_contact_count",
        "best_rmsd_reference_contact_recovery",
        "best_rmsd_shared_contacts",
        "best_rmsd_pose_pdb",
        "dock_log",
        "dock_sdf",
    ]
    with SUMMARY_TABLE.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_pymol_view(path: Path, rows: list[dict[str, object]], pose_key: str, object_suffix: str) -> None:
    residue_selection = "43+64+65+68+144+185+314+317+318+320+323+324+325+466+469+470+473+501+504+508"
    lines = [
        f"load {RECEPTOR_PDB.relative_to(ROOT)}, receptor",
        f"load {REFERENCE_LIGAND_PDB.relative_to(ROOT)}, bound_vamifeport",
        f"load {BOX_PDB.relative_to(ROOT)}, docking_box",
        "",
        "hide everything",
        "show cartoon, receptor",
        "color slate, receptor",
        "show sticks, bound_vamifeport",
        "color orange, bound_vamifeport",
        "show lines, docking_box",
        "color yellow, docking_box",
        "",
        f"select pocket_residues, receptor and chain A and resi {residue_selection}",
        "show sticks, pocket_residues",
        "color cyan, pocket_residues",
        "",
    ]
    for color, row in zip(COLORS, rows):
        object_name = f"{row['state_id']}_{object_suffix}"
        pose_path = ROOT / str(row[pose_key])
        lines.append(f"load {pose_path.relative_to(ROOT)}, {object_name}")
        lines.append(f"show sticks, {object_name}")
        lines.append(f"color {color}, {object_name}")
        lines.append(f"disable {object_name}")
    lines.extend(
        [
            "",
            "enable bound_vamifeport",
            "set stick_radius, 0.18, bound_vamifeport",
            "set line_width, 3, docking_box",
            "set cartoon_transparency, 0.35, receptor",
            "zoom bound_vamifeport, 12",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    for path in [REFERENCE_LIGAND_PDB, REFERENCE_LIGAND_H_SDF, REFERENCE_LIGAND_SDF, RECEPTOR_PDB, RECEPTOR_PDBQT, BOX_JSON]:
        require_file(path)

    DOCKING_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    box = load_box()
    receptor_atoms = [atom for atom in parse_atoms(RECEPTOR_PDB) if atom.line.startswith("ATOM")]
    reference_ligand_atoms = parse_atoms(REFERENCE_LIGAND_PDB)
    reference_contacts = residue_contacts_from_pdb_atoms(receptor_atoms, reference_ligand_atoms, CONTACT_CUTOFF)
    _, reference_noh = load_bound_reference()

    states = build_state_panel()
    all_pose_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for state in states:
        pose_rows, summary = dock_state(state, box, reference_noh, receptor_atoms, reference_contacts)
        all_pose_rows.extend(pose_rows)
        summary_rows.append(summary)

    write_pose_metrics(all_pose_rows)
    write_summary(summary_rows)
    write_pymol_view(PYMOL_TOP, summary_rows, "top_pose_pdb", "top")
    write_pymol_view(PYMOL_BEST, summary_rows, "best_rmsd_pose_pdb", "best_rmsd")

    print("Prepared and docked WT ligand-state panel:")
    print(f"  input states: {STATE_MANIFEST.relative_to(ROOT)}")
    print(f"  pose metrics: {POSE_METRICS.relative_to(ROOT)}")
    print(f"  summary: {SUMMARY_TABLE.relative_to(ROOT)}")
    print(f"  PyMOL top poses: {PYMOL_TOP.relative_to(ROOT)}")
    print(f"  PyMOL best-RMSD poses: {PYMOL_BEST.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
