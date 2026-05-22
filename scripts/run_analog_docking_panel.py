#!/usr/bin/env python
"""Dock the first-round analog library into WT, C326Y, and Y64N production receptors."""

from __future__ import annotations

import csv
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdFMCS, rdMolAlign


ROOT = Path(__file__).resolve().parents[1]

REFERENCE_LIGAND_PDB = ROOT / "data" / "prepared" / "ligands" / "8C03_SZU_bound_reference.pdb"
REFERENCE_LIGAND_H_SDF = ROOT / "data" / "prepared" / "ligands" / "8C03_SZU_bound_reference_H.sdf"
LIBRARY_MANIFEST = ROOT / "results" / "tables" / "analog_library_manifest.csv"
BOX_JSON = ROOT / "config" / "redocking_validation_box_optimized.json"

PROTEIN_DIR = ROOT / "data" / "prepared" / "proteins" / "production"
OUTPUT_DIR = ROOT / "results" / "docking" / "analog_panel"
TABLE_DIR = ROOT / "results" / "tables"

POSE_METRICS_CSV = TABLE_DIR / "analog_cross_state_pose_metrics.csv"
SUMMARY_CSV = TABLE_DIR / "analog_cross_state_summary.csv"
RANKED_HITS_CSV = TABLE_DIR / "final_ranked_hits.csv"
PYMOL_C326Y = ROOT / "scripts" / "view_analog_ranked_hits_c326y.pml"

CONTACT_CUTOFF = 4.0
SUPPORT_RECOVERY_THRESHOLD = 0.875
SUPPORT_RMSD_THRESHOLD = 2.85
SUPPORT_MODE_RANK_THRESHOLD = 10

STATE_TO_RECEPTOR = {
    "WT": PROTEIN_DIR / "8C03_WT_production.pdbqt",
    "C326Y": PROTEIN_DIR / "8C03_C326Y_production.pdbqt",
    "Y64N": PROTEIN_DIR / "8C03_Y64N_production.pdbqt",
}

STATE_TO_PDB = {
    "WT": PROTEIN_DIR / "8C03_WT_production.pdb",
    "C326Y": PROTEIN_DIR / "8C03_C326Y_production.pdb",
    "Y64N": PROTEIN_DIR / "8C03_Y64N_production.pdb",
}

COLORS = {
    "vamifeport": "orange",
    "H1_desfluoro": "tv_green",
    "H2_pyrimidyl": "cyan",
    "H3_regiofluoro": "yellow",
    "T1_benzoxazole": "marine",
    "T2_benzothiazole": "violet",
    "T3_n_methyl_benzimidazole": "salmon",
    "T4_aza_benzimidazole": "limon",
    "L1_tertiary_amine": "hotpink",
    "C1_thiazole_core": "wheat",
}


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


@dataclass(frozen=True)
class LigandEntry:
    analog_id: str
    display_name: str
    role: str
    priority: int
    edited_region: str
    hypothesis: str
    expected_effect: str
    canonical_smiles: str
    pdbqt_path: Path
    sdf_path: Path


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required input missing: {path.relative_to(ROOT)}")


def run_command(args: list[str], stdout_path: Path | None = None) -> None:
    if stdout_path is None:
        subprocess.run(args, check=True, cwd=ROOT)
        return
    with stdout_path.open("w", encoding="utf-8") as handle:
        subprocess.run(args, check=True, cwd=ROOT, stdout=handle)


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


def load_box_and_vina() -> tuple[dict[str, float], dict[str, int | float]]:
    payload = json.loads(BOX_JSON.read_text(encoding="utf-8"))
    return payload["box"], payload["vina"]


def load_bound_reference() -> tuple[Chem.Mol, Chem.Mol]:
    mol_h = Chem.SDMolSupplier(str(REFERENCE_LIGAND_H_SDF), removeHs=False)[0]
    if mol_h is None:
        raise ValueError(f"Could not load bound reference ligand from {REFERENCE_LIGAND_H_SDF.relative_to(ROOT)}")
    return mol_h, Chem.RemoveHs(mol_h)


def load_library_entries() -> list[LigandEntry]:
    entries: list[LigandEntry] = []
    with LIBRARY_MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            entries.append(
                LigandEntry(
                    analog_id=row["analog_id"],
                    display_name=row["display_name"],
                    role=row["role"],
                    priority=int(row["priority"]),
                    edited_region=row["edited_region"],
                    hypothesis=row["hypothesis"],
                    expected_effect=row["expected_effect"],
                    canonical_smiles=row["canonical_smiles"],
                    pdbqt_path=ROOT / row["pdbqt_path"],
                    sdf_path=ROOT / row["sdf_path"],
                )
            )
    return sorted(entries, key=lambda entry: entry.priority)


def residue_contacts_from_pdb_atoms(
    receptor_atoms: list[AtomRecord],
    ligand_atoms: list[AtomRecord],
    cutoff: float,
) -> set[tuple[str, int, str]]:
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


def residue_contacts_from_mol(
    receptor_atoms: list[AtomRecord],
    ligand_mol: Chem.Mol,
    cutoff: float,
) -> set[tuple[str, int, str]]:
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
    if not reference_match or not pose_match:
        raise ValueError("Could not build a heavy-atom mapping for pose RMSD.")
    atom_map = list(zip(pose_match, reference_match))
    return rdMolAlign.AlignMol(pose, reference, atomMap=atom_map)


def parse_vina_scores(out_pdbqt: Path) -> list[float]:
    scores: list[float] = []
    with out_pdbqt.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("REMARK VINA RESULT:"):
                scores.append(float(line.split()[3]))
    return scores


def write_pose_pdb(path: Path, mol: Chem.Mol) -> None:
    Chem.MolToPDBFile(mol, str(path))


def choose_best_recovery_row(pose_rows: list[dict[str, object]]) -> dict[str, object]:
    return max(
        pose_rows,
        key=lambda row: (
            float(row["reference_contact_recovery"]),
            -float(row["rmsd_to_wt_bound_reference"]),
            -float(row["score_kcal_mol"]),
            -int(row["mode_rank"]),
        ),
    )


def dock_ligand_state(
    ligand: LigandEntry,
    state_id: str,
    receptor_pdbqt: Path,
    receptor_pdb: Path,
    box: dict[str, float],
    vina_settings: dict[str, int | float],
    reference_noh: Chem.Mol,
    reference_contacts: set[tuple[str, int, str]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    receptor_atoms = [atom for atom in parse_atoms(receptor_pdb) if atom.line.startswith("ATOM")]
    result_dir = OUTPUT_DIR / ligand.analog_id / state_id
    result_dir.mkdir(parents=True, exist_ok=True)

    out_pdbqt = result_dir / f"{ligand.analog_id}_{state_id}_dock_out.pdbqt"
    out_sdf = result_dir / f"{ligand.analog_id}_{state_id}_dock_out.sdf"
    log_path = result_dir / f"{ligand.analog_id}_{state_id}_dock.log"

    if not out_pdbqt.exists() or not out_sdf.exists() or not log_path.exists():
        run_command(
            [
                "vina",
                "--receptor",
                str(receptor_pdbqt),
                "--ligand",
                str(ligand.pdbqt_path),
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
                str(vina_settings["exhaustiveness"]),
                "--num_modes",
                str(vina_settings["num_modes"]),
                "--energy_range",
                str(vina_settings["energy_range"]),
                "--seed",
                str(vina_settings["seed"]),
                "--out",
                str(out_pdbqt),
                "--verbosity",
                "1",
            ],
            stdout_path=log_path,
        )
        run_command(["obabel", str(out_pdbqt), "-O", str(out_sdf)])

    poses = [mol for mol in Chem.SDMolSupplier(str(out_sdf), removeHs=False) if mol is not None]
    scores = parse_vina_scores(out_pdbqt)
    if len(poses) != len(scores):
        raise ValueError(
            f"Pose/score count mismatch for {ligand.analog_id} {state_id}: {len(poses)} poses vs {len(scores)} scores"
        )

    pose_rows: list[dict[str, object]] = []
    best_rmsd_row: dict[str, object] | None = None
    for rank, (pose, score) in enumerate(zip(poses, scores), start=1):
        rmsd = pose_rmsd_against_reference(reference_noh, pose)
        contacts = residue_contacts_from_mol(receptor_atoms, pose, CONTACT_CUTOFF)
        shared_contacts = contacts & reference_contacts
        novel_contacts = contacts - reference_contacts
        mutation_contacts = {entry for entry in contacts if entry[1] in {64, 326}}
        pose_path = result_dir / f"pose_{rank:02d}.pdb"
        write_pose_pdb(pose_path, pose)

        row: dict[str, object] = {
            "analog_id": ligand.analog_id,
            "display_name": ligand.display_name,
            "role": ligand.role,
            "priority": ligand.priority,
            "edited_region": ligand.edited_region,
            "state_id": state_id,
            "mode_rank": rank,
            "score_kcal_mol": score,
            "rmsd_to_wt_bound_reference": rmsd,
            "contact_count": len(contacts),
            "shared_contact_count": len(shared_contacts),
            "reference_contact_recovery": len(shared_contacts) / len(reference_contacts) if reference_contacts else 0.0,
            "contacts": format_contacts(contacts),
            "shared_contacts": format_contacts(shared_contacts),
            "novel_contacts": format_contacts(novel_contacts),
            "mutation_contacts": format_contacts(mutation_contacts),
            "pose_pdb": str(pose_path.relative_to(ROOT)),
        }
        pose_rows.append(row)
        if best_rmsd_row is None or rmsd < float(best_rmsd_row["rmsd_to_wt_bound_reference"]):
            best_rmsd_row = row

    top_row = pose_rows[0]
    best_recovery_row = choose_best_recovery_row(pose_rows)
    if best_rmsd_row is None:
        raise ValueError(f"No poses parsed for {ligand.analog_id} {state_id}")

    for label, row in (("top", top_row), ("best_recovery", best_recovery_row), ("best_rmsd", best_rmsd_row)):
        src = ROOT / str(row["pose_pdb"])
        dst = result_dir / f"{label}_pose.pdb"
        dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
        row[f"{label}_pose_copy"] = str(dst.relative_to(ROOT))

    summary = {
        "analog_id": ligand.analog_id,
        "display_name": ligand.display_name,
        "role": ligand.role,
        "priority": ligand.priority,
        "edited_region": ligand.edited_region,
        "state_id": state_id,
        "top_mode_rank": int(top_row["mode_rank"]),
        "top_score_kcal_mol": float(top_row["score_kcal_mol"]),
        "top_rmsd_to_wt_bound_reference": float(top_row["rmsd_to_wt_bound_reference"]),
        "top_shared_contact_count": int(top_row["shared_contact_count"]),
        "top_reference_contact_recovery": float(top_row["reference_contact_recovery"]),
        "top_shared_contacts": top_row["shared_contacts"],
        "top_mutation_contacts": top_row["mutation_contacts"],
        "top_pose_pdb": str(top_row["top_pose_copy"]),
        "best_recovery_mode_rank": int(best_recovery_row["mode_rank"]),
        "best_recovery_score_kcal_mol": float(best_recovery_row["score_kcal_mol"]),
        "best_recovery_rmsd_to_wt_bound_reference": float(best_recovery_row["rmsd_to_wt_bound_reference"]),
        "best_recovery_shared_contact_count": int(best_recovery_row["shared_contact_count"]),
        "best_recovery_reference_contact_recovery": float(best_recovery_row["reference_contact_recovery"]),
        "best_recovery_shared_contacts": best_recovery_row["shared_contacts"],
        "best_recovery_mutation_contacts": best_recovery_row["mutation_contacts"],
        "best_recovery_pose_pdb": str(best_recovery_row["best_recovery_pose_copy"]),
        "best_rmsd_mode_rank": int(best_rmsd_row["mode_rank"]),
        "best_rmsd_score_kcal_mol": float(best_rmsd_row["score_kcal_mol"]),
        "best_rmsd_to_wt_bound_reference": float(best_rmsd_row["rmsd_to_wt_bound_reference"]),
        "best_rmsd_shared_contact_count": int(best_rmsd_row["shared_contact_count"]),
        "best_rmsd_reference_contact_recovery": float(best_rmsd_row["reference_contact_recovery"]),
        "best_rmsd_shared_contacts": best_rmsd_row["shared_contacts"],
        "best_rmsd_mutation_contacts": best_rmsd_row["mutation_contacts"],
        "best_rmsd_pose_pdb": str(best_rmsd_row["best_rmsd_pose_copy"]),
        "dock_log": str(log_path.relative_to(ROOT)),
        "dock_sdf": str(out_sdf.relative_to(ROOT)),
        "hypothesis": ligand.hypothesis,
        "expected_effect": ligand.expected_effect,
    }
    return pose_rows, summary


def write_table(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_ranked_hits(summaries: list[dict[str, object]]) -> list[dict[str, object]]:
    summary_map: dict[str, dict[str, dict[str, object]]] = {}
    ligand_meta: dict[str, dict[str, object]] = {}
    for row in summaries:
        analog_id = str(row["analog_id"])
        state_id = str(row["state_id"])
        summary_map.setdefault(analog_id, {})[state_id] = row
        ligand_meta.setdefault(
            analog_id,
            {
                "display_name": row["display_name"],
                "role": row["role"],
                "priority": row["priority"],
                "edited_region": row["edited_region"],
                "hypothesis": row["hypothesis"],
                "expected_effect": row["expected_effect"],
            },
        )

    ranked_rows: list[dict[str, object]] = []
    for analog_id, per_state in summary_map.items():
        wt = per_state["WT"]
        c326y = per_state["C326Y"]
        y64n = per_state["Y64N"]

        def is_supported(summary: dict[str, object]) -> bool:
            return (
                float(summary["best_recovery_reference_contact_recovery"]) >= SUPPORT_RECOVERY_THRESHOLD
                and float(summary["best_recovery_rmsd_to_wt_bound_reference"]) <= SUPPORT_RMSD_THRESHOLD
                and int(summary["best_recovery_mode_rank"]) <= SUPPORT_MODE_RANK_THRESHOLD
            )

        wt_supported = is_supported(wt)
        c326y_supported = is_supported(c326y)
        y64n_supported = is_supported(y64n)
        supported_states = [state for state, ok in (("WT", wt_supported), ("C326Y", c326y_supported), ("Y64N", y64n_supported)) if ok]
        recoveries = [
            float(wt["best_recovery_reference_contact_recovery"]),
            float(c326y["best_recovery_reference_contact_recovery"]),
            float(y64n["best_recovery_reference_contact_recovery"]),
        ]
        rmsds = [
            float(wt["best_recovery_rmsd_to_wt_bound_reference"]),
            float(c326y["best_recovery_rmsd_to_wt_bound_reference"]),
            float(y64n["best_recovery_rmsd_to_wt_bound_reference"]),
        ]
        mode_ranks = [
            int(wt["best_recovery_mode_rank"]),
            int(c326y["best_recovery_mode_rank"]),
            int(y64n["best_recovery_mode_rank"]),
        ]

        if c326y_supported and y64n_supported and wt_supported:
            assessment = "broadly supported lead"
        elif c326y_supported and (y64n_supported or wt_supported):
            assessment = "C326Y-tolerant lead"
        elif c326y_supported:
            assessment = "C326Y-selective lead"
        elif wt_supported and y64n_supported:
            assessment = "C326Y-sensitive analog"
        elif wt_supported or y64n_supported:
            assessment = "single-state hit"
        else:
            assessment = "low-priority scaffold edit"

        rationale = (
            f"{ligand_meta[analog_id]['edited_region']} edit; "
            f"C326Y recovery {float(c326y['best_recovery_reference_contact_recovery']):.3f}, "
            f"Y64N recovery {float(y64n['best_recovery_reference_contact_recovery']):.3f}, "
            f"WT recovery {float(wt['best_recovery_reference_contact_recovery']):.3f}"
        )

        ranked_rows.append(
            {
                "rank": 0,
                "ligand_id": analog_id,
                "series": "vamifeport_round1",
                "states_supported": ",".join(supported_states) if supported_states else "none",
                "overall_assessment": assessment,
                "rationale": rationale,
                "display_name": ligand_meta[analog_id]["display_name"],
                "role": ligand_meta[analog_id]["role"],
                "edited_region": ligand_meta[analog_id]["edited_region"],
                "wt_supported": wt_supported,
                "c326y_supported": c326y_supported,
                "y64n_supported": y64n_supported,
                "wt_recovery": float(wt["best_recovery_reference_contact_recovery"]),
                "c326y_recovery": float(c326y["best_recovery_reference_contact_recovery"]),
                "y64n_recovery": float(y64n["best_recovery_reference_contact_recovery"]),
                "wt_mode_rank": int(wt["best_recovery_mode_rank"]),
                "c326y_mode_rank": int(c326y["best_recovery_mode_rank"]),
                "y64n_mode_rank": int(y64n["best_recovery_mode_rank"]),
                "wt_rmsd": float(wt["best_recovery_rmsd_to_wt_bound_reference"]),
                "c326y_rmsd": float(c326y["best_recovery_rmsd_to_wt_bound_reference"]),
                "y64n_rmsd": float(y64n["best_recovery_rmsd_to_wt_bound_reference"]),
                "mean_recovery": sum(recoveries) / len(recoveries),
                "mean_mode_rank": sum(mode_ranks) / len(mode_ranks),
                "mean_rmsd": sum(rmsds) / len(rmsds),
                "hypothesis": ligand_meta[analog_id]["hypothesis"],
                "expected_effect": ligand_meta[analog_id]["expected_effect"],
                "c326y_best_pose_pdb": c326y["best_recovery_pose_pdb"],
            }
        )

    ranked_rows.sort(
        key=lambda row: (
            -len(str(row["states_supported"]).split(",")) if str(row["states_supported"]) != "none" else 0,
            not bool(row["c326y_supported"]),
            -float(row["c326y_recovery"]),
            int(row["c326y_mode_rank"]),
            -float(row["y64n_recovery"]),
            -float(row["wt_recovery"]),
            -float(row["mean_recovery"]),
            float(row["mean_mode_rank"]),
            float(row["mean_rmsd"]),
            int(row["role"] != "reference"),
            str(row["ligand_id"]),
        )
    )
    for rank, row in enumerate(ranked_rows, start=1):
        row["rank"] = rank
    return ranked_rows


def write_c326y_pymol_view(path: Path, ranked_rows: list[dict[str, object]], top_n: int = 5) -> None:
    top_rows = ranked_rows[:top_n]
    lines = [
        f"load {STATE_TO_PDB['C326Y'].relative_to(ROOT)}, C326Y_receptor",
        f"load {REFERENCE_LIGAND_PDB.relative_to(ROOT)}, bound_vamifeport",
        "hide everything",
        "show cartoon, C326Y_receptor",
        "color tv_green, C326Y_receptor",
        "show sticks, bound_vamifeport",
        "color orange, bound_vamifeport",
        "select pocket, C326Y_receptor and chain A and resi 43+64+65+68+144+185+314+317+318+320+323+324+325+326+466+469+470+473+501+504+508",
        "show sticks, pocket",
        "color cyan, pocket",
        "set cartoon_transparency, 0.35",
    ]
    for row in top_rows:
        analog_id = str(row["ligand_id"])
        pose_path = ROOT / str(row["c326y_best_pose_pdb"])
        object_name = f"{analog_id}_C326Y_best"
        lines.append(f"load {pose_path.relative_to(ROOT)}, {object_name}")
        lines.append(f"show sticks, {object_name}")
        lines.append(f"color {COLORS.get(analog_id, 'white')}, {object_name}")
        lines.append(f"disable {object_name}")
    lines.extend(
        [
            "set stick_radius, 0.18, bound_vamifeport",
            "zoom bound_vamifeport, 12",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    require_file(REFERENCE_LIGAND_PDB)
    require_file(REFERENCE_LIGAND_H_SDF)
    require_file(LIBRARY_MANIFEST)
    require_file(BOX_JSON)
    for path in list(STATE_TO_RECEPTOR.values()) + list(STATE_TO_PDB.values()):
        require_file(path)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    ligands = load_library_entries()
    box, vina_settings = load_box_and_vina()
    reference_ligand_atoms = parse_atoms(REFERENCE_LIGAND_PDB)
    _, reference_noh = load_bound_reference()
    reference_contacts = residue_contacts_from_pdb_atoms(
        [atom for atom in parse_atoms(STATE_TO_PDB["WT"]) if atom.line.startswith("ATOM")],
        reference_ligand_atoms,
        CONTACT_CUTOFF,
    )

    all_pose_rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    for ligand in ligands:
        for state_id in ("WT", "C326Y", "Y64N"):
            pose_rows, summary = dock_ligand_state(
                ligand,
                state_id,
                STATE_TO_RECEPTOR[state_id],
                STATE_TO_PDB[state_id],
                box,
                vina_settings,
                reference_noh,
                reference_contacts,
            )
            all_pose_rows.extend(pose_rows)
            summaries.append(summary)

    write_table(
        POSE_METRICS_CSV,
        all_pose_rows,
        [
            "analog_id",
            "display_name",
            "role",
            "priority",
            "edited_region",
            "state_id",
            "mode_rank",
            "score_kcal_mol",
            "rmsd_to_wt_bound_reference",
            "contact_count",
            "shared_contact_count",
            "reference_contact_recovery",
            "contacts",
            "shared_contacts",
            "novel_contacts",
            "mutation_contacts",
            "pose_pdb",
        ],
    )
    write_table(
        SUMMARY_CSV,
        summaries,
        [
            "analog_id",
            "display_name",
            "role",
            "priority",
            "edited_region",
            "state_id",
            "top_mode_rank",
            "top_score_kcal_mol",
            "top_rmsd_to_wt_bound_reference",
            "top_shared_contact_count",
            "top_reference_contact_recovery",
            "top_shared_contacts",
            "top_mutation_contacts",
            "top_pose_pdb",
            "best_recovery_mode_rank",
            "best_recovery_score_kcal_mol",
            "best_recovery_rmsd_to_wt_bound_reference",
            "best_recovery_shared_contact_count",
            "best_recovery_reference_contact_recovery",
            "best_recovery_shared_contacts",
            "best_recovery_mutation_contacts",
            "best_recovery_pose_pdb",
            "best_rmsd_mode_rank",
            "best_rmsd_score_kcal_mol",
            "best_rmsd_to_wt_bound_reference",
            "best_rmsd_shared_contact_count",
            "best_rmsd_reference_contact_recovery",
            "best_rmsd_shared_contacts",
            "best_rmsd_mutation_contacts",
            "best_rmsd_pose_pdb",
            "dock_log",
            "dock_sdf",
            "hypothesis",
            "expected_effect",
        ],
    )

    ranked_rows = build_ranked_hits(summaries)
    write_table(
        RANKED_HITS_CSV,
        ranked_rows,
        [
            "rank",
            "ligand_id",
            "series",
            "states_supported",
            "overall_assessment",
            "rationale",
            "display_name",
            "role",
            "edited_region",
            "wt_supported",
            "c326y_supported",
            "y64n_supported",
            "wt_recovery",
            "c326y_recovery",
            "y64n_recovery",
            "wt_mode_rank",
            "c326y_mode_rank",
            "y64n_mode_rank",
            "wt_rmsd",
            "c326y_rmsd",
            "y64n_rmsd",
            "mean_recovery",
            "mean_mode_rank",
            "mean_rmsd",
            "hypothesis",
            "expected_effect",
            "c326y_best_pose_pdb",
        ],
    )
    write_c326y_pymol_view(PYMOL_C326Y, ranked_rows)

    print("Docked analog library across WT, C326Y, and Y64N:")
    print(f"  pose metrics: {POSE_METRICS_CSV.relative_to(ROOT)}")
    print(f"  state summary: {SUMMARY_CSV.relative_to(ROOT)}")
    print(f"  ranked hits: {RANKED_HITS_CSV.relative_to(ROOT)}")
    print(f"  PyMOL view: {PYMOL_C326Y.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
