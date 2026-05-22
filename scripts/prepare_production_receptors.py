#!/usr/bin/env python
"""Build restrained-minimized production receptors for WT and selected SLC40A1 mutants."""

from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from openmm import CustomExternalForce, LangevinIntegrator, Platform, unit
from openmm.app import ForceField, Modeller, NoCutoff, PDBFile, Simulation


ROOT = Path(__file__).resolve().parents[1]
ENV_BIN = Path(sys.executable).resolve().parent
MK_PREPARE_RECEPTOR = ENV_BIN / "mk_prepare_receptor.py"

INPUT_RECEPTORS = {
    "WT": ROOT / "data" / "prepared" / "proteins" / "8C03_WT_receptor.pdb",
    "C326Y": ROOT / "data" / "prepared" / "proteins" / "8C03_C326Y_receptor.pdb",
    "Y64N": ROOT / "data" / "prepared" / "proteins" / "8C03_Y64N_receptor.pdb",
}

POCKET_RESIDUE_CSV = ROOT / "results" / "tables" / "redocking_reference_pocket_residues.csv"
OPTIMIZED_BOX_JSON = ROOT / "config" / "redocking_validation_box_optimized.json"

OUTPUT_DIR = ROOT / "data" / "prepared" / "proteins" / "production"
SUMMARY_CSV = ROOT / "results" / "tables" / "production_receptor_minimization_summary.csv"
PYMOL_VIEW = ROOT / "scripts" / "view_production_receptors.pml"

FORCEFIELD_FILES = ("amber14/protein.ff14SB.xml", "implicit/gbn2.xml")
PH = 7.4
HEAVY_ATOM_RESTRAINT = 1000.0 * unit.kilojoule_per_mole / unit.nanometer**2
MAX_MINIMIZATION_ITERATIONS = 1500
FLEXIBLE_RESIDUE_CHAIN = "A"
ADDITIONAL_FLEXIBLE_RESIDUES = {326}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
TEMPERATURE = 300 * unit.kelvin
FRICTION = 1.0 / unit.picosecond
TIMESTEP = 0.002 * unit.picoseconds
OPENMM_PLATFORM = os.environ.get("OPENMM_PLATFORM", "CPU")


@dataclass(frozen=True)
class AtomRecord:
    atom_name: str
    resname: str
    chain: str
    resid: int
    x: float
    y: float
    z: float


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
                    atom_name=line[12:16].strip(),
                    resname=line[17:20].strip(),
                    chain=line[21].strip(),
                    resid=int(line[22:26]),
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                )
            )
    return atoms


def distance(a: AtomRecord, b: AtomRecord) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def load_flexible_residues() -> set[tuple[str, int]]:
    residues: set[tuple[str, int]] = set()
    with POCKET_RESIDUE_CSV.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            residues.add((row["chain"], int(row["resid"])))
    for resid in ADDITIONAL_FLEXIBLE_RESIDUES:
        residues.add((FLEXIBLE_RESIDUE_CHAIN, resid))
    return residues


def load_box() -> dict[str, float]:
    payload = json.loads(OPTIMIZED_BOX_JSON.read_text(encoding="utf-8"))
    return payload["box"]


def run_command(args: list[str]) -> None:
    subprocess.run(args, check=True, cwd=ROOT)


def strip_hydrogens_from_pdb(input_path: Path, output_path: Path) -> None:
    with input_path.open("r", encoding="utf-8") as source, output_path.open("w", encoding="utf-8") as target:
        for line in source:
            if line.startswith(("ATOM", "HETATM")) and line[76:78].strip() == "H":
                continue
            target.write(line)


def build_restrained_system(modeller: Modeller, forcefield: ForceField, flexible_residues: set[tuple[str, int]]):
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
    )

    restraint = CustomExternalForce("0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    restraint.addGlobalParameter("k", HEAVY_ATOM_RESTRAINT)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    restrained_atoms = 0
    flexible_sidechain_atoms = 0
    backbone_restrained_atoms = 0
    positions = modeller.positions
    for atom, position in zip(modeller.topology.atoms(), positions):
        element = atom.element
        if element is None or element.symbol == "H":
            continue

        residue = atom.residue
        residue_key = (residue.chain.id, int(residue.id))
        is_flexible_residue = residue_key in flexible_residues
        is_backbone = atom.name in BACKBONE_ATOMS

        if is_flexible_residue and not is_backbone:
            flexible_sidechain_atoms += 1
            continue

        restraint.addParticle(atom.index, position.value_in_unit(unit.nanometer))
        restrained_atoms += 1
        if is_flexible_residue and is_backbone:
            backbone_restrained_atoms += 1

    system.addForce(restraint)
    return system, {
        "restrained_heavy_atoms": restrained_atoms,
        "flexible_sidechain_heavy_atoms": flexible_sidechain_atoms,
        "flexible_backbone_heavy_atoms": backbone_restrained_atoms,
    }


def heavy_atom_rmsd(topology, reference_positions, minimized_positions) -> float:
    sq_sum = 0.0
    count = 0
    for atom in topology.atoms():
        element = atom.element
        if element is None or element.symbol == "H":
            continue
        ref = reference_positions[atom.index].value_in_unit(unit.angstrom)
        cur = minimized_positions[atom.index].value_in_unit(unit.angstrom)
        sq_sum += (ref[0] - cur[0]) ** 2 + (ref[1] - cur[1]) ** 2 + (ref[2] - cur[2]) ** 2
        count += 1
    return math.sqrt(sq_sum / count) if count else 0.0


def minimize_receptor(state_id: str, input_path: Path, flexible_residues: set[tuple[str, int]], box: dict[str, float]) -> dict[str, object]:
    print(f"Minimizing {state_id} from {input_path.relative_to(ROOT)}")
    output_prefix = OUTPUT_DIR / f"8C03_{state_id}_production"
    minimized_full_pdb = OUTPUT_DIR / f"8C03_{state_id}_production_fullH.pdb"
    minimized_pdb = output_prefix.with_suffix(".pdb")
    pdbqt_prefix = OUTPUT_DIR / f"8C03_{state_id}_production"
    pdbqt_path = OUTPUT_DIR / f"8C03_{state_id}_production.pdbqt"

    pdb = PDBFile(str(input_path))
    forcefield = ForceField(*FORCEFIELD_FILES)
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield, pH=PH)

    system, restraint_stats = build_restrained_system(modeller, forcefield, flexible_residues)
    integrator = LangevinIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    platform = Platform.getPlatformByName(OPENMM_PLATFORM)
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=MAX_MINIMIZATION_ITERATIONS)
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    minimized_positions = state.getPositions()

    with minimized_full_pdb.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(modeller.topology, minimized_positions, handle, keepIds=True)
    strip_hydrogens_from_pdb(minimized_full_pdb, minimized_pdb)

    run_command(
        [
            str(MK_PREPARE_RECEPTOR),
            "--pdb",
            str(minimized_pdb),
            "-o",
            str(pdbqt_prefix),
            "--box_center",
            str(box["center_x"]),
            str(box["center_y"]),
            str(box["center_z"]),
            "--box_size",
            str(box["size_x"]),
            str(box["size_y"]),
            str(box["size_z"]),
        ]
    )

    minimized_rmsd = heavy_atom_rmsd(modeller.topology, modeller.positions, minimized_positions)
    energy_kj_mol = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    print(f"Finished {state_id}: heavy-atom RMSD {minimized_rmsd:.3f} A")

    return {
        "state_id": state_id,
        "input_pdb": str(input_path.relative_to(ROOT)),
        "output_full_pdb": str(minimized_full_pdb.relative_to(ROOT)),
        "output_pdb": str(minimized_pdb.relative_to(ROOT)),
        "output_pdbqt": str(pdbqt_path.relative_to(ROOT)),
        "heavy_atom_rmsd_angstrom": f"{minimized_rmsd:.3f}",
        "final_potential_energy_kj_mol": f"{energy_kj_mol:.2f}",
        "restrained_heavy_atoms": restraint_stats["restrained_heavy_atoms"],
        "flexible_sidechain_heavy_atoms": restraint_stats["flexible_sidechain_heavy_atoms"],
        "flexible_backbone_heavy_atoms": restraint_stats["flexible_backbone_heavy_atoms"],
    }


def write_summary(rows: list[dict[str, object]]) -> None:
    SUMMARY_CSV.parent.mkdir(parents=True, exist_ok=True)
    with SUMMARY_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "state_id",
                "input_pdb",
                "output_full_pdb",
                "output_pdb",
                "output_pdbqt",
                "heavy_atom_rmsd_angstrom",
                "final_potential_energy_kj_mol",
                "restrained_heavy_atoms",
                "flexible_sidechain_heavy_atoms",
                "flexible_backbone_heavy_atoms",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def write_pymol_view(rows: list[dict[str, object]]) -> None:
    colors = {"WT": "slate", "C326Y": "tv_green", "Y64N": "magenta"}
    box_pdb = OUTPUT_DIR / "8C03_WT_production.gpf.pdb"
    lines = [
        f"load {Path(rows[0]['output_pdb']).as_posix()}, WT_min",
        f"load {Path(rows[1]['output_pdb']).as_posix()}, C326Y_min",
        f"load {Path(rows[2]['output_pdb']).as_posix()}, Y64N_min",
    ]
    if box_pdb.exists():
        lines.append(f"load {box_pdb.relative_to(ROOT)}, docking_box")
    lines.extend(
        [
            "hide everything",
            "show cartoon, WT_min",
            "show cartoon, C326Y_min",
            "show cartoon, Y64N_min",
            "color slate, WT_min",
            "color tv_green, C326Y_min",
            "color magenta, Y64N_min",
            "set cartoon_transparency, 0.35",
            "select pocket, (WT_min or C326Y_min or Y64N_min) and chain A and resi 43+64+65+68+144+185+314+317+318+320+323+324+325+326+466+469+470+473+501+504+508",
            "show sticks, pocket",
            "color cyan, pocket",
            "select mut64, (WT_min or C326Y_min or Y64N_min) and chain A and resi 64",
            "select mut326, (WT_min or C326Y_min or Y64N_min) and chain A and resi 326",
            "show sticks, mut64 or mut326",
            "color yellow, mut64 or mut326",
            "zoom pocket, 12",
        ]
    )
    if box_pdb.exists():
        lines.extend(["show lines, docking_box", "color orange, docking_box", "set line_width, 3, docking_box"])
    PYMOL_VIEW.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    for path in INPUT_RECEPTORS.values():
        require_file(path)
    require_file(POCKET_RESIDUE_CSV)
    require_file(OPTIMIZED_BOX_JSON)
    require_file(MK_PREPARE_RECEPTOR)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    flexible_residues = load_flexible_residues()
    box = load_box()
    print(f"Using OpenMM platform: {OPENMM_PLATFORM}")

    rows: list[dict[str, object]] = []
    for state_id, input_path in INPUT_RECEPTORS.items():
        rows.append(minimize_receptor(state_id, input_path, flexible_residues, box))

    write_summary(rows)
    write_pymol_view(rows)

    print("Prepared production receptors with restrained minimization:")
    print(f"  summary: {SUMMARY_CSV.relative_to(ROOT)}")
    print(f"  PyMOL view: {PYMOL_VIEW.relative_to(ROOT)}")
    for row in rows:
        print(f"  {row['state_id']}: {row['output_pdb']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
