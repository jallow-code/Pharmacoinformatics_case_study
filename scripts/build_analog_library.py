#!/usr/bin/env python
"""Build the first-round vamifeport analog library with a fixed protonated linker state."""

from __future__ import annotations

import csv
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, rdFMCS, rdMolAlign


ROOT = Path(__file__).resolve().parents[1]
REFERENCE_LIGAND_H_SDF = ROOT / "data" / "prepared" / "ligands" / "8C03_SZU_bound_reference_H.sdf"
DESIGN_CSV = ROOT / "data" / "reference" / "analog_designs.csv"
OUTPUT_DIR = ROOT / "data" / "prepared" / "ligands" / "analog_library"
MANIFEST_CSV = ROOT / "results" / "tables" / "analog_library_manifest.csv"
BUNDLE_SDF = OUTPUT_DIR / "analog_library.sdf"

ENV_BIN = Path(sys.executable).resolve().parent
MK_PREPARE_LIGAND = ENV_BIN / "mk_prepare_ligand.py"

EMBED_SEED = 20260521


@dataclass(frozen=True)
class FragmentVariant:
    variant_id: str
    display_name: str
    smiles: str
    design_rationale: str


@dataclass(frozen=True)
class AnalogDesign:
    analog_id: str
    display_name: str
    role: str
    priority: int
    edited_region: str
    head_variant: str
    core_variant: str
    linker_variant: str
    tail_variant: str
    hypothesis: str
    expected_effect: str


HEAD_VARIANTS: dict[str, FragmentVariant] = {
    "parent": FragmentVariant(
        "parent",
        "Fluoropyridyl parent",
        "c1ncccc1F",
        "Validated parent headgroup taken from the bound vamifeport scaffold.",
    ),
    "desfluoro": FragmentVariant(
        "desfluoro",
        "Des-fluoro pyridyl",
        "c1ncccc1",
        "Minimal lipophilicity and electronics edit that removes only the fluorine.",
    ),
    "aza_pyrimidyl": FragmentVariant(
        "aza_pyrimidyl",
        "Fluoropyrimidyl",
        "c1nccnc1F",
        "Adds one aromatic nitrogen as a classic aza-aryl bioisostere.",
    ),
    "fluoro_regio": FragmentVariant(
        "fluoro_regio",
        "Fluoropyridyl regioisomer",
        "c1cncc(F)c1",
        "Repositions the ring nitrogen and fluorine while keeping size and atom count similar.",
    ),
}

CORE_VARIANTS: dict[str, FragmentVariant] = {
    "parent": FragmentVariant(
        "parent",
        "Oxazole-like parent core",
        "c1coc({linker}{tail})n1",
        "Validated parent heteroaryl core from the crystal-bound scaffold.",
    ),
    "thiazole": FragmentVariant(
        "thiazole",
        "Thiazole-like core",
        "c1csc({linker}{tail})n1",
        "Oxygen-to-sulfur core swap to probe central heteroaryl electronics.",
    ),
}

LINKER_VARIANTS: dict[str, FragmentVariant] = {
    "parent": FragmentVariant(
        "parent",
        "Secondary protonated amine",
        "CC[NH+]CC",
        "The validated protonated linker state from the WT re-docking workflow.",
    ),
    "tertiary_amine": FragmentVariant(
        "tertiary_amine",
        "Tertiary ammonium linker",
        "CC[N+](C)CC",
        "Maintains a cationic center while removing the linker NH donor.",
    ),
}

TAIL_VARIANTS: dict[str, FragmentVariant] = {
    "parent": FragmentVariant(
        "parent",
        "Benzimidazole parent tail",
        "c2nc3ccccc3[nH]2",
        "Validated fused heteroaryl tail from the crystal-bound scaffold.",
    ),
    "benzoxazole": FragmentVariant(
        "benzoxazole",
        "Benzoxazole tail",
        "c2nc3ccccc3o2",
        "O-for-N distal bioisostere that removes the donor while preserving fused aromatic size.",
    ),
    "benzothiazole": FragmentVariant(
        "benzothiazole",
        "Benzothiazole tail",
        "c2nc3ccccc3s2",
        "S-for-N distal bioisostere that shifts polarizability and hydrophobicity.",
    ),
    "n_methyl_benzimidazole": FragmentVariant(
        "n_methyl_benzimidazole",
        "N-methyl benzimidazole tail",
        "c2nc3ccccc3n2C",
        "Direct NH-blocking edit that keeps the benzimidazole ring family intact.",
    ),
    "aza_benzimidazole": FragmentVariant(
        "aza_benzimidazole",
        "Aza-benzimidazole tail",
        "c2nc3ccncc3[nH]2",
        "Adds one distal ring nitrogen as a fused aza-aryl bioisostere.",
    ),
}


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required input missing: {path.relative_to(ROOT)}")


def run_command(args: list[str]) -> None:
    subprocess.run(args, check=True, cwd=ROOT)


def load_reference_ligand() -> Chem.Mol:
    mol = Chem.SDMolSupplier(str(REFERENCE_LIGAND_H_SDF), removeHs=False)[0]
    if mol is None:
        raise ValueError(f"Could not load reference ligand from {REFERENCE_LIGAND_H_SDF.relative_to(ROOT)}")
    return mol


def load_designs() -> list[AnalogDesign]:
    designs: list[AnalogDesign] = []
    with DESIGN_CSV.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            designs.append(
                AnalogDesign(
                    analog_id=row["analog_id"],
                    display_name=row["display_name"],
                    role=row["role"],
                    priority=int(row["priority"]),
                    edited_region=row["edited_region"],
                    head_variant=row["head_variant"],
                    core_variant=row["core_variant"],
                    linker_variant=row["linker_variant"],
                    tail_variant=row["tail_variant"],
                    hypothesis=row["hypothesis"],
                    expected_effect=row["expected_effect"],
                )
            )
    return sorted(designs, key=lambda design: design.priority)


def build_smiles(design: AnalogDesign) -> str:
    head = HEAD_VARIANTS[design.head_variant].smiles
    core_template = CORE_VARIANTS[design.core_variant].smiles
    linker = LINKER_VARIANTS[design.linker_variant].smiles
    tail = TAIL_VARIANTS[design.tail_variant].smiles
    return f"O=C(NC{head}){core_template.format(linker=linker, tail=tail)}"


def prepare_conformer(smiles: str, reference_h: Chem.Mol) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse analog SMILES: {smiles}")
    if Chem.GetFormalCharge(mol) != 1:
        raise ValueError(f"Expected a +1 formal charge for {smiles}, found {Chem.GetFormalCharge(mol)}")

    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = EMBED_SEED
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    params.useRandomCoords = True
    embed_status = AllChem.EmbedMolecule(mol_h, params)
    if embed_status != 0:
        raise ValueError(f"Could not embed 3D coordinates for {smiles}")

    AllChem.UFFOptimizeMolecule(mol_h, maxIters=500)
    align_to_reference(mol_h, reference_h)
    return mol_h


def align_to_reference(mol_h: Chem.Mol, reference_h: Chem.Mol) -> float:
    reference_noh = Chem.RemoveHs(Chem.Mol(reference_h))
    pose_noh = Chem.RemoveHs(Chem.Mol(mol_h))
    mcs = rdFMCS.FindMCS(
        [reference_noh, pose_noh],
        bondCompare=rdFMCS.BondCompare.CompareAny,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        timeout=10,
    )
    query = Chem.MolFromSmarts(mcs.smartsString)
    reference_match = reference_noh.GetSubstructMatch(query)
    pose_match = pose_noh.GetSubstructMatch(query)
    if not reference_match or not pose_match:
        raise ValueError("Could not build an MCS alignment between the analog and reference ligand.")
    atom_map = list(zip(pose_match, reference_match))
    return rdMolAlign.AlignMol(mol_h, reference_h, atomMap=atom_map)


def write_manifest(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "analog_id",
        "display_name",
        "role",
        "priority",
        "edited_region",
        "head_variant",
        "core_variant",
        "linker_variant",
        "tail_variant",
        "canonical_smiles",
        "formal_charge",
        "mol_wt",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
        "mcs_rmsd_to_reference",
        "headgroup_rationale",
        "core_rationale",
        "linker_rationale",
        "tail_rationale",
        "hypothesis",
        "expected_effect",
        "sdf_path",
        "pdbqt_path",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    require_file(REFERENCE_LIGAND_H_SDF)
    require_file(DESIGN_CSV)
    require_file(MK_PREPARE_LIGAND)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    reference_h = load_reference_ligand()
    designs = load_designs()

    bundle_writer = Chem.SDWriter(str(BUNDLE_SDF))
    manifest_rows: list[dict[str, object]] = []
    for design in designs:
        smiles = build_smiles(design)
        mol_h = prepare_conformer(smiles, reference_h)
        canonical_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol_h))
        sdf_path = OUTPUT_DIR / f"{design.analog_id}.sdf"
        pdbqt_path = OUTPUT_DIR / f"{design.analog_id}.pdbqt"

        mol_h.SetProp("_Name", design.analog_id)
        mol_h.SetProp("analog_id", design.analog_id)
        mol_h.SetProp("display_name", design.display_name)
        mol_h.SetProp("canonical_smiles", canonical_smiles)
        mol_h.SetProp("edited_region", design.edited_region)
        mol_h.SetProp("role", design.role)
        Chem.MolToMolFile(mol_h, str(sdf_path))
        bundle_writer.write(mol_h)
        run_command([str(MK_PREPARE_LIGAND), "-i", str(sdf_path), "-o", str(pdbqt_path)])

        rmsd = align_to_reference(Chem.Mol(mol_h), reference_h)
        manifest_rows.append(
            {
                "analog_id": design.analog_id,
                "display_name": design.display_name,
                "role": design.role,
                "priority": design.priority,
                "edited_region": design.edited_region,
                "head_variant": design.head_variant,
                "core_variant": design.core_variant,
                "linker_variant": design.linker_variant,
                "tail_variant": design.tail_variant,
                "canonical_smiles": canonical_smiles,
                "formal_charge": Chem.GetFormalCharge(Chem.RemoveHs(mol_h)),
                "mol_wt": round(Descriptors.MolWt(Chem.RemoveHs(mol_h)), 3),
                "tpsa": round(Descriptors.TPSA(Chem.RemoveHs(mol_h)), 2),
                "hbd": Lipinski.NumHDonors(Chem.RemoveHs(mol_h)),
                "hba": Lipinski.NumHAcceptors(Chem.RemoveHs(mol_h)),
                "rotatable_bonds": Lipinski.NumRotatableBonds(Chem.RemoveHs(mol_h)),
                "mcs_rmsd_to_reference": round(rmsd, 3),
                "headgroup_rationale": HEAD_VARIANTS[design.head_variant].design_rationale,
                "core_rationale": CORE_VARIANTS[design.core_variant].design_rationale,
                "linker_rationale": LINKER_VARIANTS[design.linker_variant].design_rationale,
                "tail_rationale": TAIL_VARIANTS[design.tail_variant].design_rationale,
                "hypothesis": design.hypothesis,
                "expected_effect": design.expected_effect,
                "sdf_path": str(sdf_path.relative_to(ROOT)),
                "pdbqt_path": str(pdbqt_path.relative_to(ROOT)),
            }
        )
    bundle_writer.close()

    write_manifest(MANIFEST_CSV, manifest_rows)

    print("Built first-round analog library:")
    print(f"  design table: {DESIGN_CSV.relative_to(ROOT)}")
    print(f"  bundle sdf: {BUNDLE_SDF.relative_to(ROOT)}")
    print(f"  manifest: {MANIFEST_CSV.relative_to(ROOT)}")
    print(f"  ligands: {OUTPUT_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
