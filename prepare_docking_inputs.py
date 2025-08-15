import argparse
import os
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from pathlib import Path
from meeko.utils import run_prepare_receptor4, run_prepare_gpf4


def prepare_receptor(receptor_path: str, output_dir: str):
    receptor_pdbqt = os.path.join(output_dir, "receptor.pdbqt")
    receptor_pdb = Path(receptor_path).resolve()

    run_prepare_receptor4(
        receptor_pdb,
        output_filename=receptor_pdbqt,
        add_msol=True,           # adiciona moléculas de água conservadas
        repairs="checkhydrogens",
        verbose=True,
    )

    return receptor_pdbqt


def prepare_ligand(ligand_path: str, output_dir: str):
    preparator = MoleculePreparation()
    preparator.prepare_from_file(ligand_path)

    ligand_pdbqt = os.path.join(output_dir, "ligand.pdbqt")
    writer = PDBQTWriterLegacy(preparator)
    writer.write(ligand_pdbqt)

    return ligand_pdbqt



def generate_gpf(receptor_pdbqt: str, ligand_pdbqt: str, output_dir: str, padding: float = 5.0):
    gpf_file = os.path.join(output_dir, "receptor.gpf")

    run_prepare_gpf4(
        receptor_filename=receptor_pdbqt,
        ligand_filename=ligand_pdbqt,
        output_filename=gpf_file,
        autobox_add=padding,
        verbose=True,
    )

    return gpf_file


def main():
    parser = argparse.ArgumentParser(description="Prepara arquivos para AutoDock-GPU com Meeko")
    parser.add_argument("--receptor", required=True, help="Caminho para o arquivo receptor .pdb")
    parser.add_argument("--ligante", required=True, help="Caminho para o arquivo do ligante (.pdb, .mol2, etc.)")
    parser.add_argument("-o", "--output", default="autodock_ready", help="Diretório de saída")
    parser.add_argument("--padding", type=float, default=5.0, help="Padding para autobox")

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    receptor_pdbqt = prepare_receptor(args.receptor, args.output)
    ligand_pdbqt = prepare_ligand(args.ligante, args.output)
    gpf_path = generate_gpf(receptor_pdbqt, ligand_pdbqt, args.output, args.padding)

    print("Arquivos gerados:")
    print(" - Receptor preparado:", receptor_pdbqt)
    print(" - Ligante preparado:", ligand_pdbqt)
    print(" - Arquivo GPF:", gpf_path)


if __name__ == "__main__":
    main()
