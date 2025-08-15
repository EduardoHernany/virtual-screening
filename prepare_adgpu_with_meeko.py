#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

def which_or_die(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise SystemExit(
            f"Erro: não encontrei '{name}' no PATH.\n"
            "Instale o Meeko (pip install meeko) para ter mk_prepare_receptor.py/mk_prepare_ligand.py."
        )
    return path

def run(cmd: list[str]) -> None:
    print("$", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"Falha ao executar comando (exit {e.returncode}).") from e

def maybe_autofix_receptor(in_pdb: Path, out_pdb: Path, remove_waters: bool, strip_hetero: bool) -> bool:
    """
    Tenta corrigir o PDB com PDBFixer (se instalado):
      - completa átomos ausentes
      - remove heterogêneos (não-água) se strip_hetero=True
      - remove águas se remove_waters=True
    Retorna True se gerou out_pdb, False caso contrário.
    """
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile  # openmm>=7
    except Exception as e:
        print("[INFO] PDBFixer/OpenMM não disponíveis; pulando auto-fix.")
        return False

    print("[INFO] Executando auto-fix com PDBFixer…")
    fixer = PDBFixer(filename=str(in_pdb))

    # Remove heterogêneos conforme flags
    if remove_waters:
        fixer.removeHeterogens(True)   # remove tudo, inclusive água
    elif strip_hetero:
        fixer.removeHeterogens(False)  # remove não-água

    # Completa estrutura
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    # hidrogênios (pH ~ 7.0)
    try:
        fixer.addMissingHydrogens(7.0)
    except Exception:
        pass

    with out_pdb.open("w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    print(f"[OK] PDB “consertado” salvo em: {out_pdb}")
    return True

def main():
    parser = argparse.ArgumentParser(
        description="Prepara receptor (e opcionalmente ligante) para AutoDock-GPU com Meeko, com opções de auto-fix."
    )
    parser.add_argument("--receptor", required=True, type=Path, help="Receptor PDB (ex.: 1cjb_a.pdb)")
    parser.add_argument("--ligand",   required=True, type=Path, help="Ligante (ex.: POP.pdb)")
    parser.add_argument("-o", "--out-prefix", required=True, help="Prefixo de saída (ex.: rec_1iep)")
    parser.add_argument("--padding", type=float, default=5.0, help="Padding em Å (padrão: 5.0)")
    parser.add_argument("--prep-ligand", action="store_true", help="Também preparar ligante em .pdbqt")
    # Tratamento de problemas comuns:
    parser.add_argument("--allow-bad-res", action="store_true", help="Equivalente ao -a/--allow_bad_res (ignora resíduos problemáticos)")
    parser.add_argument("--default-altloc", type=str, default=None, help="AltLoc padrão (tipicamente 'A')")
    parser.add_argument("--wanted-altloc", action="append", default=[], help="AltLoc específico por resíduo (pode repetir). Ex.: A:229:B")
    parser.add_argument("--remove-waters", action="store_true", help="Remover águas do receptor")
    parser.add_argument("--strip-hetero", action="store_true", help="Remover HETATM não-água do receptor (antes do preparo)")
    parser.add_argument("--auto-fix", action="store_true", help="Tentar corrigir receptor com PDBFixer (se instalado)")
    args = parser.parse_args()

    if not args.receptor.exists():
        raise SystemExit(f"Receptor não encontrado: {args.receptor}")
    if not args.ligand.exists():
        raise SystemExit(f"Ligante não encontrado: {args.ligand}")

    mk_prepare_receptor = which_or_die("mk_prepare_receptor.py")
    mk_prepare_ligand  = shutil.which("mk_prepare_ligand.py") if args.prep_ligand else None

    # Pré-processamento opcional do PDB (gera um PDB “limpo” temporário)
    receptor_for_meeko = args.receptor
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        if args.auto-fix or args.strip_hetero or args.remove_waters:
            fixed_pdb = td / "receptor_fixed.pdb"
            if maybe_autofix_receptor(args.receptor, fixed_pdb, args.remove_waters, args.strip_hetero):
                receptor_for_meeko = fixed_pdb

        # Monta comando Meeko (receptor)
        cmd = [
            mk_prepare_receptor,
            "--read_pdb", str(receptor_for_meeko),
            "-o", args.out_prefix,
            "-p", "-g",
            "--box_enveloping", str(args.ligand),
            "--padding", str(args.padding),
        ]
        if args.allow_bad_res:
            cmd.append("--allow_bad_res")
        if args.default_altloc:
            cmd += ["--default_altloc", args.default_altloc]
        for w in args.wanted_altloc:
            cmd += ["--wanted_altloc", w]
        if args.remove_waters and not (args.auto-fix):
            # caso usuário não use auto-fix, delega a remoção de águas ao Meeko (se suportado)
            cmd.append("--remove_waters")

        run(cmd)

        # (Opcional) Prepara ligante
        if mk_prepare_ligand:
            lig_out = Path(args.ligand).with_suffix(".pdbqt")
            run([mk_prepare_ligand, "-i", str(args.ligand), "-o", str(lig_out)])

    # Mensagem final
    out_prefix = Path(args.out_prefix)
    print("\n[OK] Preparação concluída.")
    print(f"- Receptor PDBQT: {out_prefix.with_suffix('.pdbqt').resolve()}")
    print(f"- Arquivo GPF:    {out_prefix.with_suffix('.gpf').resolve()}")
    if args.prep_ligand:
        print(f"- Ligante PDBQT:  {Path(args.ligand).with_suffix('.pdbqt').resolve()}")
    print("\nSe ainda ocorrer falha de template, tente combinar:")
    print("  --allow-bad-res  --default-altloc A  --auto-fix  --strip-hetero  --remove-waters")

if __name__ == "__main__":
    if sys.version_info < (3, 9):
        raise SystemExit("Use Python 3.9+")
    main()
