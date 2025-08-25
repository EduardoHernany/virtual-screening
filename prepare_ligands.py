#!/usr/bin/env python3
"""
Preparador de Ligantes para Docking Molecular usando Meeko ou OpenBabel
Converte arquivos SDF em arquivos PDBQT prontos para AutoDock Vina

Uso:
    # Com Meeko (padrão)
    python prepare_ligands.py caminho/para/arquivo.sdf [--output-dir diretorio_saida]

    # Com OpenBabel (mais rápido, menos otimizado)
    python prepare_ligands.py caminho/para/arquivo.sdf --use-obabel [--output-dir diretorio_saida]
"""

import argparse
import os
import sys
import subprocess
import shutil
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy


def prepare_ligand(mol, mol_name, output_dir, pH=7.4):
    """
    Prepara um único ligante para docking.

    Args:
        mol: Molécula RDKit
        mol_name: Nome da molécula
        output_dir: Diretório de saída
        pH: pH para protonação (padrão 7.4)

    Returns:
        list: Lista de arquivos PDBQT gerados
    """
    output_files = []

    # Adiciona hidrogênios se necessário
    mol = Chem.AddHs(mol)

    # Gera conformação 3D se não existir
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol)

    # Prepara a molécula com Meeko
    preparator = MoleculePreparation(
        hydrate=False,  # Não adicionar águas
        rigid_macrocycles=True,  # Macrociclos rígidos
        min_ring_size=5,  # Tamanho mínimo do anel
        max_ring_size=7,  # Tamanho máximo do anel
        merge_these_atom_types=[],  # Tipos de átomos para mesclar
    )

    try:
        # Prepara a molécula
        mol_setups = preparator.prepare(mol)

        # Escreve arquivo PDBQT para cada configuração
        for i, setup in enumerate(mol_setups):
            pdbqt_string, success, error_msg = PDBQTWriterLegacy.write_string(setup)

            if success:
                # Define nome do arquivo de saída
                suffix = f"_{i}" if i > 0 else ""
                output_file = output_dir / f"{mol_name}{suffix}.pdbqt"

                # Escreve o arquivo
                with open(output_file, "w") as f:
                    f.write(pdbqt_string)

                output_files.append(output_file)
                print(f"  ✓ Gerado: {output_file.name}")
            else:
                print(f"  ✗ Erro ao gerar PDBQT: {error_msg}")

    except Exception as e:
        print(f"  ✗ Erro ao preparar molécula: {str(e)}")

    return output_files


def process_sdf_with_obabel(sdf_path, output_dir=None):
    """
    Processa um arquivo SDF usando OpenBabel (mais rápido, menos otimizado).

    Args:
        sdf_path: Caminho para o arquivo SDF
        output_dir: Diretório de saída (opcional)

    Returns:
        dict: Estatísticas do processamento
    """
    sdf_path = Path(sdf_path)

    # Verifica se o arquivo existe
    if not sdf_path.exists():
        print(f"Erro: Arquivo não encontrado: {sdf_path}")
        sys.exit(1)

    # Define diretório de saída
    if output_dir is None:
        output_dir = sdf_path.parent / f"{sdf_path.stem}_pdbqt"
    else:
        output_dir = Path(output_dir)

    # Cria diretório de saída
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n📁 Diretório de saída: {output_dir}")

    # Verifica se OpenBabel está instalado
    obabel_path = shutil.which("obabel")
    if not obabel_path:
        print(f"\n❌ OpenBabel não encontrado!")
        print("Instale com: conda install -c conda-forge openbabel")
        sys.exit(1)

    print(f"\n🔧 Usando OpenBabel: {obabel_path}")
    print(f"📖 Processando arquivo: {sdf_path}")

    # Comando do OpenBabel
    cmd = [
        "obabel",
        "-isdf",
        str(sdf_path),
        "-opdbqt",
        "-O",
        str(output_dir / "mol.pdbqt"),  # Template de saída
        "--split",  # Dividir moléculas em arquivos separados
        "-m",  # Gera múltiplos arquivos numerados
        "-h",  # Adiciona hidrogênios
        "--gen3d",  # Gera coordenadas 3D se não existirem
        "-p",
        "7.4",  # pH para protonação
    ]

    # Estatísticas
    stats = {"total": 0, "sucesso": 0, "falha": 0, "arquivos_gerados": []}

    print("\n🔬 Convertendo com OpenBabel...")
    print(f"   Comando: {' '.join(cmd)}")

    try:
        # Executar OpenBabel
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)

        if result.returncode != 0:
            print(f"\n⚠️  OpenBabel retornou código {result.returncode}")
            if result.stderr:
                print(f"Erro: {result.stderr}")

        # Contar arquivos gerados
        pdbqt_files = sorted(output_dir.glob("mol*.pdbqt"))

        if pdbqt_files:
            print(f"\n✅ {len(pdbqt_files)} arquivos PDBQT gerados")

            # Renomear arquivos para nomes mais descritivos
            print("\n📝 Renomeando arquivos...")

            # Tentar ler o SDF original para obter nomes das moléculas
            try:
                supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
                mol_names = []
                for mol in supplier:
                    if mol is not None and mol.HasProp("_Name"):
                        mol_names.append(mol.GetProp("_Name"))
                    else:
                        mol_names.append(None)
            except:
                mol_names = []

            # Renomear arquivos
            for i, pdbqt_file in enumerate(pdbqt_files):
                # Determinar novo nome
                if i < len(mol_names) and mol_names[i]:
                    new_name = f"{mol_names[i]}.pdbqt"
                else:
                    new_name = f"mol_{i + 1:04d}.pdbqt"

                new_path = output_dir / new_name

                # Evitar sobrescrever arquivos existentes
                if new_path.exists():
                    new_name = f"mol_{i + 1:04d}_{new_name}"
                    new_path = output_dir / new_name

                # Renomear
                pdbqt_file.rename(new_path)
                print(f"   {pdbqt_file.name} → {new_name}")

                stats["arquivos_gerados"].append(new_path)
                stats["sucesso"] += 1

            stats["total"] = len(pdbqt_files)

        else:
            print(f"\n❌ Nenhum arquivo PDBQT foi gerado")
            stats["falha"] = 1

        # Mostrar output do OpenBabel se houver
        if result.stdout:
            print(f"\n📋 Output do OpenBabel:")
            for line in result.stdout.split("\n")[:10]:
                if line.strip():
                    print(f"   {line.strip()}")

    except Exception as e:
        print(f"\n❌ Erro ao executar OpenBabel: {e}")
        stats["falha"] = 1

    return stats


def process_sdf_file(sdf_path, output_dir=None, use_obabel=False):
    """
    Processa um arquivo SDF contendo uma ou mais moléculas.

    Args:
        sdf_path: Caminho para o arquivo SDF
        output_dir: Diretório de saída (opcional)
        use_obabel: Se True, usa OpenBabel ao invés de Meeko

    Returns:
        dict: Estatísticas do processamento
    """
    # Se usar OpenBabel, delegar para função específica
    if use_obabel:
        return process_sdf_with_obabel(sdf_path, output_dir)

    # Caso contrário, usar Meeko (código original)
    sdf_path = Path(sdf_path)

    # Verifica se o arquivo existe
    if not sdf_path.exists():
        print(f"Erro: Arquivo não encontrado: {sdf_path}")
        sys.exit(1)

    # Define diretório de saída
    if output_dir is None:
        output_dir = sdf_path.parent / f"{sdf_path.stem}_pdbqt"
    else:
        output_dir = Path(output_dir)

    # Cria diretório de saída
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n📁 Diretório de saída: {output_dir}")

    # Lê moléculas do arquivo SDF
    print(f"\n📖 Lendo arquivo: {sdf_path}")
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)

    # Estatísticas
    stats = {"total": 0, "sucesso": 0, "falha": 0, "arquivos_gerados": []}

    # Processa cada molécula
    print("\n🔬 Processando moléculas...")
    for idx, mol in enumerate(supplier):
        stats["total"] += 1

        if mol is None:
            print(f"\nMolécula {idx + 1}: ✗ Erro ao ler")
            stats["falha"] += 1
            continue

        # Obtém nome da molécula ou usa índice
        mol_name = (
            mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{idx + 1:04d}"
        )
        print(f"\nMolécula {idx + 1}: {mol_name}")

        # Prepara a molécula
        output_files = prepare_ligand(mol, mol_name, output_dir)

        if output_files:
            stats["sucesso"] += 1
            stats["arquivos_gerados"].extend(output_files)
        else:
            stats["falha"] += 1

    return stats


def main():
    """Função principal."""
    parser = argparse.ArgumentParser(
        description="Prepara ligantes para docking molecular usando Meeko ou OpenBabel",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos:
  # Usar Meeko (padrão - mais preciso, mais lento)
  python prepare_ligands.py molecules.sdf
  python prepare_ligands.py molecules.sdf --output-dir prepared_ligands
  
  # Usar OpenBabel (mais rápido, menos otimizado)
  python prepare_ligands.py molecules.sdf --use-obabel
  python prepare_ligands.py molecules.sdf --use-obabel -o pdbqt_files
  
Comparação:
  Meeko:     Mais preciso, otimizações avançadas, mais lento (~1-5 seg/mol)
  OpenBabel: Mais rápido, conversão simples, menos otimizado (~0.1-0.5 seg/mol)
        """,
    )

    parser.add_argument(
        "sdf_file", help="Caminho para o arquivo SDF contendo os ligantes"
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        help="Diretório de saída para os arquivos PDBQT (padrão: nome_do_arquivo_pdbqt)",
        default=None,
    )

    parser.add_argument(
        "--use-obabel",
        "--obabel",
        action="store_true",
        help="Usar OpenBabel ao invés de Meeko (mais rápido, menos otimizado)",
    )

    args = parser.parse_args()

    # Header
    print("=" * 60)
    print("PREPARADOR DE LIGANTES PARA DOCKING MOLECULAR")
    print("=" * 60)

    # Mostrar método escolhido
    if args.use_obabel:
        print("🔧 Método: OpenBabel (conversão rápida)")
    else:
        print("🔧 Método: Meeko (preparação otimizada)")

    # Processa o arquivo
    stats = process_sdf_file(args.sdf_file, args.output_dir, args.use_obabel)

    # Exibe estatísticas
    print("\n" + "=" * 60)
    print("RESUMO DO PROCESSAMENTO")
    print("=" * 60)

    if args.use_obabel:
        # Estatísticas para OpenBabel
        print(f"📄 Arquivos PDBQT gerados: {len(stats['arquivos_gerados'])}")
        if stats["sucesso"] > 0:
            print(f"\n✅ Conversão concluída com sucesso!")
            print(f"   Método: OpenBabel (rápido)")
        else:
            print(f"\n⚠️  Nenhuma molécula foi convertida.")
    else:
        # Estatísticas para Meeko
        print(f"Total de moléculas: {stats['total']}")
        print(f"✓ Sucesso: {stats['sucesso']}")
        print(f"✗ Falha: {stats['falha']}")
        print(f"📄 Arquivos PDBQT gerados: {len(stats['arquivos_gerados'])}")

        if stats["sucesso"] > 0:
            print(f"\n✅ Processamento concluído com sucesso!")
            print(f"   Método: Meeko (otimizado)")
        else:
            print(f"\n⚠️  Nenhuma molécula foi processada com sucesso.")

    # Sugerir próximo passo
    if stats.get("arquivos_gerados"):
        print(f"\n💡 Próximo passo:")
        print(
            f"   python autodock_gpu_screening.py receptor.pdbqt {args.output_dir or 'ligands_pdbqt'}/"
        )
        sys.exit(0)


if __name__ == "__main__":
    main()

