#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
import xml.etree.ElementTree as ET

def parse_autodock_xml(xml_path: Path):
    """
    Retorna um dicionário com:
      - file: nome do arquivo
      - ligand: nome do ligante (se existir no XML)
      - min_reference_rmsd
      - binding_energy_at_min_rmsd
      - run_id_at_min_rmsd
    Se não encontrar a tabela de RMSD, retorna None.
    """
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except ET.ParseError as e:
        print(f"[WARN] Falha ao parsear {xml_path.name}: {e}")
        return None

    # tenta achar o ligante (campo opcional, mas útil)
    ligand = None
    lig_el = root.find(".//ligand")
    if lig_el is not None and lig_el.text:
        ligand = lig_el.text.strip()

    # percorre a tabela de rmsd
    runs = []
    for run_el in root.findall(".//result/rmsd_table/run"):
        try:
            ref_rmsd = float(run_el.attrib["reference_rmsd"])
            be = float(run_el.attrib["binding_energy"])
            run_id = run_el.attrib.get("run")
            cluster_rmsd = float(run_el.attrib.get("cluster_rmsd", "nan"))
        except (KeyError, ValueError) as e:
            # pulo este run se algo estiver faltando/ruim
            continue
        runs.append(
            {
                "reference_rmsd": ref_rmsd,
                "binding_energy": be,
                "run_id": run_id,
                "cluster_rmsd": cluster_rmsd,
            }
        )

    if not runs:
        return None

    # menor reference_rmsd
    min_rmsd = min(r["reference_rmsd"] for r in runs)

    # filtra todos com o menor RMSD; em caso de empate, escolhe menor binding_energy (mais negativa)
    candidates = [r for r in runs if r["reference_rmsd"] == min_rmsd]
    best = min(candidates, key=lambda r: r["binding_energy"])

    return {
        "file": xml_path.name,
        "ligand": ligand or "",
        "min_reference_rmsd": min_rmsd,
        "binding_energy_at_min_rmsd": best["binding_energy"],
        "run_id_at_min_rmsd": best["run_id"],
    }

def main():
    parser = argparse.ArgumentParser(
        description="Extrai o menor reference_rmsd e a binding_energy correspondente de arquivos AutoDock-GPU XML."
    )
    parser.add_argument("directory", type=Path, help="Caminho do diretório com arquivos .xml")
    parser.add_argument("--csv", type=Path, default=None, help="(Opcional) caminho para salvar um CSV")
    args = parser.parse_args()

    if not args.directory.is_dir():
        raise SystemExit(f"Erro: '{args.directory}' não é um diretório.")

    xml_files = sorted(args.directory.glob("*.xml"))
    if not xml_files:
        raise SystemExit("Nenhum arquivo .xml encontrado no diretório informado.")

    results = []
    for xml_path in xml_files:
        info = parse_autodock_xml(xml_path)
        if info is None:
            print(f"[WARN] Sem dados de rmsd_table em {xml_path.name} (ou XML inválido).")
            continue
        results.append(info)

    if not results:
        raise SystemExit("Nenhum resultado válido foi extraído.")

    # imprime um resumo alinhado
    header = ["file", "ligand", "min_reference_rmsd", "binding_energy_at_min_rmsd", "run_id_at_min_rmsd"]
    print("\t".join(header))
    for r in results:
        print(
            f"{r['file']}\t{r['ligand']}\t{r['min_reference_rmsd']:.2f}\t{r['binding_energy_at_min_rmsd']:.2f}\t{r['run_id_at_min_rmsd']}"
        )

    # salva CSV se pedido
    if args.csv:
        with args.csv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=header)
            writer.writeheader()
            writer.writerows(results)
        print(f"\n[OK] CSV salvo em: {args.csv}")

if __name__ == "__main__":
    main()
