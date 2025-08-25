#!/usr/bin/env python3
"""
Virtual Screening Multi-Receptor com AutoDock-GPU
Executa docking molecular de múltiplos ligantes contra múltiplos receptores

Uso:
    python autodock_gpu_screening.py --rfiles /path/to/receptors/ --lfiles /path/to/ligands/ --gpus 0,1 --runs 10
"""

import argparse
import os
import sys
import subprocess
import json
import time
import shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import pandas as pd
from itertools import product


class MultiReceptorScreening:
    """Classe para gerenciar virtual screening com múltiplos receptores."""

    def __init__(
        self,
        receptors_dir,
        ligands_dir,
        output_dir="screening_results",
        gpus="0",
        runs=10,
        batch_size=50,
    ):
        """
        Inicializa o screening multi-receptor.

        Args:
            receptors_dir: Diretório com subpastas contendo receptores e mapas
            ligands_dir: Diretório com ligantes PDBQT
            output_dir: Diretório de saída
            gpus: GPUs a usar (ex: "0,1,2")
            runs: Número de runs por ligante
            batch_size: Tamanho do lote
        """
        self.receptors_dir = Path(receptors_dir)
        self.ligands_dir = Path(ligands_dir)
        self.output_dir = Path(output_dir)
        self.gpus = gpus.split(",")
        self.runs = runs
        self.batch_size = batch_size

        # Caminho do AutoDock-GPU
        self.autodock_gpu = Path(os.environ["AUTODOCK_GPU_BIN"]).resolve()

        # Validações
        self._validate_inputs()

        # Coletar receptores e ligantes
        self.receptors = self._collect_receptors()
        self.ligands = self._collect_ligands()

        # Criar estrutura de diretórios
        self._setup_directories()

        # Resultados
        self.results = []

    def _validate_inputs(self):
        """Valida os arquivos e diretórios de entrada."""
        if not self.receptors_dir.exists():
            raise FileNotFoundError(
                f"Diretório de receptores não encontrado: {self.receptors_dir}"
            )

        if not self.ligands_dir.exists():
            raise FileNotFoundError(
                f"Diretório de ligantes não encontrado: {self.ligands_dir}"
            )

        if not self.autodock_gpu.exists():
            raise FileNotFoundError(f"AutoDock-GPU não encontrado: {self.autodock_gpu}")

    def _collect_receptors(self):
        """Coleta todos os receptores com seus mapas."""
        receptors = {}

        print("\n📂 Procurando receptores...")

        # Procurar por arquivos .fld em subdiretórios
        for subdir in self.receptors_dir.iterdir():
            if subdir.is_dir():
                fld_files = list(subdir.glob("*.maps.fld"))
                if fld_files:
                    for fld_file in fld_files:
                        receptor_name = fld_file.stem.replace(".maps", "")
                        receptors[receptor_name] = {
                            "fld": fld_file,
                            "dir": subdir,
                            "pdbqt": subdir / f"{receptor_name}.pdbqt",
                            "maps": list(subdir.glob("*.map")),
                        }
                        print(f"   ✓ {receptor_name}: {fld_file}")

        if not receptors:
            raise ValueError("Nenhum receptor encontrado com arquivos .fld")

        print(f"   Total: {len(receptors)} receptores")
        return receptors

    def _collect_ligands(self):
        """Coleta todos os ligantes PDBQT."""
        ligands = list(self.ligands_dir.glob("*.pdbqt"))

        if not ligands:
            raise ValueError(f"Nenhum ligante .pdbqt encontrado em {self.ligands_dir}")

        print(f"\n📂 Ligantes encontrados: {len(ligands)}")
        return ligands

    def _setup_directories(self):
        """Cria estrutura de diretórios para os resultados."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.docking_dir = self.output_dir / "docking"
        self.logs_dir = self.output_dir / "logs"
        self.results_dir = self.output_dir / "results"

        for dir_path in [self.docking_dir, self.logs_dir, self.results_dir]:
            dir_path.mkdir(exist_ok=True)

        # Criar estrutura de diretórios para resultados: results/receptor/ligante/
        print("\n📁 Criando estrutura de diretórios...")
        for receptor_name in self.receptors:
            receptor_dir = self.results_dir / receptor_name
            receptor_dir.mkdir(exist_ok=True)

            # Criar subdiretório para cada ligante
            for ligand in self.ligands:
                ligand_name = ligand.stem
                ligand_dir = receptor_dir / ligand_name
                ligand_dir.mkdir(exist_ok=True)

        print(f"   ✓ Estrutura criada: results/[receptor]/[ligante]/[arquivos]")

    def create_batch_files(self):
        """
        Cria arquivos de batch com todas as combinações receptor-ligante.
        Formato: 3 linhas por combinação (receptor, ligante, output)

        Returns:
            list: Lista de arquivos batch criados
        """
        batch_files = []
        batch_lines = []
        batch_id = 0
        combinations_count = 0

        print("\n📝 Criando arquivos batch...")

        # Criar todas as combinações receptor-ligante
        for receptor_name in self.receptors:
            receptor_info = self.receptors[receptor_name]

            for ligand in self.ligands:
                ligand_name = ligand.stem

                # Caminho absoluto dos arquivos
                fld_path = str(receptor_info["fld"].absolute())
                ligand_path = str(ligand.absolute())

                # Caminho de saída: results/receptor/ligante/receptor_ligante
                output_base = (
                    self.results_dir
                    / receptor_name
                    / ligand_name
                    / f"{receptor_name}_{ligand_name}"
                )
                output_path = str(output_base.absolute())

                # Adicionar 3 linhas ao batch (receptor, ligante, output)
                batch_lines.append(f"{fld_path}\n")
                batch_lines.append(f"{ligand_path}\n")
                batch_lines.append(f"{output_path}\n")
                combinations_count += 1

                # Se atingiu o tamanho do batch, salvar arquivo
                if combinations_count >= self.batch_size:
                    batch_file = self.docking_dir / f"batch_{batch_id}.txt"
                    with open(batch_file, "w") as f:
                        f.writelines(batch_lines)

                    batch_files.append(batch_file)
                    print(
                        f"   ✓ batch_{batch_id}.txt ({combinations_count} combinações)"
                    )

                    batch_lines = []
                    combinations_count = 0
                    batch_id += 1

        # Salvar último batch se houver linhas restantes
        if batch_lines:
            batch_file = self.docking_dir / f"batch_{batch_id}.txt"
            with open(batch_file, "w") as f:
                f.writelines(batch_lines)

            batch_files.append(batch_file)
            print(f"   ✓ batch_{batch_id}.txt ({combinations_count} combinações)")

        print(f"   Total: {len(batch_files)} arquivos batch")
        return batch_files

    def dock_batch(self, batch_file, batch_id, gpu_id):
        """
        Executa docking de um arquivo batch.

        Args:
            batch_file: Arquivo batch com as combinações
            batch_id: ID do batch
            gpu_id: ID da GPU

        Returns:
            list: Resultados do docking
        """
        batch_results = []

        # Comando do AutoDock-GPU
        cmd = [
            str(self.autodock_gpu),
            "--filelist",
            str(batch_file.absolute()),
            "--nrun",
            str(self.runs),
            "--lsmet",
            "ad",
            # "--devnum", str(gpu_id),
            "--gbest",
            "1",
            "--nev",
            "2500000",
            "--ngen",
            "42000",
            "--lsit",
            "300",
        ]

        # Log file
        log_file = self.logs_dir / f"batch_{batch_id}_gpu_{gpu_id}.log"

        try:
            # Contar combinações no batch (3 linhas por combinação)
            with open(batch_file, "r") as f:
                lines = f.readlines()
                num_combinations = len(lines) // 3

            print(
                f"  🚀 Batch {batch_id} na GPU {gpu_id} ({num_combinations} combinações)"
            )

            start_time = time.time()

            # Executar docking
            result = subprocess.run(
                cmd, cwd=self.docking_dir, capture_output=True, text=True, check=False
            )

            # Salvar log
            with open(log_file, "w") as log:
                log.write(f"Batch file: {batch_file}\n")
                log.write(f"Working directory: {self.docking_dir}\n")
                log.write(f"Command: {' '.join(cmd)}\n")
                log.write(f"Return code: {result.returncode}\n")
                log.write(f"\n=== STDOUT ===\n{result.stdout}\n")
                log.write(f"\n=== STDERR ===\n{result.stderr}\n")

            if result.returncode != 0:
                print(f"  ❌ Erro no batch {batch_id} (código: {result.returncode})")
                print(f"     📄 Log: {log_file}")

                # Mostrar erro
                error_output = result.stderr if result.stderr else result.stdout
                if error_output:
                    error_lines = error_output.split("\n")[:5]
                    for line in error_lines:
                        if line.strip():
                            print(f"     ⚠️  {line.strip()[:100]}")
            else:
                elapsed = time.time() - start_time
                print(f"  ✅ Batch {batch_id} concluído em {elapsed:.1f}s")

                # Verificar se os arquivos foram gerados
                print(f"     🔍 Verificando resultados...")

                # Processar resultados de cada combinação no batch
                batch_results = self._parse_batch_results(batch_file, batch_id)

                if not batch_results:
                    print(
                        f"     ⚠️  Nenhum resultado válido encontrado no batch {batch_id}"
                    )
                    # Listar alguns arquivos para debug
                    print(f"     📁 Verificando diretório de resultados...")
                    for receptor_name in list(self.receptors.keys())[
                        :2
                    ]:  # Verificar apenas 2 receptores
                        receptor_dir = self.results_dir / receptor_name
                        if receptor_dir.exists():
                            dlg_files = list(receptor_dir.glob("*/*.dlg"))[:3]
                            if dlg_files:
                                print(
                                    f"        Encontrados em {receptor_name}: {len(dlg_files)} arquivos DLG"
                                )

        except Exception as e:
            print(f"  ❌ Erro no batch {batch_id}: {e}")
            import traceback

            traceback.print_exc()

        return batch_results

    def _parse_batch_results(self, batch_file, batch_id):
        """
        Extrai resultados de um batch.

        Args:
            batch_file: Arquivo batch processado
            batch_id: ID do batch

        Returns:
            list: Resultados extraídos
        """
        results = []

        # Ler combinações do batch (3 linhas por combinação)
        with open(batch_file, "r") as f:
            lines = f.readlines()

        # Processar cada combinação (grupos de 3 linhas)
        for i in range(0, len(lines), 3):
            if i + 2 < len(lines):
                fld_path = lines[i].strip()
                ligand_path = lines[i + 1].strip()
                output_base = lines[i + 2].strip()

                # Extrair nomes
                receptor_name = Path(fld_path).stem.replace(".maps", "")
                ligand_name = Path(ligand_path).stem

                # Verificar arquivo DLG
                dlg_file = Path(f"{output_base}.dlg")

                if dlg_file.exists():
                    energy = self._extract_energy(dlg_file)
                    if energy is not None:
                        # Salvar caminho relativo se possível, senão apenas o nome
                        try:
                            dlg_relative = dlg_file.relative_to(self.output_dir)
                        except ValueError:
                            # Se não for possível fazer relativo, usar apenas o nome do arquivo
                            dlg_relative = (
                                f"results/{receptor_name}/{ligand_name}/{dlg_file.name}"
                            )

                        results.append(
                            {
                                "receptor": receptor_name,
                                "ligand": ligand_name,
                                "best_energy": energy,
                                "dlg_file": str(dlg_relative),
                                "batch_id": batch_id,
                            }
                        )
                    else:
                        print(f"     ⚠️  Sem energia válida: {dlg_file.name}")
                else:
                    # Tentar caminho alternativo se o arquivo não existir no caminho absoluto
                    alt_dlg_file = (
                        self.results_dir
                        / receptor_name
                        / ligand_name
                        / f"{receptor_name}_{ligand_name}.dlg"
                    )
                    if alt_dlg_file.exists():
                        energy = self._extract_energy(alt_dlg_file)
                        if energy is not None:
                            results.append(
                                {
                                    "receptor": receptor_name,
                                    "ligand": ligand_name,
                                    "best_energy": energy,
                                    "dlg_file": f"results/{receptor_name}/{ligand_name}/{receptor_name}_{ligand_name}.dlg",
                                    "batch_id": batch_id,
                                }
                            )
                    else:
                        print(
                            f"     ⚠️  Arquivo não encontrado: {receptor_name}_{ligand_name}.dlg"
                        )

        print(f"     📊 {len(results)} resultados extraídos do batch {batch_id}")
        return results

    def _extract_energy(self, dlg_file):
        """
        Extrai a melhor energia de um arquivo DLG.

        Args:
            dlg_file: Caminho para o arquivo DLG

        Returns:
            float: Melhor energia ou None
        """
        best_energy = None

        try:
            with open(dlg_file, "r") as f:
                for line in f:
                    if "Estimated Free Energy of Binding" in line:
                        try:
                            # Formato: = -X.XX kcal/mol
                            parts = line.split("=")
                            if len(parts) > 1:
                                energy_str = parts[1].split("kcal")[0].strip()
                                energy = float(energy_str)
                                if best_energy is None or energy < best_energy:
                                    best_energy = energy
                        except (IndexError, ValueError):
                            pass
        except Exception as e:
            print(f"        ⚠️  Erro ao ler {dlg_file.name}: {e}")

        return best_energy

    def run_screening(self):
        """Executa o screening de todos ligantes contra todos receptores."""
        print(f"\n{'=' * 60}")
        print(f"MULTI-RECEPTOR VIRTUAL SCREENING")
        print(f"{'=' * 60}")
        print(f"🧬 Receptores: {len(self.receptors)}")
        print(f"💊 Ligantes: {len(self.ligands)}")
        print(f"🖥️  GPUs: {', '.join(self.gpus)}")
        print(f"🔄 Runs por ligante: {self.runs}")
        print(f"📦 Tamanho do batch: {self.batch_size}")

        total_combinations = len(self.receptors) * len(self.ligands)
        print(f"📊 Total de combinações: {total_combinations}")

        start_time = time.time()

        # Criar arquivos batch
        batch_files = self.create_batch_files()

        print(
            f"\n📁 Resultados serão salvos em: {self.results_dir}/[receptor]/[ligante]/"
        )
        print(f"\n🚀 Iniciando screening...")

        # Executar batches em paralelo
        with ThreadPoolExecutor(max_workers=len(self.gpus)) as executor:
            futures = []

            for batch_id, batch_file in enumerate(batch_files):
                gpu_id = self.gpus[batch_id % len(self.gpus)]
                future = executor.submit(self.dock_batch, batch_file, batch_id, gpu_id)
                futures.append(future)

            # Coletar resultados
            completed = 0
            for future in as_completed(futures):
                batch_results = future.result()
                self.results.extend(batch_results)
                completed += 1
                print(
                    f"   📊 Progresso: {completed}/{len(batch_files)} batches concluídos"
                )

        elapsed = time.time() - start_time

        print(f"\n✅ Screening concluído em {elapsed:.1f}s")
        print(f"⚡ Velocidade: {total_combinations / elapsed:.2f} combinações/s")

        # Salvar resultados
        self._save_results()

    def _save_results(self):
        """Salva os resultados em formato CSV."""
        if not self.results:
            print("⚠️  Nenhum resultado para salvar")
            return

        # Criar DataFrame
        df = pd.DataFrame(self.results)

        # Salvar resultado geral
        csv_path = self.output_dir / "all_results.csv"
        df.to_csv(csv_path, index=False)
        print(f"\n📁 Resultados salvos:")
        print(f"   ✓ {csv_path}")

        # Criar matriz de resultados (receptores x ligantes)
        print(f"\n📊 Criando matriz de resultados...")
        pivot_table = df.pivot_table(
            index="ligand", columns="receptor", values="best_energy", aggfunc="min"
        )

        matrix_csv = self.output_dir / "energy_matrix.csv"
        pivot_table.to_csv(matrix_csv)
        print(f"   ✓ {matrix_csv}")

        # Contar arquivos DLG gerados
        dlg_count = 0
        for receptor_name in self.receptors:
            receptor_dir = self.results_dir / receptor_name
            if receptor_dir.exists():
                for ligand_dir in receptor_dir.iterdir():
                    if ligand_dir.is_dir():
                        dlg_files = list(ligand_dir.glob("*.dlg"))
                        dlg_count += len(dlg_files)

        print(f"\n📄 Arquivos DLG gerados: {dlg_count}")
        print(f"📁 Localização: {self.results_dir}/[receptor]/[ligante]/")

        # Estatísticas por receptor
        print(f"\n📊 Resumo por Receptor:")
        print(f"{'Receptor':<15} {'Processados':<12} {'Melhor Energia':<15}")
        print(f"{'-' * 42}")

        for receptor_name in self.receptors:
            receptor_df = df[df["receptor"] == receptor_name]
            if not receptor_df.empty:
                count = len(receptor_df)
                best = receptor_df["best_energy"].min()
                print(f"{receptor_name:<15} {count:<12} {best:<15.2f}")

        # Encontrar melhores combinações globais
        df_sorted = df.sort_values("best_energy")

        print(f"\n🌟 TOP 10 MELHORES COMBINAÇÕES:")
        print(f"{'=' * 75}")
        print(
            f"{'Rank':<6} {'Receptor':<15} {'Ligante':<20} {'Energia (kcal/mol)':<20}"
        )
        print(f"{'-' * 75}")

        for idx, (_, row) in enumerate(df_sorted.head(10).iterrows(), 1):
            print(
                f"{idx:<6} {row['receptor']:<15} {row['ligand']:<20} {row['best_energy']:<20.2f}"
            )

        print(f"\n✅ Análise concluída!")
        print(f"   • Total de combinações processadas: {len(df)}")
        print(f"   • Melhor energia global: {df['best_energy'].min():.2f} kcal/mol")
        print(f"   • Média de energia: {df['best_energy'].mean():.2f} kcal/mol")


def main():
    """Função principal."""
    parser = argparse.ArgumentParser(
        description="Multi-Receptor Virtual Screening com AutoDock-GPU",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos:
  # Screening básico
  python autodock_gpu_screening.py \\
    --rfiles /path/to/receptors/ \\
    --lfiles /path/to/ligands/
  
  # Com múltiplas GPUs
  python autodock_gpu_screening.py \\
    --rfiles /path/to/receptors/ \\
    --lfiles /path/to/ligands/ \\
    --gpus 0,1,2,3 \\
    --runs 20
  
  # Com configurações customizadas
  python autodock_gpu_screening.py \\
    --rfiles /path/to/receptors/ \\
    --lfiles /path/to/ligands/ \\
    --gpus 0,1 \\
    --runs 10 \\
    --batch-size 100 \\
    --output-dir multi_screening_results
        """,
    )

    parser.add_argument(
        "--rfiles",
        required=True,
        help="Diretório com subpastas contendo receptores e mapas .fld",
    )

    parser.add_argument("--lfiles", required=True, help="Diretório com ligantes PDBQT")

    parser.add_argument(
        "-o",
        "--output-dir",
        default="multi_screening_results",
        help="Diretório de saída (padrão: multi_screening_results)",
    )

    parser.add_argument(
        "-g",
        "--gpus",
        default="0",
        help="GPUs a usar, separadas por vírgula (padrão: 0)",
    )

    parser.add_argument(
        "-r",
        "--runs",
        type=int,
        default=10,
        help="Número de runs por ligante (padrão: 10)",
    )

    parser.add_argument(
        "-b",
        "--batch-size",
        type=int,
        default=50,
        help="Ligantes por lote (padrão: 50)",
    )

    args = parser.parse_args()

    try:
        # Criar objeto de screening
        screener = MultiReceptorScreening(
            receptors_dir=args.rfiles,
            ligands_dir=args.lfiles,
            output_dir=args.output_dir,
            gpus=args.gpus,
            runs=args.runs,
            batch_size=args.batch_size,
        )

        # Executar screening
        screener.run_screening()

    except Exception as e:
        print(f"\n❌ Erro: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

