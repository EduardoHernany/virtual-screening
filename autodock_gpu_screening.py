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
    
    def __init__(self, receptors_dir, ligands_dir, output_dir="screening_results", 
                 gpus="0", runs=10, batch_size=50):
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
        self.gpus = gpus.split(',')
        self.runs = runs
        self.batch_size = batch_size
        
        # Caminho do AutoDock-GPU
        self.autodock_gpu = Path.home() / "AutoDock-GPU/bin/autodock_gpu_128wi"
        
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
            raise FileNotFoundError(f"Diretório de receptores não encontrado: {self.receptors_dir}")
        
        if not self.ligands_dir.exists():
            raise FileNotFoundError(f"Diretório de ligantes não encontrado: {self.ligands_dir}")
        
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
                        receptor_name = fld_file.stem.replace('.maps', '')
                        receptors[receptor_name] = {
                            'fld': fld_file,
                            'dir': subdir,
                            'pdbqt': subdir / f"{receptor_name}.pdbqt",
                            'maps': list(subdir.glob("*.map"))
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
        
        # Criar subdiretórios para cada receptor
        for receptor_name in self.receptors:
            receptor_dir = self.docking_dir / receptor_name
            receptor_dir.mkdir(exist_ok=True)
    
    def prepare_receptor_files(self, receptor_name):
        """
        Copia os arquivos do receptor para o diretório de trabalho.
        
        Args:
            receptor_name: Nome do receptor
        
        Returns:
            Path: Diretório de trabalho do receptor
        """
        receptor_info = self.receptors[receptor_name]
        work_dir = self.docking_dir / receptor_name
        
        # Copiar arquivo .fld
        fld_dest = work_dir / receptor_info['fld'].name
        if not fld_dest.exists():
            shutil.copy(receptor_info['fld'], fld_dest)
        
        # Copiar arquivo .pdbqt se existir
        if receptor_info['pdbqt'].exists():
            pdbqt_dest = work_dir / receptor_info['pdbqt'].name
            if not pdbqt_dest.exists():
                shutil.copy(receptor_info['pdbqt'], pdbqt_dest)
        
        # Copiar todos os arquivos .map
        for map_file in receptor_info['maps']:
            map_dest = work_dir / map_file.name
            if not map_dest.exists():
                shutil.copy(map_file, map_dest)
        
        # Copiar arquivo .xyz se existir
        xyz_file = receptor_info['dir'] / receptor_info['fld'].name.replace('.fld', '.xyz')
        if xyz_file.exists():
            xyz_dest = work_dir / xyz_file.name
            if not xyz_dest.exists():
                shutil.copy(xyz_file, xyz_dest)
        
        return work_dir
    
    def create_batch_file(self, receptor_name, ligand_batch, batch_id):
        """
        Cria arquivo de lote para docking.
        
        Args:
            receptor_name: Nome do receptor
            ligand_batch: Lista de ligantes
            batch_id: ID do lote
        
        Returns:
            tuple: (batch_file_path, work_dir)
        """
        work_dir = self.docking_dir / receptor_name
        batch_file = work_dir / f"batch_{receptor_name}_{batch_id}.txt"
        
        receptor_info = self.receptors[receptor_name]
        
        with open(batch_file, 'w') as f:
            for ligand in ligand_batch:
                # Usar caminhos absolutos
                fld_path = str(receptor_info['fld'].absolute())
                ligand_path = str(ligand.absolute())
                
                f.write(f"{fld_path}\n")
                f.write(f"{ligand_path}\n")
        
        return batch_file, work_dir
    
    def dock_batch(self, receptor_name, ligand_batch, batch_id, gpu_id):
        """
        Executa docking de um lote.
        
        Args:
            receptor_name: Nome do receptor
            ligand_batch: Lista de ligantes
            batch_id: ID do lote
            gpu_id: ID da GPU
        
        Returns:
            list: Resultados do docking
        """
        batch_results = []
        
        # Preparar arquivos do receptor
        work_dir = self.prepare_receptor_files(receptor_name)
        
        # Criar arquivo de lote
        batch_file, work_dir = self.create_batch_file(receptor_name, ligand_batch, batch_id)
        
        # Comando do AutoDock-GPU
        cmd = [
            str(self.autodock_gpu),
            "--filelist", str(batch_file.absolute()),
            "--nrun", str(self.runs),
            "--lsmet", "ad",
            # "--devnum", str(gpu_id),
            "--nev", "2500000",
            "--ngen", "42000",
            "--lsit", "300",
            "--resnam", f"{receptor_name}_batch_{batch_id}",
        ]
        
        # Log file
        log_file = self.logs_dir / f"{receptor_name}_batch_{batch_id}_gpu_{gpu_id}.log"
        
        try:
            print(f"  🚀 [{receptor_name}] Lote {batch_id} na GPU {gpu_id} ({len(ligand_batch)} ligantes)")
            
            start_time = time.time()
            
            # Executar docking
            result = subprocess.run(cmd, cwd=work_dir, 
                                  capture_output=True, text=True, check=False)
            
            # Salvar log
            with open(log_file, 'w') as log:
                log.write(f"Receptor: {receptor_name}\n")
                log.write(f"Working directory: {work_dir}\n")
                log.write(f"Command: {' '.join(cmd)}\n")
                log.write(f"Return code: {result.returncode}\n")
                log.write(f"\n=== STDOUT ===\n{result.stdout}\n")
                log.write(f"\n=== STDERR ===\n{result.stderr}\n")
            
            if result.returncode != 0:
                print(f"  ❌ [{receptor_name}] Erro no lote {batch_id}")
                print(f"     📄 Log: {log_file}")
                
                # Mostrar erro
                error_output = result.stderr if result.stderr else result.stdout
                if error_output:
                    error_lines = error_output.split('\n')[:3]
                    for line in error_lines:
                        if line.strip():
                            print(f"     ⚠️  {line.strip()}")
            else:
                elapsed = time.time() - start_time
                print(f"  ✅ [{receptor_name}] Lote {batch_id} concluído em {elapsed:.1f}s")
                
                # Processar resultados
                batch_results = self._parse_results(receptor_name, batch_id, ligand_batch, work_dir)
            
        except Exception as e:
            print(f"  ❌ [{receptor_name}] Erro no lote {batch_id}: {e}")
        
        return batch_results
    
    def _parse_results(self, receptor_name, batch_id, ligand_batch, work_dir):
        """
        Extrai resultados do arquivo DLG.
        
        Args:
            receptor_name: Nome do receptor
            batch_id: ID do lote
            ligand_batch: Lista de ligantes
            work_dir: Diretório de trabalho
        
        Returns:
            list: Resultados parseados
        """
        results = []
        dlg_file = work_dir / f"{receptor_name}_batch_{batch_id}.dlg"
        
        if not dlg_file.exists():
            # Tentar encontrar qualquer arquivo .dlg
            dlg_files = list(work_dir.glob("*.dlg"))
            if dlg_files:
                dlg_file = dlg_files[-1]  # Pegar o mais recente
            else:
                return results
        
        current_ligand = None
        best_energy = float('inf')
        
        try:
            with open(dlg_file, 'r') as f:
                for line in f:
                    if "Ligand:" in line:
                        # Salvar resultado anterior
                        if current_ligand and best_energy < float('inf'):
                            ligand_name = Path(current_ligand).stem
                            results.append({
                                'receptor': receptor_name,
                                'ligand': ligand_name,
                                'best_energy': best_energy,
                                'batch_id': batch_id
                            })
                        
                        # Novo ligante
                        current_ligand = line.split()[-1]
                        best_energy = float('inf')
                    
                    elif "Estimated Free Energy of Binding" in line:
                        try:
                            # Formato: = -X.XX kcal/mol
                            parts = line.split('=')
                            if len(parts) > 1:
                                energy_str = parts[1].split('kcal')[0].strip()
                                energy = float(energy_str)
                                best_energy = min(best_energy, energy)
                        except (IndexError, ValueError):
                            pass
            
            # Salvar último resultado
            if current_ligand and best_energy < float('inf'):
                ligand_name = Path(current_ligand).stem
                results.append({
                    'receptor': receptor_name,
                    'ligand': ligand_name,
                    'best_energy': best_energy,
                    'batch_id': batch_id
                })
        
        except Exception as e:
            print(f"     ⚠️  Erro ao processar {dlg_file}: {e}")
        
        return results
    
    def run_screening(self):
        """Executa o screening de todos ligantes contra todos receptores."""
        print(f"\n{'='*60}")
        print(f"MULTI-RECEPTOR VIRTUAL SCREENING")
        print(f"{'='*60}")
        print(f"🧬 Receptores: {len(self.receptors)}")
        print(f"💊 Ligantes: {len(self.ligands)}")
        print(f"🖥️  GPUs: {', '.join(self.gpus)}")
        print(f"🔄 Runs por ligante: {self.runs}")
        print(f"📦 Tamanho do lote: {self.batch_size}")
        
        total_combinations = len(self.receptors) * len(self.ligands)
        print(f"📊 Total de combinações: {total_combinations}")
        
        start_time = time.time()
        
        # Criar tarefas para cada receptor
        all_tasks = []
        
        for receptor_name in self.receptors:
            print(f"\n🧬 Preparando tarefas para {receptor_name}...")
            
            # Dividir ligantes em lotes
            ligand_batches = [self.ligands[i:i+self.batch_size] 
                            for i in range(0, len(self.ligands), self.batch_size)]
            
            # Criar tarefas
            for batch_id, batch in enumerate(ligand_batches):
                gpu_id = self.gpus[len(all_tasks) % len(self.gpus)]
                all_tasks.append((receptor_name, batch, batch_id, gpu_id))
        
        print(f"\n📋 Total de tarefas: {len(all_tasks)}")
        print(f"\n🚀 Iniciando screening paralelo...")
        
        # Executar tarefas em paralelo
        with ThreadPoolExecutor(max_workers=len(self.gpus)) as executor:
            futures = []
            
            for receptor_name, batch, batch_id, gpu_id in all_tasks:
                future = executor.submit(self.dock_batch, receptor_name, batch, batch_id, gpu_id)
                futures.append(future)
            
            # Coletar resultados
            for future in as_completed(futures):
                batch_results = future.result()
                self.results.extend(batch_results)
        
        elapsed = time.time() - start_time
        
        print(f"\n✅ Screening concluído em {elapsed:.1f}s")
        print(f"⚡ Velocidade: {total_combinations/elapsed:.2f} combinações/s")
        
        # Salvar resultados
        self._save_results()
    
    def _save_results(self):
        """Salva os resultados em diferentes formatos."""
        if not self.results:
            print("⚠️  Nenhum resultado para salvar")
            return
        
        # Criar DataFrame
        df = pd.DataFrame(self.results)
        
        # Salvar resultado geral
        csv_path = self.output_dir / "all_results.csv"
        df.to_csv(csv_path, index=False)
        print(f"\n📁 Resultados gerais salvos em: {csv_path}")
        
        # Salvar resultados por receptor
        for receptor_name in self.receptors:
            receptor_df = df[df['receptor'] == receptor_name].copy()
            if not receptor_df.empty:
                receptor_df = receptor_df.sort_values('best_energy')
                
                # Salvar CSV
                receptor_csv = self.results_dir / f"{receptor_name}_results.csv"
                receptor_df.to_csv(receptor_csv, index=False)
                
                # Mostrar top 5
                print(f"\n🏆 Top 5 para {receptor_name}:")
                print(f"{'Rank':<6} {'Ligante':<20} {'Energia (kcal/mol)':<15}")
                print(f"{'-'*45}")
                
                for idx, (_, row) in enumerate(receptor_df.head(5).iterrows(), 1):
                    print(f"{idx:<6} {row['ligand']:<20} {row['best_energy']:<15.2f}")
        
        # Criar matriz de resultados (receptores x ligantes)
        print(f"\n📊 Criando matriz de resultados...")
        pivot_table = df.pivot_table(
            index='ligand',
            columns='receptor',
            values='best_energy',
            aggfunc='min'
        )
        
        matrix_csv = self.output_dir / "energy_matrix.csv"
        pivot_table.to_csv(matrix_csv)
        print(f"   ✓ Matriz salva em: {matrix_csv}")
        
        # Encontrar melhores combinações globais
        df_sorted = df.sort_values('best_energy')
        
        print(f"\n🌟 TOP 10 MELHORES COMBINAÇÕES GLOBAIS:")
        print(f"{'='*60}")
        print(f"{'Rank':<6} {'Receptor':<15} {'Ligante':<20} {'Energia':<15}")
        print(f"{'-'*60}")
        
        for idx, (_, row) in enumerate(df_sorted.head(10).iterrows(), 1):
            print(f"{idx:<6} {row['receptor']:<15} {row['ligand']:<20} {row['best_energy']:<15.2f}")


def main():
    """Função principal."""
    parser = argparse.ArgumentParser(
        description='Multi-Receptor Virtual Screening com AutoDock-GPU',
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
        """
    )
    
    parser.add_argument('--rfiles', required=True,
                       help='Diretório com subpastas contendo receptores e mapas .fld')
    
    parser.add_argument('--lfiles', required=True,
                       help='Diretório com ligantes PDBQT')
    
    parser.add_argument('-o', '--output-dir', default='multi_screening_results',
                       help='Diretório de saída (padrão: multi_screening_results)')
    
    parser.add_argument('-g', '--gpus', default='0',
                       help='GPUs a usar, separadas por vírgula (padrão: 0)')
    
    parser.add_argument('-r', '--runs', type=int, default=10,
                       help='Número de runs por ligante (padrão: 10)')
    
    parser.add_argument('-b', '--batch-size', type=int, default=50,
                       help='Ligantes por lote (padrão: 50)')
    
    args = parser.parse_args()
    
    try:
        # Criar objeto de screening
        screener = MultiReceptorScreening(
            receptors_dir=args.rfiles,
            ligands_dir=args.lfiles,
            output_dir=args.output_dir,
            gpus=args.gpus,
            runs=args.runs,
            batch_size=args.batch_size
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