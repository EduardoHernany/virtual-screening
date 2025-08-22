#!/usr/bin/env python3
"""
Preparador de Macromoléculas para AutoDock
Prepara receptores e gera mapas de grid para todos os tipos de átomos

Uso:
    python prepare_macromolecules.py -r receptor.pdb --gsize 60 60 60 --gcenter 15.5 20.3 10.2
    python prepare_macromolecules.py -r receptor.pdb --gsize 60 60 60 --gcenter 15.5 20.3 10.2 -l ligand.pdb
"""

import argparse
import os
import sys
import subprocess
import shutil
from pathlib import Path
import time

from util import textfld


class MacromoleculePreparator:
    """Classe para preparar macromoléculas e gerar mapas de grid."""
    
    def __init__(self):
        """Inicializa os caminhos das ferramentas."""
        home = Path.home()
        
        # Caminhos das ferramentas MGLTools
        self.pythonsh_path = home / "mgltools_x86_64Linux2_1.5.7/bin/pythonsh"
        self.prepare_receptor_path = home / "mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
        self.prepare_ligand_path = home / "mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
        self.prepare_gpf_path = home / "mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py"
        
        # Caminhos do AutoDock
        self.ad4_parameters = home / "x86_64Linux2/AD4_parameters.dat"
        self.autogrid4_path = home / "x86_64Linux2/autogrid4"
        
        # Grupos de tipos de ligantes para cada grid
        self.ligand_groups = [
            "C,A,N,NA,NS,OA,OS,SA,S,H,HD",
            "HS,P,Br,BR,Ca,CA,Cl,CL,F,Fe,FE",
            "I,Mg,MG,Mn,MN,Zn,ZN,He,Li,Be",
            "B,Ne,Na,Al,Si,K,Sc,Ti,V,Co",
            "Ni,Cu,Ga,Ge,As,Se,Kr,Rb,Sr,Y",
            "Zr,Nb,Cr,Tc,Ru,Rh,Pd,Ag,Cd,In",
            "Sn,Sb,Te,Xe,Cs,Ba,La,Ce,Pr,Nd",
            "Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb",
            "Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg",
            "Tl,Pb,Bi,Po,At,Rn,Fr,Ra,Ac,Th",
            "Pa,U,Np,Pu,Am,Cm,Bk,Cf,E,Fm"
        ]
        
    def validate_tools(self):
        """Valida se todas as ferramentas necessárias estão instaladas."""
        tools = {
            "pythonsh": self.pythonsh_path,
            "prepare_receptor4": self.prepare_receptor_path,
            "prepare_ligand4": self.prepare_ligand_path,
            "prepare_gpf4": self.prepare_gpf_path,
            "AD4_parameters": self.ad4_parameters,
            "autogrid4": self.autogrid4_path
        }
        
        print("\n🔍 Verificando ferramentas...")
        missing = []
        
        for name, path in tools.items():
            if path.exists():
                print(f"   ✅ {name}: {path}")
            else:
                print(f"   ❌ {name}: NÃO ENCONTRADO em {path}")
                missing.append(name)
        
        if missing:
            print(f"\n❌ Ferramentas faltando: {', '.join(missing)}")
            print("\n📦 Instale o MGLTools e AutoDock4:")
            print("   1. Download MGLTools: http://mgltools.scripps.edu/downloads")
            print("   2. Download AutoDock4: http://autodock.scripps.edu/downloads")
            return False
        
        return True
    
    def prepare_receptor(self, receptor_path, output_dir):
        """
        Prepara o receptor convertendo PDB para PDBQT.
        
        Args:
            receptor_path: Caminho do arquivo PDB do receptor
            output_dir: Diretório de saída
            
        Returns:
            Path: Caminho do arquivo PDBQT gerado
        """
        receptor_name = receptor_path.stem
        pdbqt_file = output_dir / f"{receptor_name}.pdbqt"
        
        print(f"\n🧬 Preparando receptor: {receptor_name}")
        
        # Copiar arquivo PDB para o diretório de saída primeiro
        receptor_copy = output_dir / receptor_path.name
        shutil.copy2(receptor_path, receptor_copy)
        
        cmd = [
            str(self.pythonsh_path),
            str(self.prepare_receptor_path),
            "-r", receptor_copy.name,  # Usar apenas o nome do arquivo
            "-o", f"{receptor_name}.pdbqt",  # Nome do arquivo de saída
            "-A", "hydrogens"  # Adicionar hidrogênios
        ]
        
        try:
            print(f"   Executando: prepare_receptor4.py")
            # Executar no diretório de saída
            result = subprocess.run(cmd, cwd=output_dir, capture_output=True, text=True, check=True)
            
            if pdbqt_file.exists():
                print(f"   ✅ Receptor preparado: {pdbqt_file.name}")
                return pdbqt_file
            else:
                print(f"   ❌ Erro: arquivo PDBQT não foi gerado")
                if result.stderr:
                    print(f"      Erro: {result.stderr[:200]}")
                return None
                
        except subprocess.CalledProcessError as e:
            print(f"   ❌ Erro ao preparar receptor: {e}")
            if e.stderr:
                print(f"      Erro: {e.stderr[:200]}")
            return None
    
    def prepare_gpf_files(self, pdbqt_file, grid_center, grid_size, output_dir):
        """
        Gera arquivos GPF para todos os grupos de ligantes.
        
        Args:
            pdbqt_file: Arquivo PDBQT do receptor
            grid_center: Centro do grid (x, y, z)
            grid_size: Tamanho do grid (x, y, z)
            output_dir: Diretório de saída
            
        Returns:
            list: Lista de arquivos GPF gerados
        """
        gpf_files = []
        
        print(f"\n📐 Gerando arquivos GPF...")
        print(f"   Centro: {grid_center[0]}, {grid_center[1]}, {grid_center[2]}")
        print(f"   Tamanho: {grid_size[0]}, {grid_size[1]}, {grid_size[2]}")
        
        for i, ligand_types in enumerate(self.ligand_groups, 1):
            gpf_file = output_dir / f"grid_{i}.gpf"
            
            # Usar apenas o nome do arquivo PDBQT, não o caminho completo
            receptor_name = pdbqt_file.name
            
            cmd = [
                str(self.pythonsh_path),
                str(self.prepare_gpf_path),
                "-r", receptor_name,  # Usar apenas o nome do arquivo
                "-o", f"grid_{i}.gpf",  # Nome do arquivo de saída
                "-p", f"gridcenter={grid_center[0]},{grid_center[1]},{grid_center[2]}",
                "-p", f"npts={grid_size[0]},{grid_size[1]},{grid_size[2]}",
                "-p", f"ligand_types={ligand_types}"
            ]
            
            try:
                # Executar no diretório onde está o receptor.pdbqt
                result = subprocess.run(cmd, cwd=output_dir, capture_output=True, text=True, check=True)
                
                if gpf_file.exists():
                    # Adicionar linha de parameter_file no início do arquivo
                    self.add_parameter_file(gpf_file)
                    gpf_files.append(gpf_file)
                    print(f"   ✅ grid_{i}.gpf ({ligand_types[:20]}...)")
                else:
                    print(f"   ❌ Falha ao gerar grid_{i}.gpf")
                    
            except subprocess.CalledProcessError as e:
                print(f"   ❌ Erro ao gerar grid_{i}.gpf: {e}")
                if e.stderr:
                    print(f"      Erro: {e.stderr[:200]}")
        
        return gpf_files
    
    def add_parameter_file(self, gpf_file):
        """
        Adiciona a linha parameter_file no início do arquivo GPF.
        
        Args:
            gpf_file: Arquivo GPF a modificar
        """
        with open(gpf_file, 'r') as f:
            content = f.read()
        
        with open(gpf_file, 'w') as f:
            f.write(f"parameter_file {self.ad4_parameters}\n")
            f.write(content)
    
    def run_autogrid(self, gpf_files, output_dir):
        """
        Executa autogrid4 para gerar os mapas.
        
        Args:
            gpf_files: Lista de arquivos GPF
            output_dir: Diretório de trabalho
            
        Returns:
            Path: Arquivo .maps.fld gerado
        """
        print(f"\n🗺️  Gerando mapas de grid com autogrid4...")
        
        for i, gpf_file in enumerate(gpf_files, 1):
            log_file = output_dir / f"grid_{i}.glg"
            
            cmd = [
                str(self.autogrid4_path),
                "-p", str(gpf_file.name),
                "-l", str(log_file.name)
            ]
            
            try:
                print(f"   Processando grid_{i}.gpf...")
                result = subprocess.run(cmd, cwd=output_dir, 
                                      capture_output=True, text=True, check=True)
                print(f"   ✅ Grid {i} concluído")
                
            except subprocess.CalledProcessError as e:
                print(f"   ❌ Erro no grid {i}: {e}")
                if e.stderr:
                    print(f"      {e.stderr[:100]}")
        
        # Procurar arquivo .maps.fld gerado
        fld_files = list(output_dir.glob("*.maps.fld"))
        if fld_files:
            return fld_files[0]
        return None
    
    def process_fld_file(self, fld_file, receptor_name):
        """
        Processa o arquivo .maps.fld removendo linhas desnecessárias.
        
        Args:
            fld_file: Arquivo .maps.fld
            receptor_name: Nome do receptor
        """
        print(f"\n📝 Processando arquivo .maps.fld...")
        
        
        line_number = 23
        
        with open(fld_file, 'r') as file:
            lines = file.readlines()


        if len(lines) >= line_number:
            new_lines = lines[:line_number]
            with open(fld_file, 'w') as file:
                file.writelines(new_lines)
        
        texto = textfld()
        novo_texto = texto.replace("kakakakaka", receptor_name)
        
        # Salvar arquivo processado
        with open(fld_file, "a") as file:
            file.write(novo_texto)
        
        print(f"   ✅ Arquivo processado: {fld_file.name}")
    
    def calculate_ligand_center(self, ligand_path):
        """
        Calcula o centro do ligante para redocking.
        
        Args:
            ligand_path: Caminho do arquivo PDB do ligante
            
        Returns:
            tuple: Centro (x, y, z) ou None
        """
        print(f"\n📍 Calculando centro do ligante...")
        
        try:
            coords = []
            with open(ligand_path, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
            
            if coords:
                import numpy as np
                center = np.mean(coords, axis=0)
                print(f"   Centro calculado: {center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}")
                return tuple(center)
            else:
                print(f"   ❌ Nenhuma coordenada encontrada no ligante")
                return None
                
        except Exception as e:
            print(f"   ❌ Erro ao calcular centro: {e}")
            return None
    
    def run(self, receptor_path, grid_size, grid_center=None, ligand_path=None):
        """
        Executa o pipeline completo de preparação.
        
        Args:
            receptor_path: Caminho do receptor PDB
            grid_size: Tamanho do grid (x, y, z)
            grid_center: Centro do grid (x, y, z) ou None
            ligand_path: Caminho do ligante PDB (opcional, para redocking)
        """
        receptor_path = Path(receptor_path)
        receptor_name = receptor_path.stem
        
        # Criar diretório de saída
        output_dir = Path("macros_preparadas") / receptor_name
        output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"\n{'='*60}")
        print(f"PREPARAÇÃO DE MACROMOLÉCULA")
        print(f"{'='*60}")
        print(f"📁 Receptor: {receptor_path.name}")
        print(f"📂 Saída: {output_dir}")
        
        # Se ligante fornecido e sem centro, calcular centro do ligante
        if ligand_path and not grid_center:
            ligand_path = Path(ligand_path)
            if ligand_path.exists():
                grid_center = self.calculate_ligand_center(ligand_path)
                if not grid_center:
                    print("❌ Não foi possível calcular o centro do ligante")
                    return
            else:
                print(f"❌ Ligante não encontrado: {ligand_path}")
                return
        elif not grid_center:
            print("❌ Centro do grid não especificado e sem ligante para calcular")
            return
        
        # 1. Preparar receptor
        pdbqt_file = self.prepare_receptor(receptor_path, output_dir)
        if not pdbqt_file:
            print("❌ Falha na preparação do receptor")
            return
        
        # 2. Gerar arquivos GPF
        gpf_files = self.prepare_gpf_files(pdbqt_file, grid_center, grid_size, output_dir)
        if not gpf_files:
            print("❌ Falha na geração dos arquivos GPF")
            return
        
        # 3. Executar autogrid
        fld_file = self.run_autogrid(gpf_files, output_dir)
        if fld_file:
            # 4. Processar arquivo .maps.fld
            self.process_fld_file(fld_file, receptor_name)
            print(f"\n✅ Preparação concluída com sucesso!")
            print(f"📄 Arquivo de mapas: {fld_file}")
            print(f"📁 Todos os arquivos em: {output_dir}")
            
            # Se redocking, preparar ligante também
            if ligand_path:
                self.prepare_ligand_for_redocking(ligand_path, output_dir)
        else:
            print("❌ Falha na geração dos mapas")
    
    def prepare_ligand_for_redocking(self, ligand_path, output_dir):
        """
        Prepara o ligante para redocking.
        
        Args:
            ligand_path: Caminho do ligante PDB
            output_dir: Diretório de saída
        """
        ligand_name = ligand_path.stem
        ligand_pdbqt = output_dir / f"{ligand_name}.pdbqt"
        
        print(f"\n💊 Preparando ligante para redocking: {ligand_name}")
        
        # Copiar ligante PDB para o diretório de saída
        ligand_copy = output_dir / ligand_path.name
        if not ligand_copy.exists():
            shutil.copy2(ligand_path, ligand_copy)
        
        cmd = [
            str(self.pythonsh_path),
            str(self.prepare_ligand_path),
            "-l", ligand_copy.name,  # Usar apenas o nome do arquivo
            "-o", f"{ligand_name}.pdbqt"  # Nome do arquivo de saída
        ]
        
        try:
            # Executar no diretório de saída
            result = subprocess.run(cmd, cwd=output_dir, capture_output=True, text=True, check=True)
            
            if ligand_pdbqt.exists():
                print(f"   ✅ Ligante preparado: {ligand_pdbqt.name}")
            else:
                print(f"   ❌ Erro: arquivo PDBQT do ligante não foi gerado")
                if result.stderr:
                    print(f"      Erro: {result.stderr[:200]}")
                
        except subprocess.CalledProcessError as e:
            print(f"   ❌ Erro ao preparar ligante: {e}")
            if e.stderr:
                print(f"      Erro: {e.stderr[:200]}")


def main():
    """Função principal."""
    parser = argparse.ArgumentParser(
        description='Preparador de Macromoléculas para AutoDock',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos:
  # Preparar receptor com centro e tamanho especificados
  python prepare_macromolecules.py -r receptor.pdb --gsize 60 60 60 --gcenter 15.5 20.3 10.2
  
  # Preparar com redocking (calcula centro baseado no ligante)
  python prepare_macromolecules.py -r receptor.pdb --gsize 60 60 60 -l ligand.pdb
  
  # Preparar com centro e tamanho customizados
  python prepare_macromolecules.py -r protein.pdb --gsize 80 80 80 --gcenter 10 15 20
        """
    )
    
    parser.add_argument('-r', '--receptor', required=True,
                       help='Arquivo PDB do receptor')
    
    parser.add_argument('--gsize', nargs=3, type=int, required=True,
                       help='Tamanho do grid box (X Y Z)')
    
    parser.add_argument('--gcenter', nargs=3, type=float,
                       help='Centro do grid box (X Y Z)')
    
    parser.add_argument('-l', '--ligand',
                       help='Arquivo PDB do ligante (para redocking)')
    
    args = parser.parse_args()
    
    # Criar preparador
    preparator = MacromoleculePreparator()
    
    # Validar ferramentas
    if not preparator.validate_tools():
        sys.exit(1)
    
    # Executar preparação
    try:
        preparator.run(
            receptor_path=args.receptor,
            grid_size=args.gsize,
            grid_center=args.gcenter,
            ligand_path=args.ligand
        )
    except Exception as e:
        print(f"\n❌ Erro: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()