import os
import shutil
import subprocess
from pathlib import Path
from datetime import datetime
import csv

from util import textfld

# Defina os diretórios base e de saída
base_dir = os.path.expanduser("~/preparacao/arquivos")
output_base_dir = os.path.expanduser("~/preparacao/arquivos_preparados")

pythonsh_path = os.path.expanduser("~/mgltools_x86_64Linux2_1.5.7/bin/pythonsh")
prep_gpf_path = os.path.expanduser("~/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py")
autogrid_path = os.path.expanduser("~/x86_64Linux2/autogrid4")
ad4_parameters_path = os.path.expanduser("~/x86_64Linux2/AD4_parameters.dat")
autodock_gpu_path = os.path.expanduser("~/AutoDock-GPU/bin/autodock_gpu_128wi")

ligand_types_list = [
    "ligand_types=C,A,N,NA,NS,OA,OS,SA,S,H,HD",
    "ligand_types=HS,P,Br,BR,Ca,CA,Cl,CL,F,Fe,FE",
    "ligand_types=I,Mg,MG,Mn,MN,Zn,ZN,He,Li,Be",
    "ligand_types=B,Ne,Na,Al,Si,K,Sc,Ti,V,Co",
    "ligand_types=Ni,Cu,Ga,Ge,As,Se,Kr,Rb,Sr,Y",
    "ligand_types=Zr,Nb,Cr,Tc,Ru,Rh,Pd,Ag,Cd,In",
    "ligand_types=Sn,Sb,Te,Xe,Cs,Ba,La,Ce,Pr,Nd",
    "ligand_types=Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb",
    "ligand_types=Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg",
    "ligand_types=Tl,Pb,Bi,Po,At,Rn,Fr,Ra,Ac,Th",
    "ligand_types=Pa,U,Np,Pu,Am,Cm,Bk,Cf,E,Fm"
]
def log_message(title='', message=''):
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] - [{title}] - {message}")

# Função para extrair o menor RMSD
def extrair_menor_rmsd(dpf_file_path):
    menor_rmsd = float('inf')
    energia_rmsd = None

    try:
        with open(dpf_file_path, 'r') as arquivo:
            for linha in arquivo:
                if 'RANKING' in linha:
                    dados = linha.split()
                    if len(dados) >= 6:
                        rmsd_value = float(dados[5])
                        binding_energy = float(dados[3])

                        if rmsd_value < menor_rmsd:
                            menor_rmsd = rmsd_value
                            energia_rmsd = binding_energy
    except FileNotFoundError:
        return None, None
    except Exception as e:
        log_message("Erro", f"Erro ao processar o arquivo .dpf: {str(e)}")
        return None, None

    return menor_rmsd, energia_rmsd

# Função para preparar o receptor
def preparar_receptor(receptor_path, ligant_filename, work_dir, output_dir, receptor):
    print("==================================================================================================================")
    log_message(receptor,f"Iniciando preparacao para receptor: {receptor} e ligante: {ligant_filename}")

    # Criar o diretório de saída se não existir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Ler os parâmetros do arquivo gridbox_param.txt e obter o nome da macromolécula
    param_file_path = os.path.join(work_dir, "gridbox_param.txt")
    with open(param_file_path, 'r') as param_file:
        lines = param_file.readlines()
        gridcenter = f'gridcenter={lines[0].strip()}'
        gridsize = f'npts={lines[1].strip()}'
        macromolecula_name = lines[2].strip()

    # Preparar os arquivos .gpf
    log_message(receptor,"Preparando arquivos .gpf")

    for i, ligand_types in enumerate(ligand_types_list, start=1):
        output_gpf = output_dir / f'gridbox{i}.gpf'
        
        command = [
            pythonsh_path, prep_gpf_path, "-r", receptor_path, "-o", str(output_gpf), 
            "-p", gridcenter, "-p", gridsize, "-p", ligand_types
        ]
        
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=work_dir)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            log_message(receptor,f"Erro ao executar o comando: {stderr.decode()}")
        
        
        with open(output_dir / 'comandos.txt', 'a') as log_file:
            log_file.write(" ".join(command) + "\n")

    log_message(receptor,"Arquivos .gpf preparados com sucesso")

    # Adicionar o parâmetro ao início dos arquivos .gpf e executar o AutoGrid
    log_message(receptor,"Adicionando parâmetros e executando AutoGrid")

    for i in range(1, len(ligand_types_list) + 1):
        gpf_file = output_dir / f'gridbox{i}.gpf'
        sed_command = f"sed -i '1i\\parameter_file {os.path.abspath(ad4_parameters_path)}' {gpf_file}"
        # sed_command = f"sed -i '1i\\parameter_file ad4_parameters' gridbox1.gpf"
        
        subprocess.run(sed_command, shell=True, check=True)
        
        autogrid_command = f"{autogrid_path} -p gridbox{i}.gpf -l gridbox{i}.glg"
        subprocess.run(autogrid_command, shell=True, check=True, cwd=output_dir)
        log_message(receptor,f"Grid {i+1} concluido.")
        
        
    log_message(receptor,"Parâmetros e execução do AutoGrid concluido com sucesso")

    # Manipular o arquivo receptor.maps.fld
    receptor_name = Path(receptor_path).stem
    fld_file_path = os.path.join(output_dir, f'{receptor_name}.maps.fld')
    line_number = 23

    with open(fld_file_path, 'r') as file:
        lines = file.readlines()

    log_message(receptor,"modificando arquivo .fld")

    if len(lines) >= line_number:
        new_lines = lines[:line_number]
        with open(fld_file_path, 'w') as file:
            file.writelines(new_lines)
        log_message(receptor,f"Linhas abaixo da linha {line_number} foram removidas do arquivo {fld_file_path}")

    # Substituir texto no arquivo fld
    texto = textfld()
    novo_texto = texto.replace("kakakakaka", receptor_name)

    with open(fld_file_path, "a") as file:
        file.write(novo_texto)

    log_message(receptor,"Arquivo .fld modificado com sucesso")
    
    fld_file_path = os.path.join(output_dir, f'{Path(receptor_path).stem}.maps.fld')
    autodock_command = [
        autodock_gpu_path, 
        '--ffile', str(fld_file_path), 
        '--lfile', str(ligant_filename),
        '--gbest', '1',
        '--nrun', "100"
    ]

    log_message(receptor,f"Executando AutoDock-GPU.")

    process = subprocess.Popen(autodock_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=output_dir)
    stdout, stderr = process.communicate()
    
    if process.returncode != 0:
        log_message(receptor,f"Erro ao executar AutoDock-GPU: {stderr.decode()}")
    else:
        log_message(receptor,f"AutoDock-GPU executado com sucesso.")

    # Extrair o menor RMSD e energia do arquivo .dpf
    dpf_file_path = output_dir / f'{Path(ligant_filename).stem}.dlg'
    log_message(receptor,f"dpf: {dpf_file_path}")
    
    menor_rmsd, energia_rmsd = extrair_menor_rmsd(dpf_file_path)
    log_message(receptor,f"Energia: {energia_rmsd}, RMSD: {menor_rmsd}")

    # Formatar o nome do rec_fld
    nome_receptor_sem_a = Path(receptor_path).stem.replace("_A", "")
    # rec_fld = f"macromoleculas/vivax/comRedocking/{nome_receptor_sem_a}/{receptor_path}.maps.fld"
    rec_fld = f"macromoleculas/falciparum/comRedocking/{fld_file_path.replace('/home/eduardo/preparacao/arquivos_preparados/','')}"
    log_message(receptor,f"fld: {rec_fld}")
    

    return {
        "id": None,  # Será preenchido sequencialmente
        "nome": macromolecula_name,
        "rec": Path(receptor_path).stem,
        "rec_fld": rec_fld,
        "ligante_original": Path(ligant_filename).stem,
        "rmsd_redoking": menor_rmsd,
        "energia_orinal": energia_rmsd,
        "gridsize": gridsize.replace("npts=",""),
        "gridcenter": gridcenter.replace("gridcenter=","")
    }

# Função para gerar CSV
def gerar_csv(dados, output_file):
    # Adicionar ID sequencial aos dados
    for idx, row in enumerate(dados, start=1):
        row['id'] = idx

    keys = dados[0].keys()
    with open(output_file, 'w', newline='') as output_csv:
        dict_writer = csv.DictWriter(output_csv, fieldnames=keys)
        dict_writer.writeheader()
        dict_writer.writerows(dados)

# Função principal
def processar_diretorios():
    log_message("Iniciando a preparacao e geração do CSV.")

    # Copiar todas as pastas da pasta "arquivos" para "arquivos_preparados"
    if os.path.exists(output_base_dir):
        shutil.rmtree(output_base_dir)  # Remover o diretório de saída se já existir
    shutil.copytree(base_dir, output_base_dir)

    resultados = []

    # Iterar sobre os diretórios em "arquivos_preparados" e listar os nomes dos receptores e ligantes
    for root, dirs, files in os.walk(output_base_dir):
        receptor_files = [f for f in files if f.endswith(("_A.pdbqt", "_a.pdbqt", "_b.pdbqt", "_ab.pdbqt", "_bd.pdbqt", "_A_ativo.pdbqt", "_A_alosterico.pdbqt"))]

        ligant_files = [f for f in files if f.endswith(".pdbqt") and f not in receptor_files and not f.endswith(".txt")]

        for receptor in receptor_files:
            receptor_path = os.path.join(root, receptor)
            work_dir = root
            output_dir = Path(root)

            for ligant in ligant_files:
                ligant_path = os.path.join(root, ligant)
                resultado = preparar_receptor(receptor_path, ligant_files[0], work_dir, output_dir, receptor)
                resultados.append(resultado)

    # Gerar o arquivo CSV final
    gerar_csv(resultados, os.path.expanduser("~/preparacao/resultados_macromoleculas.csv"))

    log_message("Processamento concluído e CSV gerado com sucesso.")

processar_diretorios()