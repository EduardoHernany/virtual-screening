# 1. Instalar ambiente
conda env create -f environment.yml

~/miniconda3/bin/conda init fish

conda activate meeko-prep

-c "from rdkit import Chem; print(Chem.MolFromSmiles('CCO'))"

# 2. Preparar ligantes
python prepare_ligands.py ligands.sdf -o docking_ready

# 3. Os arquivos PDBQT estarão em ./docking_ready/

# Uso básico - processa um arquivo SDF
python prepare_ligands.py molecules.sdf

# Especificar diretório de saída
python prepare_ligands.py molecules.sdf --output-dir prepared_ligands

# Forma abreviada
python prepare_ligands.py molecules.sdf -o pdbqt_files