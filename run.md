# Pipeline de Virtual Screening com AutoDock-GPU

## 📋 Índice

1. [Visão Geral](#visão-geral)
2. [prepare_ligands.py](#prepare_ligandspy)
   - [Descrição](#descrição)
   - [Instalação](#instalação)
   - [Parâmetros de Entrada](#parâmetros-de-entrada)
   - [Formatos Suportados](#formatos-suportados)
   - [Exemplos de Uso](#exemplos-de-uso)
   - [Estrutura de Saída](#estrutura-de-saída)
3. [autodock_gpu_screening.py](#autodock_gpu_screeningpy)
   - [Descrição](#descrição-1)
   - [Requisitos](#requisitos)
   - [Parâmetros de Entrada](#parâmetros-de-entrada-1)
   - [Estrutura de Entrada Esperada](#estrutura-de-entrada-esperada)
   - [Exemplos de Uso](#exemplos-de-uso-1)
   - [Estrutura de Saída](#estrutura-de-saída-1)
   - [Arquivos Gerados](#arquivos-gerados)
4. [Pipeline Completo](#pipeline-completo)
5. [Troubleshooting](#troubleshooting)

---

## Visão Geral

Este pipeline consiste em dois scripts principais para realizar virtual screening molecular usando AutoDock-GPU:

1. **`prepare_ligands.py`** - Prepara ligantes convertendo diversos formatos moleculares para PDBQT
2. **`autodock_gpu_screening.py`** - Executa docking molecular de múltiplos ligantes contra múltiplos receptores

---

## prepare_ligands.py

### Descrição

Script para preparação de ligantes para docking molecular. Converte arquivos moleculares (SDF, PDB, MOL, MOL2) para o formato PDBQT requerido pelo AutoDock.

### Instalação

```bash
# Criar ambiente conda
conda create -n meeko-prep python=3.9
conda activate meeko-prep

# Instalar dependências
conda install -c conda-forge rdkit openbabel pandas numpy
pip install meeko prody
```

### Parâmetros de Entrada

| Parâmetro | Obrigatório | Descrição | Valor Padrão |
|-----------|-------------|-----------|--------------|
| `input_file` | ✅ Sim | Arquivo de entrada (SDF, PDB, MOL, MOL2) | - |
| `-o, --output-dir` | ❌ Não | Diretório de saída para arquivos PDBQT | `{nome_arquivo}_pdbqt` |
| `--use-obabel` | ❌ Não | Usar OpenBabel ao invés de Meeko (mais rápido, menos preciso) | `False` (usa Meeko) |
| `--pH` | ❌ Não | pH para protonação | `7.4` |

### Formatos Suportados

| Formato | Extensão | Método Disponível | Moléculas |
|---------|----------|-------------------|-----------|
| SDF | `.sdf` | Meeko ✅ OpenBabel ✅ | Múltiplas |
| PDB | `.pdb` | Meeko ✅ | Única |
| MOL | `.mol`, `.mdl` | Meeko ✅ | Única |
| MOL2 | `.mol2` | Meeko ✅ | Única |

### Exemplos de Uso

#### Preparação com Meeko (Padrão - Mais Preciso)
```bash
# Arquivo SDF com múltiplas moléculas
python prepare_ligands.py molecules.sdf

# Arquivo PDB único
python prepare_ligands.py ligand.pdb -o prepared_ligands

# Com pH customizado
python prepare_ligands.py compound.mol --pH 6.5
```

#### Preparação com OpenBabel (Mais Rápido)
```bash
# Apenas para arquivos SDF
python prepare_ligands.py library.sdf --use-obabel -o quick_prep

# Com diretório de saída específico
python prepare_ligands.py compounds.sdf --use-obabel --output-dir screening_lib
```

### Estrutura de Saída

```
output_dir/
├── molecule1.pdbqt
├── molecule2.pdbqt
├── molecule3.pdbqt
└── ...
```

**Observações:**
- Nomes dos arquivos são preservados do arquivo original quando possível
- Moléculas sem nome recebem nomes sequenciais: `mol_0001.pdbqt`, `mol_0002.pdbqt`, etc.

---

## autodock_gpu_screening.py

### Descrição

Script para virtual screening em larga escala usando AutoDock-GPU. Executa docking de múltiplos ligantes contra múltiplos receptores em paralelo, utilizando uma ou mais GPUs.

### Requisitos

- **AutoDock-GPU** instalado em `~/AutoDock-GPU/bin/autodock_gpu_128wi`
- **CUDA** configurado e GPU(s) disponível(is)
- **Python 3.9+** com pandas

### Parâmetros de Entrada

| Parâmetro | Obrigatório | Descrição | Valor Padrão |
|-----------|-------------|-----------|--------------|
| `--rfiles` | ✅ Sim | Diretório com subpastas contendo receptores e mapas `.fld` | - |
| `--lfiles` | ✅ Sim | Diretório com ligantes PDBQT | - |
| `-o, --output-dir` | ❌ Não | Diretório de saída | `multi_screening_results` |
| `-g, --gpus` | ❌ Não | GPUs a usar (separadas por vírgula) | `0` |
| `-r, --runs` | ❌ Não | Número de runs por ligante | `10` |
| `-b, --batch-size` | ❌ Não | Número de combinações por batch | `50` |

### Estrutura de Entrada Esperada

#### Diretório de Receptores (`--rfiles`)
```
receptors_dir/
├── 1cjb/
│   ├── 1cjb_a.maps.fld    # Arquivo de mapas (obrigatório)
│   ├── 1cjb_a.pdbqt       # Receptor PDBQT (obrigatório)
│   ├── *.map              # Arquivos de mapa (obrigatório)
│   └── 1cjb_a.maps.xyz    # Coordenadas do grid (opcional)
├── 1g1g/
│   ├── 1g1g_a.maps.fld
│   ├── 1g1g_a.pdbqt
│   └── *.map
└── 2i7c/
    ├── 2i7c_a.maps.fld
    ├── 2i7c_a.pdbqt
    └── *.map
```

#### Diretório de Ligantes (`--lfiles`)
```
ligands_dir/
├── compound1.pdbqt
├── compound2.pdbqt
├── molecule1.pdbqt
└── ...
```

### Exemplos de Uso

#### Screening Básico
```bash
python autodock_gpu_screening.py \
    --rfiles /path/to/receptors/ \
    --lfiles /path/to/ligands/
```

#### Com Múltiplas GPUs
```bash
python autodock_gpu_screening.py \
    --rfiles /path/to/receptors/ \
    --lfiles /path/to/ligands/ \
    --gpus 0,1,2,3 \
    --runs 20 \
    --batch-size 100
```

#### Configuração Completa
```bash
python autodock_gpu_screening.py \
    --rfiles /home/user/receptors/ \
    --lfiles /home/user/ligands/ \
    --gpus 0,1 \
    --runs 10 \
    --batch-size 50 \
    --output-dir screening_results
```

### Estrutura de Saída

```
multi_screening_results/
├── all_results.csv          # Todos os resultados consolidados
├── energy_matrix.csv        # Matriz receptor × ligante
├── docking/                 # Arquivos de trabalho
│   ├── batch_0.txt         # Arquivos batch para AutoDock-GPU
│   ├── batch_1.txt
│   └── ...
├── logs/                    # Logs de execução
│   ├── batch_0_gpu_0.log
│   ├── batch_1_gpu_1.log
│   └── ...
└── results/                 # Resultados do docking
    ├── 1cjb_a/
    │   ├── compound1/
    │   │   ├── 1cjb_a_compound1.dlg    # Resultado do docking
    │   │   └── 1cjb_a_compound1.xml    # Arquivo XML (se gerado)
    │   ├── compound2/
    │   │   └── 1cjb_a_compound2.dlg
    │   └── ...
    ├── 1g1g_a/
    │   └── ...
    └── ...
```

### Arquivos Gerados

#### 1. **all_results.csv**
Arquivo CSV com todos os resultados de docking:

| Coluna | Descrição |
|--------|-----------|
| `receptor` | Nome do receptor |
| `ligand` | Nome do ligante |
| `best_energy` | Melhor energia de ligação (kcal/mol) |
| `dlg_file` | Caminho relativo do arquivo DLG |
| `batch_id` | ID do batch processado |

#### 2. **energy_matrix.csv**
Matriz pivotada com energias de ligação:

| ligand | 1cjb_a | 1g1g_a | 2i7c_a | ... |
|--------|--------|--------|--------|-----|
| compound1 | -8.5 | -7.2 | -9.1 | ... |
| compound2 | -7.8 | -8.3 | -7.5 | ... |
| ... | ... | ... | ... | ... |

#### 3. **Arquivos DLG**
Arquivos detalhados do AutoDock-GPU contendo:
- Todas as poses geradas
- Energias de cada run
- Coordenadas das melhores poses
- Informações de clustering

#### 4. **Logs**
Logs completos de cada batch contendo:
- Comando executado
- Output do AutoDock-GPU
- Erros (se houver)
- Tempo de execução

---

## Pipeline Completo

### Exemplo de Workflow Típico

```bash
# 1. Preparar ligantes de uma biblioteca SDF
python prepare_ligands.py compounds.sdf -o ligands_ready

# 2. Executar virtual screening
python autodock_gpu_screening.py \
    --rfiles /path/to/receptors/ \
    --lfiles ligands_ready/ \
    --gpus 0,1,2,3 \
    --runs 10 \
    --batch-size 100 \
    --output-dir screening_results

# 3. Analisar resultados
cd screening_results
head all_results.csv        # Ver primeiros resultados
python analyze_matrix.py energy_matrix.csv  # Script de análise customizado
```

### Workflow para Grande Escala

```bash
# 1. Preparação rápida com OpenBabel (para triagem inicial)
python prepare_ligands.py large_library.sdf --use-obabel -o quick_prep

# 2. Screening inicial rápido
python autodock_gpu_screening.py \
    --rfiles receptors/ \
    --lfiles quick_prep/ \
    --runs 5 \
    --batch-size 200 \
    --output-dir initial_screening

# 3. Selecionar top hits (exemplo: top 100)
python select_top_hits.py initial_screening/all_results.csv --top 100 > top_hits.sdf

# 4. Re-preparar top hits com Meeko (mais preciso)
python prepare_ligands.py top_hits.sdf -o refined_ligands

# 5. Re-docking com mais precisão
python autodock_gpu_screening.py \
    --rfiles receptors/ \
    --lfiles refined_ligands/ \
    --runs 50 \
    --batch-size 20 \
    --output-dir final_screening
```

---

## Troubleshooting

### Problemas Comuns e Soluções

#### 1. **Erro: AutoDock-GPU não encontrado**
```bash
# Verificar instalação
ls ~/AutoDock-GPU/bin/

# Se não existir, instalar AutoDock-GPU
git clone https://github.com/ccsb-scripps/AutoDock-GPU.git
cd AutoDock-GPU
make DEVICE=CUDA NUMWI=128
```

#### 2. **Erro: Arquivo PDB não pode ser lido**
```bash
# Tentar conversão com OpenBabel primeiro
obabel input.pdb -O output.sdf
python prepare_ligands.py output.sdf
```

#### 3. **Erro: GPU não disponível**
```bash
# Verificar GPUs
nvidia-smi

# Executar sem GPU (não recomendado - muito lento)
# Modificar o script para usar CPU
```

#### 4. **Erro: Sem resultados salvos**
```bash
# Verificar logs
cat multi_screening_results/logs/batch_0_gpu_0.log

# Verificar se arquivos DLG foram gerados
ls multi_screening_results/results/*/*/*.dlg

# Verificar arquivo batch
cat multi_screening_results/docking/batch_0.txt
```

#### 5. **Memória insuficiente na GPU**
```bash
# Reduzir batch size
python autodock_gpu_screening.py \
    --rfiles receptors/ \
    --lfiles ligands/ \
    --batch-size 10  # Menor batch size
```

### Verificação de Resultados

```bash
# Contar combinações processadas
wc -l all_results.csv

# Ver melhores energias
sort -t',' -k3 -n all_results.csv | head -20

# Verificar matriz
python -c "import pandas as pd; df=pd.read_csv('energy_matrix.csv'); print(df.describe())"
```

### Performance Esperada

| Configuração | Ligantes | Receptores | GPUs | Tempo Estimado |
|-------------|----------|------------|------|----------------|
| Pequena | 100 | 3 | 1 | ~15 min |
| Média | 1000 | 5 | 2 | ~2 horas |
| Grande | 10000 | 10 | 4 | ~24 horas |

**Fatores que afetam performance:**
- Tamanho dos receptores
- Número de átomos rotativos nos ligantes
- Número de runs (`--runs`)
- Velocidade da GPU
- Batch size

---

## 📧 Suporte

Para problemas ou sugestões, verificar:
- Logs em `logs/`
- Documentação do AutoDock-GPU
- Documentação do Meeko

---

*Última atualização: Agosto 2025*