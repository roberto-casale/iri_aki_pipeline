# Magic Solution — In Silico Pipeline
## Progettazione computazionale di una soluzione di perfusione multi-target per IR-AKI

### Struttura del progetto

```
magic_solution/
├── config/
│   └── datasets.py          # Metadati di tutti i dataset GEO
├── src/
│   ├── download.py           # Download e caching dei dataset
│   ├── preprocess.py         # Preprocessing per piattaforma
│   ├── differential.py       # Analisi espressione differenziale
│   └── meta_analysis.py      # Meta-analisi cross-dataset
├── notebooks/
│   ├── 01_fase1_omics.py     # Fase 1: Meta-analisi omici (Jupyter-compatibile)
│   ├── 02_fase2_network.py   # Fase 2: Network medicine
│   ├── 03_fase3_drugs.py     # Fase 3: Drug repurposing
│   └── 04_fase4_silico.py    # Fase 4: Valutazione in silico
├── output/                   # Risultati generati
│   ├── fase1/
│   ├── fase2/
│   ├── fase3/
│   └── fase4/
├── data/                     # Dataset scaricati (non versionato)
│   └── geo/
├── requirements.txt
└── README.md
```

### Dataset utilizzati

| Accession   | Piattaforma            | Tipo       | Campioni | Ruolo        |
|-------------|------------------------|------------|----------|--------------|
| GSE43974    | Illumina HumanHT-12   | Microarray | 554      | Discovery    |
| GSE90861    | Illumina NextSeq 500  | RNA-seq    | 46       | Validazione  |
| GSE126805   | Illumina HiSeq 3000   | RNA-seq    | 163      | Validazione  |
| GSE142077   | Illumina HiSeq 2500   | RNA-seq    | 15       | Validazione  |
| GSE293480   | Illumina NovaSeq 6000 | RNA-seq    | 7        | Validazione  |
| GSE54888    | Affymetrix Gene 1.0ST | Microarray | 54       | Supporto     |
| GSE37838    | Affymetrix U133+2     | Microarray | 78       | Supporto     |
| GSE53769    | Affymetrix Gene 2.0ST | Microarray | 36       | Supporto     |

### Setup

```bash
# Crea ambiente virtuale
python -m venv .venv
source .venv/bin/activate  # Linux/Mac

# Installa dipendenze
pip install -r requirements.txt

# Per usare in Jupyter:
pip install ipykernel
python -m ipykernel install --user --name magic_solution
```

### Uso con Jupyter

I file in `notebooks/` usano il formato `# %%` (percent format).
Apri con Jupyter o VS Code/Cursor — ogni `# %%` diventa una cella.

### Uso con Cursor

Consigliato per lo sviluppo. Apri la cartella del progetto in Cursor,
usa i notebooks per esplorazione e i moduli in `src/` per codice riutilizzabile.
