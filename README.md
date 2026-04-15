# IRI-AKI Pipeline

**Computational pipeline for the identification of molecular targets in ischemia-reperfusion acute kidney injury (IR-AKI) in kidney transplantation**

## Overview

This project implements a two-phase in silico pipeline to identify genes and pathways consistently altered in ischemia-reperfusion injury (IRI) during kidney transplantation, and to characterize their functional role through co-expression network analysis.

- **Phase 1** (Python): Per-dataset preprocessing and differential expression analysis
- **Phase 2** (R): Cross-dataset meta-analysis, pathway enrichment, independent validation, co-expression network inference, cross-platform integration, and edge set enrichment analysis

Phase 2 leverages [MUUMI](https://github.com/fhaive/muumi) (Inkala et al., BMC Bioinformatics, 2026), an R package for statistical and network-based meta-analysis for multi-omics data integration.

## Datasets

### Category A -- IRI vs Control (meta-analysis)

| Dataset   | Platform              | Samples | Comparison                      |
|-----------|-----------------------|---------|---------------------------------|
| GSE43974  | Microarray (Illumina) | 206     | Deceased donor T3 vs Living T1  |
| GSE142077 | RNA-seq (HiSeq 2500)  | 10      | Reperfusion vs Pre-ischemia     |
| GSE126805 | RNA-seq               | 83      | Baseline post vs Baseline pre   |

### Category B -- Validation (DGF vs non-DGF)

| Dataset  | Platform                 | Samples | Comparison            |
|----------|--------------------------|---------|-----------------------|
| GSE54888 | Microarray (Affymetrix)  | 54      | DGF vs non-DGF        |

All datasets were retrieved from NCBI Gene Expression Omnibus (GEO).

## Project structure

```
iri_aki_pipeline/
|-- notebooks/
|   |-- fase1_GSE43974_raw.ipynb        # Phase 1: GSE43974 preprocessing + DE
|   |-- fase1_GSE142077_raw.ipynb       # Phase 1: GSE142077 preprocessing + DE
|   |-- fase1_GSE126805_raw.ipynb       # Phase 1: GSE126805 preprocessing + DE
|   |-- fase1_GSE54888_raw.ipynb        # Phase 1: GSE54888 preprocessing + DE
|   |-- fase2_meta_analysis_MUUMI.ipynb # Phase 2: Meta-analysis + networks (R)
|
|-- output/
|   |-- fase1/
|   |   |-- GSE43974/                   # DE results, expression matrix, sample labels
|   |   |-- GSE142077/
|   |   |-- GSE126805/
|   |   |-- GSE54888/
|   |-- fase2/
|       |-- meta_analysis_ranked_gene_list.csv
|       |-- seed_genes.csv
|       |-- seed_genes_high_confidence.csv
|       |-- validation_catB_summary.csv
|       |-- dotplot_ORA_up.jpg
|       |-- dotplot_ORA_down.jpg
|       |-- volcano_meta_analysis.jpg
|       |-- heatmap_top_seed_genes.jpg
|       |-- platform_contributions.jpg
|       |-- bubbleplot_*.jpg
|       |-- network_*.rds
|       |-- pathway_*/
|       |-- esea_*_results.csv
|
|-- data/                               # Reference data (GMT files, etc.)
|-- requirements.txt                    # Python dependencies
|-- README.md
```

## Methods

### Phase 1 -- Preprocessing and differential expression

Each dataset is processed independently in Python:

- **Microarray** (GSE43974, GSE54888): Log2 transformation, quantile normalization, probe-to-gene mapping, Welch t-test with Benjamini-Hochberg FDR correction
- **RNA-seq** (GSE142077, GSE126805): Count matrix loading, differential expression with pyDESeq2, Benjamini-Hochberg correction

Output: one CSV per dataset with columns `gene, log2FoldChange, pvalue, padj, comparison, dataset`.

### Phase 2 -- Meta-analysis and network integration

Implemented in R using MUUMI (all four modules):

1. **Statistical meta-analysis** (Module 1): Ensemble method combining effect size and Fisher p-value approaches via `run_ensembl_metanalysis()`. Produces a consensus gene ranking across all Cat. A datasets.

2. **GSEA and seed gene selection**: Pre-ranked GSEA on Reactome pathways determines a biologically informed threshold. Seed genes are selected by: (i) rank within the GSEA threshold, (ii) |weighted log2FC| > 0.379 (FC >= 1.3), and (iii) direction consistency > 0.5 across datasets. A high-confidence subset requires consistency = 1.0 (concordant in all 3 datasets). Overrepresentation analysis (ORA) is performed separately on up- and down-regulated seed genes.

3. **Independent validation** (Cat. B): Seed genes are tested on the GSE54888 ranking (DGF vs non-DGF) via pre-ranked GSEA, confirming their relevance to injury severity.

4. **Co-expression network inference** (Module 2): Networks are inferred separately for IRI and Control conditions on two platforms (GSE43974 microarray, GSE126805 RNA-seq) using CLR with Pearson correlation. Community detection is performed with the walktrap algorithm, followed by Reactome pathway annotation.

5. **Cross-platform late integration** (Module 3): Communities from microarray and RNA-seq networks are aggregated using Non-negative Matrix Factorization (NMF) via `aggregate_communities()`, producing integrated meta-modules and quantifying each platform's contribution.

6. **Edge Set Enrichment Analysis** (Module 4): ESEA via `compute_esea()` evaluates functional enrichment at the edge level using KEGG pathway edge maps, identifying pathways with coordinated gene-gene interactions.

## Key results

| Metric | Value |
|--------|-------|
| Genes in meta-analysis | 12,127 |
| Seed genes (full set) | 650 (544 up, 106 down) |
| Seed genes (high-confidence) | 269 (229 up, 40 down) |
| Validation NES (full set) | 1.73 (p = 9.4 x 10^-11) |
| Validation NES (HC set) | 2.39 (p = 2.6 x 10^-18) |

Top enriched pathways in IRI seed genes include Toll-like receptor cascades (TLR2, TLR4), interleukin signaling (IL-4/IL-13, IL-10), interferon alpha/beta signaling, NGF-stimulated transcription, and regulated necrosis. Down-regulated seed genes are enriched in biological oxidations, peroxisomal protein import, and tight junction interactions.

Cross-platform late integration confirmed that IL-4/IL-13 signaling is the top enriched pathway in both microarray and RNA-seq networks independently, demonstrating reproducibility across platforms.

## Requirements

### Python (Phase 1)

Python >= 3.12. Install dependencies:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Key libraries: `pandas`, `numpy`, `scipy`, `statsmodels`, `GEOparse`, `pyDESeq2`.

### R (Phase 2)

R >= 4.4. Install dependencies:

```r
# Bioconductor
install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(c("clusterProfiler", "ReactomePA", "org.Hs.eg.db",
                        "enrichplot", "fgsea", "minet", "RankProd"))

# CRAN
install.packages(c("tidyverse", "metap", "esc", "igraph", "TopKLists",
                    "pamr", "foreach", "doParallel", "plyr", "gridExtra",
                    "abind", "stringr", "nmfbin", "SNFtool", "msigdbr",
                    "ggrepel", "reshape2"))

# MUUMI
remotes::install_github("fhaive/muumi")
```

## Reproducibility

1. Run Phase 1 notebooks in order (`fase1_GSE*.ipynb`) to generate per-dataset DE results in `output/fase1/`
2. Run the Phase 2 notebook (`fase2_meta_analysis_MUUMI.ipynb`) to perform the full meta-analysis, validation, and network analysis
3. All figures are saved as 300 DPI JPG in `output/fase2/`

## References

- Inkala S, Fratello M, del Giudice G, et al. MUUMI: an R package for statistical and network-based meta-analysis for multi-omics data integration. *BMC Bioinformatics*. 2026;27:56. doi:10.1186/s12859-026-06394-3

## License

This project is intended for academic research purposes.
