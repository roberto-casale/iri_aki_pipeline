# %% [markdown]
# # Fase 1 — Meta-analisi dei Dati Omici (IR-AKI)
#
# Pipeline per identificare geni e pathway alterati nel danno renale
# da ischemia-riperfusione, integrando 8 dataset GEO umani.
#
# **Struttura:**
# 1. Download dataset (con cache locale)
# 2. Preprocessing per piattaforma
# 3. Espressione differenziale per dataset
# 4. Meta-analisi cross-dataset
# 5. Pathway enrichment (GSEA + ORA)

# %% Setup e imports
import os
import sys
import warnings
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
from collections import defaultdict

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(".").resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

DATA_DIR = PROJECT_ROOT / "data" / "geo"
OUTPUT_DIR = PROJECT_ROOT / "output" / "fase1"
DATA_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"Project root: {PROJECT_ROOT}")
print(f"Data dir: {DATA_DIR}")
print(f"Output dir: {OUTPUT_DIR}")

import GEOparse


# %% Utility: cache helper
def _cache_path(accession: str, kind: str) -> Path:
    """Restituisce il path del file pickle per un dataset."""
    return DATA_DIR / accession / f"{kind}.pkl"


def _is_cached(accession: str, *kinds: str) -> bool:
    """Controlla se tutti i file pickle richiesti esistono."""
    return all(_cache_path(accession, k).exists() for k in kinds)


def _load_cached(accession: str, kind: str) -> pd.DataFrame:
    """Carica un dataframe dal pickle cache."""
    path = _cache_path(accession, kind)
    print(f"  [CACHE] Loading {kind} from {path.name}")
    return pd.read_pickle(path)


def _save_cache(df: pd.DataFrame, accession: str, kind: str):
    """Salva un dataframe come pickle cache."""
    path = _cache_path(accession, kind)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_pickle(path)


# %% [markdown]
# ## 1. Download dei Dataset
#
# Ogni dataset viene scaricato solo la prima volta. Le esecuzioni successive
# caricano i dati dalla cache locale (file .pkl in data/geo/).

# %% 1.1 — GSE43974 (Discovery, Illumina HumanHT-12, n=554)
def download_GSE43974():
    acc = "GSE43974"

    if _is_cached(acc, "expression", "metadata"):
        print(f"[{acc}] Loading from cache...")
        return _load_cached(acc, "expression"), _load_cached(acc, "metadata")

    print(f"[{acc}] Downloading Series Matrix...")
    dest = DATA_DIR / acc
    dest.mkdir(exist_ok=True)
    gse = GEOparse.get_GEO(geo=acc, destdir=str(dest), silent=True)

    expr = gse.pivot_samples("VALUE")
    meta = gse.phenotype_data

    meta["donor_type"] = meta["title"].str.extract(r"Kidney_donor_(\w+)_T")
    meta["timepoint"] = meta["title"].str.extract(r"_T(\d)_")
    meta["condition"] = meta["timepoint"].map({
        "1": "control",
        "2": "cold_isch",
        "3": "IRI",
    })

    print(f"  Expression: {expr.shape[0]} probes x {expr.shape[1]} samples")
    print(f"  Donor types: {meta['donor_type'].value_counts().to_dict()}")
    print(f"  Timepoints: {meta['timepoint'].value_counts().to_dict()}")

    _save_cache(expr, acc, "expression")
    _save_cache(meta, acc, "metadata")
    return expr, meta

expr_43974, meta_43974 = download_GSE43974()


# %% 1.2 — GSE90861 (Validation, RNA-seq, n=46)
def download_GSE90861():
    acc = "GSE90861"

    if _is_cached(acc, "metadata"):
        meta = _load_cached(acc, "metadata")
        print(f"[{acc}] Loading from cache... ({len(meta)} samples)")
        return meta

    print(f"[{acc}] Downloading metadata...")
    dest = DATA_DIR / acc
    dest.mkdir(exist_ok=True)
    gse = GEOparse.get_GEO(geo=acc, destdir=str(dest), silent=True)
    meta = gse.phenotype_data

    meta["outcome"] = meta["title"].str.extract(r"^(IGF|DGF)")
    meta["biopsy_time"] = meta["title"].apply(
        lambda t: "post" if "b1" in t else "pre"
    )
    meta["patient_id"] = meta["title"].str.extract(r"_(\d+)")
    meta["condition"] = meta["biopsy_time"].map({"pre": "control", "post": "IRI"})

    print(f"  Outcomes: {meta['outcome'].value_counts().to_dict()}")
    print(f"  Biopsy times: {meta['biopsy_time'].value_counts().to_dict()}")

    counts_url = (
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE90861"
        "&format=file&file=GSE90861_RAW.tar"
    )
    print(f"  [INFO] Scarica manualmente i counts da: {counts_url}")
    print(f"  [INFO] Dopo il download, i file .txt nel tar contengono counts per campione")

    raw_tar = dest / "GSE90861_RAW.tar"
    if raw_tar.exists() and not _is_cached(acc, "counts"):
        import tarfile
        with tarfile.open(raw_tar) as tar:
            tar.extractall(dest / "raw")
        count_files = list((dest / "raw").glob("*.txt*"))
        counts_list = []
        for f in count_files:
            comp = "gzip" if str(f).endswith(".gz") else None
            df = pd.read_csv(f, sep="\t", index_col=0, compression=comp)
            counts_list.append(df)
        if counts_list:
            counts = pd.concat(counts_list, axis=1)
            _save_cache(counts, acc, "counts")
            print(f"  Counts matrix: {counts.shape}")

    _save_cache(meta, acc, "metadata")
    return meta

meta_90861 = download_GSE90861()


# %% 1.3 — GSE126805 (Validation, RNA-seq, n=163)
def download_GSE126805():
    acc = "GSE126805"

    if _is_cached(acc, "metadata"):
        meta = _load_cached(acc, "metadata")
        print(f"[{acc}] Loading metadata from cache... ({len(meta)} samples)")
        # Controlla anche counts
        counts_file = DATA_DIR / acc / "GSE126805_all.gene_counts.txt.gz"
        if counts_file.exists() and not _is_cached(acc, "counts"):
            counts = pd.read_csv(counts_file, sep="\t", index_col=0, compression="gzip")
            _save_cache(counts, acc, "counts")
            print(f"  Counts matrix loaded: {counts.shape}")
        elif _is_cached(acc, "counts"):
            print(f"  Counts also cached.")
        return meta

    print(f"[{acc}] Downloading...")
    dest = DATA_DIR / acc
    dest.mkdir(exist_ok=True)
    gse = GEOparse.get_GEO(geo=acc, destdir=str(dest), silent=True)
    meta = gse.phenotype_data

    counts_file = dest / "GSE126805_all.gene_counts.txt.gz"
    if counts_file.exists():
        counts = pd.read_csv(counts_file, sep="\t", index_col=0, compression="gzip")
        _save_cache(counts, acc, "counts")
        print(f"  Counts matrix: {counts.shape}")
    else:
        counts_url = (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126805/suppl/"
            "GSE126805_all.gene_counts.txt.gz"
        )
        print(f"  [INFO] Scarica il counts file da: {counts_url}")
        print(f"  [INFO] Salvalo in: {counts_file}")

    meta["patient_id"] = meta["title"].str.extract(r"kidney_(\d+)")
    meta["timepoint"] = meta["title"].apply(lambda t:
        "baseline_pre" if "baseline_pre" in t else
        "baseline_post" if "baseline_post" in t else
        "3months" if "3month" in t else
        "1year" if "1year" in t else "unknown"
    )
    meta["condition"] = meta["timepoint"].map({
        "baseline_pre": "control",
        "baseline_post": "IRI",
        "3months": "recovery_3m",
        "1year": "recovery_12m",
    })

    print(f"  Timepoints: {meta['timepoint'].value_counts().to_dict()}")
    _save_cache(meta, acc, "metadata")
    return meta

meta_126805 = download_GSE126805()


# %% 1.4 — GSE142077 (Validation, RNA-seq, n=15)
def download_GSE142077():
    acc = "GSE142077"

    if _is_cached(acc, "metadata"):
        meta = _load_cached(acc, "metadata")
        print(f"[{acc}] Loading from cache... ({len(meta)} samples)")
        return meta

    print(f"[{acc}] Downloading...")
    dest = DATA_DIR / acc
    dest.mkdir(exist_ok=True)
    gse = GEOparse.get_GEO(geo=acc, destdir=str(dest), silent=True)
    meta = gse.phenotype_data

    meta["condition_raw"] = meta["title"].str.extract(r"^(\w+)-IRI")
    meta["patient_id"] = meta["title"].str.extract(r"IRI-h-(\d+)")
    meta["condition"] = meta["condition_raw"].map({
        "P": "control",
        "I": "ischemia",
        "R10": "IRI",
    })

    print(f"  Conditions: {meta['condition'].value_counts().to_dict()}")
    print(f"  [INFO] RAW.tar da: https://www.ncbi.nlm.nih.gov/geo/download/?acc={acc}&format=file")

    _save_cache(meta, acc, "metadata")
    return meta

meta_142077 = download_GSE142077()


# %% 1.5 — GSE293480 (Validation, RNA-seq, n=7)
def download_GSE293480():
    acc = "GSE293480"

    if _is_cached(acc, "metadata"):
        meta = _load_cached(acc, "metadata")
        print(f"[{acc}] Loading from cache... ({len(meta)} samples)")
        # Controlla counts
        counts_file = DATA_DIR / acc / "GSE293480_rawCounts.csv.gz"
        if counts_file.exists() and not _is_cached(acc, "counts"):
            counts = pd.read_csv(counts_file, index_col=0, compression="gzip")
            _save_cache(counts, acc, "counts")
            print(f"  Counts loaded: {counts.shape}")
        return meta

    print(f"[{acc}] Downloading...")
    dest = DATA_DIR / acc
    dest.mkdir(exist_ok=True)
    gse = GEOparse.get_GEO(geo=acc, destdir=str(dest), silent=True)
    meta = gse.phenotype_data

    meta["outcome"] = meta["title"].str.extract(r"^(PDGF|IGF)")
    meta["condition"] = meta["outcome"].map({"PDGF": "IRI", "IGF": "control"})

    counts_file = dest / "GSE293480_rawCounts.csv.gz"
    if counts_file.exists():
        counts = pd.read_csv(counts_file, index_col=0, compression="gzip")
        _save_cache(counts, acc, "counts")
        print(f"  Counts: {counts.shape}")
    else:
        print("  [INFO] Scarica rawCounts.csv.gz dal supplementary GEO")

    _save_cache(meta, acc, "metadata")
    return meta

meta_293480 = download_GSE293480()


# %% 1.6 — GSE54888 (Support, Affymetrix Gene 1.0ST, n=54) — FIXED
def download_GSE54888():
    acc = "GSE54888"

    if _is_cached(acc, "expression", "metadata"):
        print(f"[{acc}] Loading from cache...")
        return _load_cached(acc, "expression"), _load_cached(acc, "metadata")

    print(f"[{acc}] Downloading...")
    dest = DATA_DIR / acc
    dest.mkdir(exist_ok=True)
    gse = GEOparse.get_GEO(geo=acc, destdir=str(dest), silent=True)

    expr = gse.pivot_samples("VALUE")
    meta = gse.phenotype_data

    # FIX: La colonna 'characteristics_ch1.2.dgf occurrence' ha valori 'DGF' / 'without DGF'
    # Usare confronto ESATTO, non substring match
    dgf_col = [c for c in meta.columns if "dgf occurrence" in c.lower()]
    if dgf_col:
        col = dgf_col[0]
        meta["condition"] = meta[col].apply(
            lambda x: "IRI" if str(x).strip() == "DGF" else "control"
        )
    else:
        # Fallback generico
        meta["condition"] = "unknown"
        for col in meta.columns:
            if "dgf" in col.lower():
                uniq = meta[col].unique()
                if len(uniq) >= 2 and len(uniq) <= 5:
                    meta["condition"] = meta[col].apply(
                        lambda x: "IRI" if str(x).strip() == "DGF" else "control"
                    )
                    break

    print(f"  Expression: {expr.shape}")
    print(f"  Conditions: {meta['condition'].value_counts().to_dict()}")

    _save_cache(expr, acc, "expression")
    _save_cache(meta, acc, "metadata")
    return expr, meta

expr_54888, meta_54888 = download_GSE54888()


# %% 1.7 — GSE37838 (Support, Affymetrix U133+2, n=78) — FIXED
def download_GSE37838():
    acc = "GSE37838"

    if _is_cached(acc, "expression", "metadata"):
        print(f"[{acc}] Loading from cache...")
        return _load_cached(acc, "expression"), _load_cached(acc, "metadata")

    print(f"[{acc}] Downloading...")
    dest = DATA_DIR / acc
    dest.mkdir(exist_ok=True)
    gse = GEOparse.get_GEO(geo=acc, destdir=str(dest), silent=True)

    expr = gse.pivot_samples("VALUE")
    meta = gse.phenotype_data

    # FIX: Cerca 'dgf status' (IGF/DGF/nephrectomy) o 'diagnosis' (good/bad)
    dgf_col = [c for c in meta.columns if "dgf status" in c.lower()]
    diag_col = [c for c in meta.columns if "diagnosis" in c.lower()]

    if dgf_col:
        col = dgf_col[0]
        print(f"  Using column: '{col}' -> {meta[col].value_counts().to_dict()}")
        meta["condition"] = meta[col].apply(
            lambda x: "IRI" if str(x).strip() == "DGF"
            else "control" if str(x).strip() in ("IGF", "nephrectomy")
            else "exclude"
        )
    elif diag_col:
        col = diag_col[0]
        print(f"  Using column: '{col}' -> {meta[col].value_counts().to_dict()}")
        meta["condition"] = meta[col].apply(
            lambda x: "IRI" if str(x).strip() == "bad"
            else "control" if str(x).strip() in ("good", "nephrectomy")
            else "exclude"
        )
    else:
        print("  [WARN] No suitable condition column found")
        meta["condition"] = "unknown"

    meta = meta[meta["condition"] != "exclude"]

    print(f"  Expression: {expr.shape}")
    print(f"  Conditions: {meta['condition'].value_counts().to_dict()}")

    _save_cache(expr, acc, "expression")
    _save_cache(meta, acc, "metadata")
    return expr, meta

expr_37838, meta_37838 = download_GSE37838()


# %% 1.8 — GSE53769 (Support, Affymetrix Gene 2.0ST, n=36)
def download_GSE53769():
    acc = "GSE53769"

    if _is_cached(acc, "expression", "metadata"):
        print(f"[{acc}] Loading from cache...")
        return _load_cached(acc, "expression"), _load_cached(acc, "metadata")

    print(f"[{acc}] Downloading...")
    dest = DATA_DIR / acc
    dest.mkdir(exist_ok=True)
    gse = GEOparse.get_GEO(geo=acc, destdir=str(dest), silent=True)

    expr = gse.pivot_samples("VALUE")
    meta = gse.phenotype_data

    meta["condition"] = meta["title"].apply(
        lambda t: "IRI" if "AKI" in str(t).upper() else "control"
    )

    print(f"  Expression: {expr.shape}")
    print(f"  Conditions: {meta['condition'].value_counts().to_dict()}")

    _save_cache(expr, acc, "expression")
    _save_cache(meta, acc, "metadata")
    return expr, meta

expr_53769, meta_53769 = download_GSE53769()


# %% [markdown]
# ## 2. Analisi di Espressione Differenziale
#
# Per ogni dataset, confrontiamo IRI vs control.
# - **Microarray**: t-test con correzione FDR (Benjamini-Hochberg)
# - **RNA-seq**: pyDESeq2


# %% 2.1 — DE per Microarray
def de_analysis_microarray(expr, meta, group_col="condition",
                            case="IRI", control="control",
                            dataset_name=""):
    print(f"\n[DE - {dataset_name}] Microarray: {case} vs {control}")

    case_samples = [s for s in meta[meta[group_col] == case].index if s in expr.columns]
    ctrl_samples = [s for s in meta[meta[group_col] == control].index if s in expr.columns]

    print(f"  Case ({case}): {len(case_samples)} samples")
    print(f"  Control ({control}): {len(ctrl_samples)} samples")

    if len(case_samples) < 2 or len(ctrl_samples) < 2:
        print("  [SKIP] Non abbastanza campioni per DE analysis")
        return pd.DataFrame()

    case_data = expr[case_samples]
    ctrl_data = expr[ctrl_samples]

    if expr.max().max() > 100:
        print("  Applying log2 transform...")
        case_data = np.log2(case_data + 1)
        ctrl_data = np.log2(ctrl_data + 1)

    results = []
    for gene in case_data.index:
        case_vals = case_data.loc[gene].dropna().values.astype(float)
        ctrl_vals = ctrl_data.loc[gene].dropna().values.astype(float)
        if len(case_vals) < 2 or len(ctrl_vals) < 2:
            continue
        tstat, pval = stats.ttest_ind(case_vals, ctrl_vals, equal_var=False)
        log2fc = case_vals.mean() - ctrl_vals.mean()
        results.append({
            "gene": gene,
            "log2FoldChange": log2fc,
            "pvalue": pval,
            "t_statistic": tstat,
            "mean_case": case_vals.mean(),
            "mean_ctrl": ctrl_vals.mean(),
        })

    df = pd.DataFrame(results).set_index("gene")

    from statsmodels.stats.multitest import multipletests
    _, df["padj"], _, _ = multipletests(df["pvalue"].fillna(1), method="fdr_bh")

    sig = df[(df["padj"] < 0.05) & (df["log2FoldChange"].abs() > 1)]
    print(f"  Total genes tested: {len(df)}")
    print(f"  Significant (|log2FC|>1, FDR<0.05): {len(sig)}")
    print(f"    Up: {(sig['log2FoldChange'] > 0).sum()}, Down: {(sig['log2FoldChange'] < 0).sum()}")

    df["dataset"] = dataset_name
    return df


# %% 2.2 — DE per RNA-seq (pyDESeq2)
def de_analysis_rnaseq(counts, meta, group_col="condition",
                        case="IRI", control="control",
                        dataset_name=""):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    print(f"\n[DE - {dataset_name}] RNA-seq: {case} vs {control}")

    case_samples = [s for s in meta[meta[group_col] == case].index if s in counts.columns]
    ctrl_samples = [s for s in meta[meta[group_col] == control].index if s in counts.columns]

    print(f"  Case ({case}): {len(case_samples)} samples")
    print(f"  Control ({control}): {len(ctrl_samples)} samples")

    if len(case_samples) < 2 or len(ctrl_samples) < 2:
        print("  [SKIP] Non abbastanza campioni")
        return pd.DataFrame()

    all_samples = ctrl_samples + case_samples
    conditions = ["control"] * len(ctrl_samples) + ["IRI"] * len(case_samples)

    count_matrix = counts[all_samples].T.astype(int)
    count_matrix = count_matrix.loc[:, count_matrix.sum() >= 10]

    clinical = pd.DataFrame({"condition": conditions}, index=all_samples)

    dds = DeseqDataSet(counts=count_matrix, metadata=clinical, design="~condition")
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=["condition", "IRI", "control"])
    stat_res.summary()

    results = stat_res.results_df.copy()
    results["dataset"] = dataset_name

    sig = results[(results["padj"] < 0.05) & (results["log2FoldChange"].abs() > 1)]
    print(f"  Total genes tested: {len(results)}")
    print(f"  Significant (|log2FC|>1, FDR<0.05): {len(sig)}")
    return results


# %% 2.3 — Esegui DE per tutti i dataset
all_de_results = {}

print("=" * 60)
print("ANALISI ESPRESSIONE DIFFERENZIALE")
print("=" * 60)

all_de_results["GSE43974"] = de_analysis_microarray(
    expr_43974, meta_43974, dataset_name="GSE43974"
)

all_de_results["GSE54888"] = de_analysis_microarray(
    expr_54888, meta_54888, dataset_name="GSE54888"
)

# GSE37838 — ora ha condition column
if "condition" in meta_37838.columns and "IRI" in meta_37838["condition"].values:
    all_de_results["GSE37838"] = de_analysis_microarray(
        expr_37838, meta_37838, dataset_name="GSE37838"
    )
else:
    print("\n[GSE37838] SKIP — no IRI samples found")

all_de_results["GSE53769"] = de_analysis_microarray(
    expr_53769, meta_53769, dataset_name="GSE53769"
)

# RNA-seq (se counts disponibili)
if _is_cached("GSE126805", "counts"):
    counts_126805 = _load_cached("GSE126805", "counts")
    all_de_results["GSE126805"] = de_analysis_rnaseq(
        counts_126805, meta_126805, dataset_name="GSE126805"
    )

# Salva
for name, df in all_de_results.items():
    if not df.empty:
        df.to_csv(OUTPUT_DIR / f"de_results_{name}.csv")
        print(f"  Saved: de_results_{name}.csv")

print(f"\nDataset con risultati DE: {[k for k,v in all_de_results.items() if not v.empty]}")

# %%
from IPython.display import display
print("expr_37838")
display(expr_37838)
print("meta_37838")
display(meta_37838)
print("expr_43974")
display(expr_43974)
print("meta_43974")
display(meta_43974)
print("expr_54888")
display(expr_54888)
print("meta_54888")
display(meta_54888)
print("expr_53769")
display(expr_53769)
print("meta_53769")
display(meta_53769)
print("meta_126805")
display(meta_126805)
print("meta_142077")
display(meta_142077)
print("meta_293480")
display(meta_293480)
print("meta_90861")
display(meta_90861)
print("counts_126805")
display(counts_126805)
print("all_de_results")
display(all_de_results)
print("df")
display(df)


# %% [markdown]
# ## 3. Meta-analisi Cross-Dataset

# %% 3.1 — Mappatura probe IDs → gene symbols (FIXED)
def _get_gpl_annotation(accession: str, gpl_id: str) -> dict:
    """
    Scarica/carica la platform annotation e restituisce un dict probe_id -> gene_symbol.
    Cerca prima i file GPL già scaricati, poi prova via GEOparse.
    """
    dest = DATA_DIR / accession

    # 1. Cerca file GPL già presenti
    gpl_files = sorted(
        list(dest.glob("GPL*.txt")) +
        list(dest.glob("GPL*.gz")) +
        list(dest.glob("GPL*.annot*")),
        key=lambda f: f.stat().st_size,
        reverse=True  # Il più grande per primo (più probabilmente l'annotation completa)
    )

    for gpl_f in gpl_files:
        try:
            comp = "gzip" if str(gpl_f).endswith(".gz") else None
            ann = pd.read_csv(gpl_f, sep="\t", comment="#", low_memory=False, compression=comp)

            symbol_cols = [c for c in ann.columns if any(
                kw in c.lower() for kw in ["symbol", "gene_symbol", "ilmn_gene"]
            ) and "synonym" not in c.lower() and "source" not in c.lower()]

            id_cols = [c for c in ann.columns if c in ("ID", "ID_REF", "PROBE_ID")]

            if symbol_cols and id_cols:
                id_col, sym_col = id_cols[0], symbol_cols[0]
                clean = ann[[id_col, sym_col]].dropna()
                clean = clean[clean[sym_col].astype(str).str.strip().ne("")]
                clean = clean[~clean[sym_col].astype(str).str.strip().isin(["nan", "---", "NA"])]
                mapping = dict(zip(clean[id_col].astype(str), clean[sym_col].astype(str).str.strip()))
                if len(mapping) > 100:
                    print(f"    File: {gpl_f.name}, col: {sym_col} -> {len(mapping)} mappings")
                    return mapping
        except Exception:
            continue

    # 2. Scarica via GEOparse
    print(f"    Downloading GPL {gpl_id} via GEOparse...")
    try:
        gpl = GEOparse.get_GEO(geo=gpl_id, destdir=str(dest), silent=True)
        ann = gpl.table

        symbol_cols = [c for c in ann.columns if any(
            kw in c.lower() for kw in ["symbol", "gene_symbol", "ilmn_gene"]
        ) and "synonym" not in c.lower() and "source" not in c.lower()]

        if symbol_cols:
            sym_col = symbol_cols[0]
            clean = ann[["ID", sym_col]].dropna()
            clean = clean[~clean[sym_col].astype(str).str.strip().isin(["", "nan", "---", "NA"])]
            mapping = dict(zip(clean["ID"].astype(str), clean[sym_col].astype(str).str.strip()))
            print(f"    GPL table col: {sym_col} -> {len(mapping)} mappings")
            return mapping
    except Exception as e:
        print(f"    [WARN] GPL download failed: {e}")

    return {}


# Mapping: dataset -> GPL ID
GPL_MAP = {
    "GSE43974": "GPL10558",   # Illumina HumanHT-12 V4.0
    "GSE37838": "GPL570",     # Affymetrix U133 Plus 2.0
    "GSE54888": "GPL6244",    # Affymetrix Gene 1.0 ST
    "GSE53769": "GPL16686",   # Affymetrix Gene 2.0 ST
}


def map_probes_to_genes(de_results: dict) -> dict:
    """Mappa tutti i probe IDs a gene symbols."""
    mapped = {}

    for name, df in de_results.items():
        if df.empty:
            continue

        sample_idx = str(df.index[0])

        # Tutti i microarray hanno probe IDs (ILMN_, numeri, *_at, *_st)
        needs_mapping = (
            sample_idx.startswith("ILMN_") or
            sample_idx.replace("_", "").replace(".", "").isnumeric() or
            sample_idx.endswith("_at") or
            sample_idx.endswith("_st") or
            sample_idx.endswith("_s_at") or
            sample_idx.endswith("_x_at")
        )

        if needs_mapping:
            print(f"\n  [{name}] Mapping probes (example: {sample_idx})...")

            gpl_id = GPL_MAP.get(name, "")
            if not gpl_id:
                print(f"    [SKIP] No GPL mapping configured for {name}")
                continue

            mapping = _get_gpl_annotation(name, gpl_id)

            if not mapping:
                print(f"    [SKIP] No annotation found for {name}")
                continue

            df = df.copy()
            df["gene_symbol"] = df.index.astype(str).map(mapping)

            n_before = len(df)
            df = df.dropna(subset=["gene_symbol"])
            df = df[~df["gene_symbol"].isin(["", "nan", "---", "NA"])]
            n_mapped = len(df)

            if n_mapped == 0:
                print(f"    [SKIP] No probes mapped to gene symbols")
                continue

            # Multiple probes per gene: tieni quella con p-value più basso
            df = df.sort_values("pvalue").groupby("gene_symbol").first()
            print(f"    Mapped {n_mapped}/{n_before} probes -> {len(df)} unique genes")

        else:
            # RNA-seq: già gene symbols
            df = df.copy()
            df["gene_symbol"] = df.index
            print(f"\n  [{name}] Already gene symbols ({len(df)} genes)")

        mapped[name] = df

    print(f"\n  Datasets mappati: {list(mapped.keys())}")
    return mapped


# %% 3.2 — Meta-analisi con metodo di Fisher (FIXED)
def meta_analysis_fisher(de_results: dict, min_datasets: int = 2) -> pd.DataFrame:
    print("\n" + "=" * 60)
    print(f"META-ANALISI (Fisher, min {min_datasets} datasets)")
    print("=" * 60)

    # Diagnostica
    for name, df in de_results.items():
        if not df.empty:
            examples = df.index[:3].tolist()
            print(f"  {name}: {len(df)} genes, examples: {examples}")

    # Raccogli p-value e log2FC per ogni gene symbol
    gene_data = defaultdict(lambda: {"pvalues": [], "log2fcs": [], "datasets": []})

    for name, df in de_results.items():
        if df.empty:
            continue
        for gene in df.index:
            if "gene_symbol" in df.columns:
                symbol = str(df.loc[gene, "gene_symbol"])
            else:
                symbol = str(gene)

            if symbol in ("", "nan", "None", "---", "NaN"):
                continue

            pval = df.loc[gene, "pvalue"] if "pvalue" in df.columns else None
            l2fc = df.loc[gene, "log2FoldChange"] if "log2FoldChange" in df.columns else 0

            if pd.notna(pval) and 0 < pval < 1:
                gene_data[symbol]["pvalues"].append(pval)
                gene_data[symbol]["log2fcs"].append(l2fc)
                gene_data[symbol]["datasets"].append(name)

    # Conta geni con dati in >= min_datasets
    multi_ds = {g: d for g, d in gene_data.items() if len(d["pvalues"]) >= min_datasets}
    print(f"\n  Geni totali raccolti: {len(gene_data)}")
    print(f"  Geni in >= {min_datasets} datasets: {len(multi_ds)}")

    # Fallback se nessun gene in comune
    if len(multi_ds) == 0:
        print("\n  [WARN] Nessun gene in >= 2 dataset!")
        print("  Fallback: uso solo GSE43974...")
        multi_ds = {g: d for g, d in gene_data.items() if len(d["pvalues"]) >= 1}
        if len(multi_ds) == 0:
            print("  [ERROR] Nessun gene disponibile!")
            return pd.DataFrame()

    # Combina
    meta_results = []
    for gene, data in multi_ds.items():
        if len(data["pvalues"]) >= 2:
            _, combined_p = stats.combine_pvalues(data["pvalues"], method="fisher")
        else:
            combined_p = data["pvalues"][0]

        weights = [-np.log10(p + 1e-300) for p in data["pvalues"]]
        weighted_l2fc = np.average(data["log2fcs"], weights=weights)

        directions = [1 if l > 0 else -1 for l in data["log2fcs"]]
        direction_consistency = abs(sum(directions)) / len(directions)

        meta_results.append({
            "gene": gene,
            "meta_pvalue": combined_p,
            "meta_log2FC": weighted_l2fc,
            "n_datasets": len(data["pvalues"]),
            "direction_consistency": direction_consistency,
            "datasets": ",".join(data["datasets"]),
        })

    meta_df = pd.DataFrame(meta_results).set_index("gene")

    from statsmodels.stats.multitest import multipletests
    _, meta_df["meta_padj"], _, _ = multipletests(
        meta_df["meta_pvalue"].fillna(1), method="fdr_bh"
    )
    meta_df = meta_df.sort_values("meta_pvalue")

    sig = meta_df[
        (meta_df["meta_padj"] < 0.05) &
        (meta_df["meta_log2FC"].abs() > 0.5) &
        (meta_df["direction_consistency"] > 0.5)
    ]
    print(f"\nGeni totali nella meta-analisi: {len(meta_df)}")
    print(f"Significativi (FDR<0.05, |log2FC|>0.5, consistent): {len(sig)}")
    print(f"  Up in IRI: {(sig['meta_log2FC'] > 0).sum()}")
    print(f"  Down in IRI: {(sig['meta_log2FC'] < 0).sum()}")
    print(f"\nTop 20 geni:")
    print(sig.head(20)[["meta_log2FC", "meta_padj", "n_datasets", "direction_consistency"]])

    return meta_df


# %% 3.3 — Esegui mapping + meta-analisi
mapped_results = map_probes_to_genes(all_de_results)
meta_df = meta_analysis_fisher(mapped_results, min_datasets=2)

if len(meta_df) > 0:
    meta_df.to_csv(OUTPUT_DIR / "meta_analysis_results.csv")

    seed_genes = meta_df[
        (meta_df["meta_padj"] < 0.05) &
        (meta_df["meta_log2FC"].abs() > 0.5) &
        (meta_df["direction_consistency"] > 0.5)
    ].copy()
    seed_genes.to_csv(OUTPUT_DIR / "seed_genes.csv")
    print(f"\n==> {len(seed_genes)} SEED GENES identificati e salvati")
else:
    print("\n[WARN] Meta-analisi vuota. Controlla probe mapping.")
    seed_genes = pd.DataFrame()


# %% [markdown]
# ## 4. Pathway Enrichment Analysis

# %% 4.1 — GSEA (pre-ranked)
import gseapy as gp

def run_gsea_preranked(meta_df, gene_sets=None):
    if gene_sets is None:
        gene_sets = [
            "KEGG_2021_Human",
            "Reactome_2022",
            "GO_Biological_Process_2023",
            "MSigDB_Hallmark_2020",
        ]

    if len(meta_df) == 0:
        print("[GSEA] Skipping — no meta-analysis results")
        return pd.DataFrame()

    ranking = meta_df["meta_log2FC"].dropna().sort_values(ascending=False)

    all_results = []
    for gs_name in gene_sets:
        try:
            print(f"\n[GSEA] Running {gs_name}...")
            res = gp.prerank(
                rnk=ranking,
                gene_sets=gs_name,
                outdir=str(OUTPUT_DIR / "gsea" / gs_name),
                seed=42,
                permutation_num=1000,
                no_plot=True,
                min_size=10,
                max_size=500,
            )
            res_df = res.res2d.copy()
            res_df["library"] = gs_name
            all_results.append(res_df)
            sig = res_df[res_df["FDR q-val"] < 0.05]
            print(f"  -> {len(sig)} significant terms (FDR < 0.05)")
        except Exception as e:
            print(f"  [WARN] {gs_name} failed: {e}")

    if all_results:
        return pd.concat(all_results)
    return pd.DataFrame()

gsea_results = run_gsea_preranked(meta_df)
if not gsea_results.empty:
    gsea_results.to_csv(OUTPUT_DIR / "gsea_results.csv")
    print(f"\nGSEA results saved ({len(gsea_results)} terms)")


# %% 4.2 — ORA (Over-Representation Analysis)
def run_ora(seed_genes, gene_sets=None):
    if gene_sets is None:
        gene_sets = ["KEGG_2021_Human", "Reactome_2022",
                      "GO_Biological_Process_2023"]

    if len(seed_genes) == 0:
        print("[ORA] Skipping — no seed genes")
        return pd.DataFrame()

    up_genes = seed_genes[seed_genes["meta_log2FC"] > 0].index.tolist()
    down_genes = seed_genes[seed_genes["meta_log2FC"] < 0].index.tolist()

    results = []
    for direction, genes in [("up", up_genes), ("down", down_genes)]:
        if len(genes) < 5:
            print(f"[ORA] {direction}: only {len(genes)} genes, skipping (need >= 5)")
            continue
        try:
            print(f"\n[ORA] {direction}-regulated: {len(genes)} genes")
            res = gp.enrichr(
                gene_list=genes,
                gene_sets=gene_sets,
                outdir=str(OUTPUT_DIR / "ora" / direction),
                no_plot=True,
            )
            res_df = res.results.copy()
            res_df["direction"] = direction
            results.append(res_df)
        except Exception as e:
            print(f"  [WARN] ORA failed: {e}")

    if results:
        return pd.concat(results)
    return pd.DataFrame()

ora_results = run_ora(seed_genes)
if not ora_results.empty:
    ora_results.to_csv(OUTPUT_DIR / "ora_results.csv")
    print(f"\nORA results saved ({len(ora_results)} terms)")


# %% [markdown]
# ## 5. Riepilogo Fase 1

# %% Riepilogo
print("\n" + "=" * 60)
print("RIEPILOGO FASE 1")
print("=" * 60)

print(f"\nDataset processati: {len(all_de_results)}")
for name, df in all_de_results.items():
    if not df.empty:
        sig_count = len(df[df.get("padj", df.get("pvalue", pd.Series(dtype=float))) < 0.05])
        print(f"  {name}: {len(df)} genes tested, {sig_count} significant (p<0.05)")

if len(meta_df) > 0:
    print(f"\nMeta-analisi: {len(meta_df)} genes")
    print(f"Seed genes: {len(seed_genes)}")

if not gsea_results.empty:
    sig_pathways = gsea_results[gsea_results["FDR q-val"] < 0.05]
    print(f"GSEA significant pathways: {len(sig_pathways)}")

if not ora_results.empty:
    sig_ora = ora_results[ora_results["Adjusted P-value"] < 0.05]
    print(f"ORA significant terms: {len(sig_ora)}")

print(f"\nFile salvati in: {OUTPUT_DIR}")
print("Pronti per Fase 2 (Network Medicine)")