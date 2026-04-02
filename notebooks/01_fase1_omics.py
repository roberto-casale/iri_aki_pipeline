# %% [markdown]
# # Fase 1 — Meta-analisi Omici IR-AKI
#
# ## Obiettivo
# Identificare geni e pathway consistentemente alterati nel danno renale
# da ischemia-riperfusione (IR-AKI), attraverso una meta-analisi multi-dataset.
#
# ## Strategia
# Per ogni dataset → un confronto principale **IRI vs Control**, dove i gruppi
# sono definiti in modo biologicamente coerente. I risultati vengono poi
# combinati cross-dataset con il metodo di Fisher.
#
# ---
#
# ## GSE43974 — Damman et al., Transplantation 2015
#
# **Piattaforma:** Illumina HumanHT-12 V4 (GPL10558), 47323 sonde
#
# **Design dello studio:**
# 554 biopsie renali da donatori viventi (Living), dopo morte cerebrale (DBD)
# e dopo arresto cardiaco (DCD), a 3 timepoint:
# - **T1** = prima del retrieval (Living n≈37, DBD n≈82, DCD n≈38)
# - **T2** = fine ischemia fredda (DBD n≈110, DCD n≈53)
# - **T3** = 45-60 min dopo riperfusione (DBD n≈105, DCD n≈64, Living n≈34)
#
# **Selezione campioni per la nostra meta-analisi:**
#
# | Confronto | Case (IRI) | Control | Razionale |
# |-----------|-----------|---------|-----------|
# | **Principale** | DBD_T3 + DCD_T3 (riperfusione) | Living_T1 (baseline sano) | Cattura l'intero danno IR: ischemia + riperfusione vs rene sano |
# | **Supporto** | DBD_T1 + DCD_T1 (post-ischemia, pre-retrieval) | Living_T1 (baseline sano) | Cattura la componente ischemica sola (brain death / cardiac arrest) |
#
# Il confronto principale (T3 vs Living_T1) è quello più comparabile con gli
# altri dataset (GSE54888, GSE37838, GSE53769, ecc.) che confrontano tutti
# biopsie post-riperfusione con controlli sani o pre-danno.
#
# **Dati di espressione:**
# I dati da GEO sono già processati con GeneSpring: log2 + median-shift al
# 75° percentile + baseline transformation (mediana di tutti i campioni
# sottratta). I valori sono centrati attorno a 0.
# - FC ≥ 1.3 nell'articolo → |log2FC| ≥ log2(1.3) ≈ 0.379
# - Noi usiamo |log2FC| > 0.379, FDR < 0.05

# %% Setup
import os, sys, warnings
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
from collections import defaultdict
from IPython.display import display

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(".").resolve().parent
DATA_DIR = PROJECT_ROOT / "data" / "geo"
OUTPUT_DIR = PROJECT_ROOT / "output" / "fase1"
DATA_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"Project root: {PROJECT_ROOT}")

# Soglie globali per la DE
FC_THRESHOLD = 0.379      # |log2FC| ≥ log2(1.3)
FDR_THRESHOLD = 0.05


# %% [markdown]
# ## 1. Caricamento GSE43974

import GEOparse

def load_GSE43974():
    acc = "GSE43974"
    dest = DATA_DIR / acc
    dest.mkdir(exist_ok=True)
    expr_pkl, meta_pkl = dest / "expression.pkl", dest / "metadata.pkl"

    if expr_pkl.exists() and meta_pkl.exists():
        print(f"[{acc}] Loading from cache...")
        expr = pd.read_pickle(expr_pkl)
        meta = pd.read_pickle(meta_pkl)
    else:
        print(f"[{acc}] Downloading...")
        gse = GEOparse.get_GEO(geo=acc, destdir=str(dest), silent=True)
        expr = gse.pivot_samples("VALUE")
        meta = gse.phenotype_data
        expr.to_pickle(expr_pkl)
        meta.to_pickle(meta_pkl)

    # Parse metadata
    meta["donor_type"] = meta["title"].str.extract(r"Kidney_donor_(\w+)_T")
    meta["timepoint"] = meta["title"].str.extract(r"_T(\d)_").astype(float)

    # DGF column (per analisi future)
    # Troviamo tutte le colonne che parlano di DGF
    dgf_cols = [c for c in meta.columns if "delayed graft function" in c.lower()]
    if dgf_cols:
        # 1. Fonde le colonne trovate riempiendo i buchi a vicenda
        meta["dgf"] = meta[dgf_cols].bfill(axis=1).iloc[:, 0]
        # 2. Converte in numeri. Testi strani o celle vuote diventano automaticamente NaN grazie a errors="coerce"
        meta["dgf"] = pd.to_numeric(meta["dgf"], errors="coerce")
        # (Opzionale ma consigliato) Assicuriamoci che ci siano SOLO 0 e 1. 
        # Qualsiasi altro numero anomalo viene trasformato in NaN.
        meta.loc[~meta["dgf"].isin([0.0, 1.0]), "dgf"] = np.nan
        # 3. Calcola i totali per il resoconto
        totale_pazienti = len(meta)
        dati_validi = meta["dgf"].notna().sum()
        dati_mancanti = totale_pazienti - dati_validi
        # 4. Se mancano dei dati, lancia un warning visibile a schermo
        if dati_mancanti > 0:
            print(f"⚠️ WARNING: Trovati {dati_mancanti} pazienti senza un dato valido (0 o 1) per la DGF.")
            print("   Le celle corrispondenti sono state impostate su 'NaN' (Not a Number).")
        # 5. Stampa il riassunto finale richiesto
        print(f"✅ Inserimento DGF completato: inseriti {dati_validi} dati validi su un totale di {totale_pazienti} pazienti.")
    else:
        print("⚠️ Nessuna colonna relativa alla DGF trovata nei metadati.")

    print(f"  Expression: {expr.shape[0]} probes × {expr.shape[1]} samples")
    print(f"  Gruppi:")
    for (dt, tp), n in meta.groupby(["donor_type", "timepoint"]).size().items():
        print(f"    {dt}_T{int(tp)}: {n}")

    return expr, meta

expr_43974, meta_43974 = load_GSE43974()

# Verifica scala dati
print(f"\n  Scala dati: min={expr_43974.values.min():.2f}, "
      f"max={expr_43974.values.max():.2f}, mean={expr_43974.values.mean():.3f}")
print(f"  → Dati baseline-transformed (centrati ~0), NO log2 aggiuntivo")


# %% [markdown]
# ## 2. Probe → Gene Symbol Mapping

# %%
def get_probe_gene_mapping(accession, gpl_id, data_dir):
    """Carica mapping probe_id → gene_symbol dalla GPL annotation."""
    dest = data_dir / accession

    gpl_files = sorted(
        list(dest.glob("GPL*.txt")) + list(dest.glob("GPL*.gz")),
        key=lambda f: f.stat().st_size, reverse=True
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
                clean = clean[~clean[sym_col].astype(str).str.strip().isin(["", "nan", "---", "NA"])]
                mapping = dict(zip(clean[id_col].astype(str), clean[sym_col].astype(str).str.strip()))
                if len(mapping) > 100:
                    print(f"  [{accession}] GPL file: {gpl_f.name} → {len(mapping)} mappings")
                    return mapping
        except Exception:
            continue

    # Fallback: download
    print(f"  [{accession}] Downloading GPL {gpl_id}...")
    try:
        gpl = GEOparse.get_GEO(geo=gpl_id, destdir=str(dest), silent=True)
        ann = gpl.table
        symbol_cols = [c for c in ann.columns if "symbol" in c.lower()
                       and "synonym" not in c.lower()]
        if symbol_cols:
            clean = ann[["ID", symbol_cols[0]]].dropna()
            clean = clean[~clean[symbol_cols[0]].astype(str).str.strip().isin(["", "nan", "---"])]
            mapping = dict(zip(clean["ID"].astype(str), clean[symbol_cols[0]].astype(str).str.strip()))
            print(f"  [{accession}] GPL table → {len(mapping)} mappings")
            return mapping
    except Exception as e:
        print(f"  [ERROR] {e}")
    return {}


probe_to_gene_43974 = get_probe_gene_mapping("GSE43974", "GPL10558", DATA_DIR)
mapped = set(expr_43974.index.astype(str)) & set(probe_to_gene_43974.keys())
print(f"  Sonde mappabili: {len(mapped)}/{len(expr_43974)} ({100*len(mapped)/len(expr_43974):.1f}%)")


# %% [markdown]
# ## 3. Espressione Differenziale

# %%
def de_microarray(expr_df, samples_case, samples_ctrl, comparison_name=""):
    """
    DE con Welch t-test per microarray. Nessun log2 aggiuntivo
    (i dati sono già baseline-transformed).
    """
    case = [s for s in samples_case if s in expr_df.columns]
    ctrl = [s for s in samples_ctrl if s in expr_df.columns]

    print(f"\n[DE] {comparison_name}")
    print(f"  Case: {len(case)}, Control: {len(ctrl)}")

    if len(case) < 2 or len(ctrl) < 2:
        print("  [SKIP] <2 campioni per gruppo")
        return pd.DataFrame()

    case_data = expr_df[case]
    ctrl_data = expr_df[ctrl]

    results = []
    for probe in expr_df.index:
        cv = case_data.loc[probe].dropna().values.astype(float)
        tv = ctrl_data.loc[probe].dropna().values.astype(float)
        if len(cv) < 2 or len(tv) < 2:
            continue
        tstat, pval = stats.ttest_ind(cv, tv, equal_var=False)
        results.append({
            "probe": probe,
            "log2FoldChange": cv.mean() - tv.mean(),
            "pvalue": pval,
            "t_statistic": tstat,
            "mean_case": cv.mean(),
            "mean_ctrl": tv.mean(),
        })

    df = pd.DataFrame(results).set_index("probe")

    from statsmodels.stats.multitest import multipletests
    _, df["padj"], _, _ = multipletests(df["pvalue"].fillna(1), method="fdr_bh")

    sig = df[(df["padj"] < FDR_THRESHOLD) & (df["log2FoldChange"].abs() > FC_THRESHOLD)]
    print(f"  Tested: {len(df)} probes")
    print(f"  Significant (|log2FC|>{FC_THRESHOLD}, FDR<{FDR_THRESHOLD}): "
          f"{len(sig)} ({(sig['log2FoldChange']>0).sum()}↑ {(sig['log2FoldChange']<0).sum()}↓)")

    df["comparison"] = comparison_name
    df["dataset"] = "GSE43974"
    return df


# %% Definisci gruppi
m = meta_43974
living_t1 = m[(m["donor_type"] == "Living") & (m["timepoint"] == 1)].index.tolist()
dbd_t1    = m[(m["donor_type"] == "BD")     & (m["timepoint"] == 1)].index.tolist()
dcd_t1    = m[(m["donor_type"] == "DCD")    & (m["timepoint"] == 1)].index.tolist()
dbd_t3    = m[(m["donor_type"] == "BD")     & (m["timepoint"] == 3)].index.tolist()
dcd_t3    = m[(m["donor_type"] == "DCD")    & (m["timepoint"] == 3)].index.tolist()

deceased_t1 = dbd_t1 + dcd_t1    # ischemia (pre-retrieval)
deceased_t3 = dbd_t3 + dcd_t3    # riperfusione

print(f"Living_T1 (control):     {len(living_t1)}")
print(f"Deceased_T1 (ischemia):  {len(deceased_t1)} (DBD:{len(dbd_t1)}, DCD:{len(dcd_t1)})")
print(f"Deceased_T3 (IRI):       {len(deceased_t3)} (DBD:{len(dbd_t3)}, DCD:{len(dcd_t3)})")


# %% Confronto principale: Deceased_T3 vs Living_T1
de_main = de_microarray(
    expr_43974, deceased_t3, living_t1,
    comparison_name="IRI_vs_Control (T3_deceased vs T1_living)"
)

# Confronto supporto: Deceased_T1 vs Living_T1
de_support = de_microarray(
    expr_43974, deceased_t1, living_t1,
    comparison_name="Ischemia_vs_Control (T1_deceased vs T1_living)"
)


# %% [markdown]
# ## 4. Probe → Gene + Aggregazione per gene

# %%
def probes_to_genes(de_df, probe_gene_map):
    """Converte DE da probe-level a gene-level (best probe per gene)."""
    if de_df.empty:
        return de_df
    df = de_df.copy()
    df["gene_symbol"] = df.index.astype(str).map(probe_gene_map)
    n0 = len(df)
    df = df.dropna(subset=["gene_symbol"])
    df = df[~df["gene_symbol"].isin(["", "nan", "---", "NA"])]
    # Multi-probe per gene: tieni quella con p-value più basso
    df = df.sort_values("pvalue").groupby("gene_symbol").first()
    print(f"  Mapped: {len(df)}/{n0} probes → {len(df)} unique genes")
    return df


de_main_genes = probes_to_genes(de_main, probe_to_gene_43974)
de_support_genes = probes_to_genes(de_support, probe_to_gene_43974)


# %% [markdown]
# ## 5. Combinazione confronti intra-dataset (Fisher)
#
# Combiniamo il confronto principale e quello di supporto per ottenere
# un ranking di geni robusto per GSE43974.
# Questo ranking sarà poi combinato con quelli degli altri dataset.

# %%
def combine_de_results_fisher(de_list, names, min_comparisons=1):
    """
    Combina più DE results (gene-level) con Fisher's method.
    Restituisce un DataFrame con p-value combinato e log2FC pesato.
    """
    gene_data = defaultdict(lambda: {"pvalues": [], "log2fcs": [], "sources": []})

    for de_df, name in zip(de_list, names):
        if de_df.empty:
            continue
        for gene in de_df.index:
            pval = de_df.loc[gene, "pvalue"]
            l2fc = de_df.loc[gene, "log2FoldChange"]
            if pd.notna(pval) and 0 < pval < 1:
                gene_data[gene]["pvalues"].append(pval)
                gene_data[gene]["log2fcs"].append(l2fc)
                gene_data[gene]["sources"].append(name)

    results = []
    for gene, d in gene_data.items():
        n = len(d["pvalues"])
        if n < min_comparisons:
            continue

        if n >= 2:
            _, combined_p = stats.combine_pvalues(d["pvalues"], method="fisher")
        else:
            combined_p = d["pvalues"][0]

        weights = [-np.log10(p + 1e-300) for p in d["pvalues"]]
        weighted_l2fc = np.average(d["log2fcs"], weights=weights)
        directions = [1 if l > 0 else -1 for l in d["log2fcs"]]
        consistency = abs(sum(directions)) / len(directions)

        results.append({
            "gene": gene,
            "meta_pvalue": combined_p,
            "meta_log2FC": weighted_l2fc,
            "n_comparisons": n,
            "direction_consistency": consistency,
            "sources": "|".join(d["sources"]),
        })

    meta_df = pd.DataFrame(results).set_index("gene")

    from statsmodels.stats.multitest import multipletests
    _, meta_df["meta_padj"], _, _ = multipletests(
        meta_df["meta_pvalue"].fillna(1), method="fdr_bh"
    )
    meta_df = meta_df.sort_values("meta_pvalue")
    return meta_df


combined_43974 = combine_de_results_fisher(
    [de_main_genes, de_support_genes],
    ["IRI_vs_Control", "Ischemia_vs_Control"],
    min_comparisons=1
)

# Seed genes da GSE43974
seed_43974 = combined_43974[
    (combined_43974["meta_padj"] < FDR_THRESHOLD) &
    (combined_43974["meta_log2FC"].abs() > FC_THRESHOLD) &
    (combined_43974["direction_consistency"] >= 0.5)
].copy()

print(f"\nGSE43974 — Risultati combinati:")
print(f"  Geni totali: {len(combined_43974)}")
print(f"  Seed genes: {len(seed_43974)}")
print(f"    Up in IRI:   {(seed_43974['meta_log2FC'] > 0).sum()}")
print(f"    Down in IRI: {(seed_43974['meta_log2FC'] < 0).sum()}")
print(f"    In 2/2 confronti: {(seed_43974['n_comparisons'] == 2).sum()}")

print(f"\nTop 20 seed genes:")
print(seed_43974.head(20)[["meta_log2FC", "meta_padj", "n_comparisons",
                            "direction_consistency"]].to_string())


# %% [markdown]
# ## 6. Salvataggio risultati intermedi
#
# Salviamo tutto ciò che serve per la meta-analisi cross-dataset.
# Quando aggiungiamo un nuovo dataset, produrremo lo stesso formato
# (un DataFrame gene-level con log2FC, pvalue, padj) e lo combineremo.

# %%
# Salva DE gene-level (serviranno per la meta-analisi cross-dataset)
de_main_genes.to_csv(OUTPUT_DIR / "de_GSE43974_IRI_vs_Control.csv")
de_support_genes.to_csv(OUTPUT_DIR / "de_GSE43974_Ischemia_vs_Control.csv")

# Salva combinato intra-dataset
combined_43974.to_csv(OUTPUT_DIR / "combined_GSE43974.csv")
seed_43974.to_csv(OUTPUT_DIR / "seed_genes_GSE43974.csv")

print(f"Salvati in {OUTPUT_DIR}/")


# %% [markdown]
# ## 7. Pathway Enrichment (GSEA + ORA)

# %%
import gseapy as gp

def run_gsea(ranking, output_subdir, gene_sets=None):
    """GSEA pre-ranked."""
    if gene_sets is None:
        gene_sets = ["KEGG_2021_Human", "Reactome_2022",
                     "GO_Biological_Process_2023", "MSigDB_Hallmark_2020"]

    ranking = ranking.dropna().sort_values(ascending=False)
    ranking = ranking[~ranking.index.duplicated(keep='first')]

    all_res = []
    for gs in gene_sets:
        try:
            print(f"  GSEA: {gs}...")
            res = gp.prerank(
                rnk=ranking, gene_sets=gs,
                outdir=str(OUTPUT_DIR / output_subdir / gs),
                seed=42, permutation_num=1000,
                no_plot=True, min_size=10, max_size=500,
            )
            df = res.res2d.copy()
            df["library"] = gs
            all_res.append(df)
            n_sig = (df["FDR q-val"] < 0.05).sum()
            print(f"    → {n_sig} terms FDR<0.05")
        except Exception as e:
            print(f"    [WARN] {e}")
    return pd.concat(all_res, ignore_index=True) if all_res else pd.DataFrame()


def run_ora(up_genes, down_genes, output_subdir, gene_sets=None):
    """ORA con Enrichr."""
    if gene_sets is None:
        gene_sets = ["KEGG_2021_Human", "Reactome_2022",
                     "GO_Biological_Process_2023"]
    results = []
    for direction, genes in [("up", up_genes), ("down", down_genes)]:
        if len(genes) < 5:
            print(f"  ORA {direction}: {len(genes)} geni (skip)")
            continue
        try:
            print(f"  ORA {direction}: {len(genes)} geni...")
            res = gp.enrichr(gene_list=genes, gene_sets=gene_sets,
                             outdir=str(OUTPUT_DIR / output_subdir / direction),
                             no_plot=True)
            df = res.results.copy()
            df["direction"] = direction
            results.append(df)
            n_sig = (df["Adjusted P-value"] < 0.05).sum()
            print(f"    → {n_sig} terms FDR<0.05")
        except Exception as e:
            print(f"    [WARN] {e}")
    return pd.concat(results, ignore_index=True) if results else pd.DataFrame()


# %% GSEA sul ranking combinato
print("=" * 60)
print("GSEA — GSE43974")
print("=" * 60)

gsea_results = run_gsea(combined_43974["meta_log2FC"], "gsea_GSE43974")

if not gsea_results.empty:
    gsea_results.to_csv(OUTPUT_DIR / "gsea_GSE43974.csv", index=False)
    sig_gsea = gsea_results[gsea_results["FDR q-val"] < 0.05].sort_values("FDR q-val")
    print(f"\nTop pathways (FDR<0.05):")
    for _, r in sig_gsea.head(15).iterrows():
        arrow = "↑" if float(r.get("NES", r.get("ES", 0))) > 0 else "↓"
        print(f"  {arrow} {r['Term']} (FDR={r['FDR q-val']:.4f}, {r['library']})")

# %% ORA sui seed genes
print("\n" + "=" * 60)
print("ORA — GSE43974")
print("=" * 60)

if len(seed_43974) > 0:
    up = seed_43974[seed_43974["meta_log2FC"] > 0].index.tolist()
    down = seed_43974[seed_43974["meta_log2FC"] < 0].index.tolist()
    ora_results = run_ora(up, down, "ora_GSE43974")
    if not ora_results.empty:
        ora_results.to_csv(OUTPUT_DIR / "ora_GSE43974.csv", index=False)


# %% [markdown]
# ## 8. Riepilogo
#
# Questo script produce per GSE43974:
# - `de_GSE43974_IRI_vs_Control.csv` — DE gene-level del confronto principale
# - `de_GSE43974_Ischemia_vs_Control.csv` — DE gene-level del confronto supporto
# - `combined_GSE43974.csv` — Risultato combinato Fisher (ranking per GSEA)
# - `seed_genes_GSE43974.csv` — Geni significativi (per ORA e fasi successive)
# - `gsea_GSE43974.csv` — Risultati GSEA
# - `ora_GSE43974.csv` — Risultati ORA
#
# **Per aggiungere un nuovo dataset** (es. GSE54888):
# 1. Carica dati + metadata
# 2. Definisci campioni IRI e Control
# 3. Esegui `de_microarray()` → `probes_to_genes()`
# 4. Il risultato gene-level si combina cross-dataset con `combine_de_results_fisher()`

# %%
print("\n" + "=" * 70)
print("RIEPILOGO FASE 1 — GSE43974")
print("=" * 70)
print(f"Campioni totali: {len(meta_43974)}")
print(f"  Control (Living_T1):     {len(living_t1)}")
print(f"  Case IRI (Deceased_T3):  {len(deceased_t3)}")
print(f"  Support (Deceased_T1):   {len(deceased_t1)}")
print(f"Probes: {len(expr_43974)}, Gene mappabili: {len(mapped)}")
print(f"Seed genes: {len(seed_43974)} "
      f"({(seed_43974['meta_log2FC']>0).sum()}↑, "
      f"{(seed_43974['meta_log2FC']<0).sum()}↓)")
print(f"\nFile in: {OUTPUT_DIR}/")
print("Pronto per aggiungere il prossimo dataset.")