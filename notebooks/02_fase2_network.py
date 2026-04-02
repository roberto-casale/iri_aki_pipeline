"""
MAGIC SOLUTION - Fase 2: Network Medicine e Target Prioritization
==================================================================
Costruisce una rete PPI, esegue network diffusion con seed genes
dalla Fase 1, e classifica i target per druggability.

Requisiti:
    pip install networkx numpy scipy pandas requests pymed

Input:  seed_genes.csv (dalla Fase 1)
Output: ranked_targets.csv, ppi_network.graphml
"""

import pandas as pd
import numpy as np
import networkx as nx
import requests
import time
from pathlib import Path

# ============================================================
# CONFIGURAZIONE
# ============================================================
INPUT_DIR = Path("../output/fase1")
OUTPUT_DIR = Path("../output/fase2")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# STRING parameters
STRING_SPECIES = 9606  # Homo sapiens
STRING_SCORE_THRESHOLD = 700  # High confidence

# RWR parameters
RWR_RESTART_PROB = 0.7
RWR_MAX_ITER = 100
RWR_TOL = 1e-6

# PageRank
PAGERANK_ALPHA = 0.85

# DGIdb druggable categories
DRUGGABLE_CATEGORIES = [
    "KINASE", "ENZYME", "ION CHANNEL", "TRANSPORTER",
    "G PROTEIN COUPLED RECEPTOR", "NUCLEAR HORMONE RECEPTOR",
    "PROTEASE", "PHOSPHATASE",
]


# ============================================================
# STEP 1: Download interazioni PPI
# ============================================================
def get_string_interactions(genes: list[str], species: int = STRING_SPECIES,
                            score_threshold: int = STRING_SCORE_THRESHOLD) -> pd.DataFrame:
    """Scarica interazioni da STRING per una lista di geni."""
    print(f"[Fase2] Querying STRING for {len(genes)} genes...")

    # STRING accetta max ~2000 proteine per query
    # Split in batch se necessario
    batch_size = 500
    all_interactions = []

    for i in range(0, len(genes), batch_size):
        batch = genes[i:i + batch_size]
        url = "https://string-db.org/api/json/network"
        params = {
            "identifiers": "%0d".join(batch),
            "species": species,
            "required_score": score_threshold,
            "caller_identity": "magic_solution_project",
        }

        try:
            response = requests.get(url, params=params, timeout=60)
            response.raise_for_status()
            data = response.json()
            if data:
                all_interactions.append(pd.DataFrame(data))
            print(f"  -> Batch {i//batch_size + 1}: {len(data)} interactions")
        except Exception as e:
            print(f"  [WARN] STRING query failed for batch {i//batch_size + 1}: {e}")

        time.sleep(1)  # Rate limiting

    if all_interactions:
        return pd.concat(all_interactions).drop_duplicates(
            subset=["preferredName_A", "preferredName_B"]
        )
    return pd.DataFrame()


def get_biogrid_interactions(genes: list[str],
                              biogrid_file: str = None) -> pd.DataFrame:
    """Scarica/carica interazioni da BioGRID.

    BioGRID offre download completo (tab-delimited) o API REST.
    Per dataset grandi, preferire il file completo:
    https://downloads.thebiogrid.org/BioGRID
    """
    if biogrid_file and Path(biogrid_file).exists():
        print(f"[Fase2] Loading BioGRID from {biogrid_file}...")
        bg = pd.read_csv(biogrid_file, sep="\t", low_memory=False)
        # Filtra per geni di interesse e specie
        gene_set = set(genes)
        mask = (
            bg["Official Symbol Interactor A"].isin(gene_set) |
            bg["Official Symbol Interactor B"].isin(gene_set)
        )
        filtered = bg[mask][["Official Symbol Interactor A",
                              "Official Symbol Interactor B",
                              "Score"]]
        return filtered.rename(columns={
            "Official Symbol Interactor A": "gene_a",
            "Official Symbol Interactor B": "gene_b",
        })
    else:
        print("[Fase2] BioGRID file not provided. Using STRING only.")
        return pd.DataFrame()


# ============================================================
# STEP 2: Costruzione grafo PPI
# ============================================================
def build_ppi_network(string_df: pd.DataFrame,
                       biogrid_df: pd.DataFrame = None) -> nx.Graph:
    """Costruisce un grafo PPI unificato da STRING + BioGRID."""
    print("[Fase2] Building PPI network...")
    G = nx.Graph()

    # Aggiungi archi STRING
    if not string_df.empty:
        for _, row in string_df.iterrows():
            G.add_edge(
                row["preferredName_A"],
                row["preferredName_B"],
                weight=row["score"] / 1000,
                source="STRING",
            )

    # Aggiungi archi BioGRID
    if biogrid_df is not None and not biogrid_df.empty:
        for _, row in biogrid_df.iterrows():
            if G.has_edge(row["gene_a"], row["gene_b"]):
                # Se l'arco esiste gia', aggiorna il peso (max)
                existing_weight = G[row["gene_a"]][row["gene_b"]]["weight"]
                G[row["gene_a"]][row["gene_b"]]["weight"] = max(
                    existing_weight, 0.7  # default BioGRID weight
                )
            else:
                G.add_edge(row["gene_a"], row["gene_b"],
                           weight=0.7, source="BioGRID")

    print(f"  -> Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Rimuovi nodi isolati
    isolates = list(nx.isolates(G))
    G.remove_nodes_from(isolates)
    print(f"  -> After removing isolates: {G.number_of_nodes()} nodes")

    return G


# ============================================================
# STEP 3: Network diffusion algorithms
# ============================================================
def personalized_pagerank(G: nx.Graph, seed_genes: list[str],
                           alpha: float = PAGERANK_ALPHA) -> dict:
    """PageRank personalizzato con seed genes come sorgente."""
    print("[Fase2] Running Personalized PageRank...")
    seeds_in_network = [g for g in seed_genes if g in G]
    if not seeds_in_network:
        print("  [WARN] No seed genes found in network!")
        return {}

    personalization = {g: 1.0 / len(seeds_in_network) for g in seeds_in_network}
    scores = nx.pagerank(G, alpha=alpha, personalization=personalization,
                          weight="weight")
    print(f"  -> Scores computed for {len(scores)} nodes")
    return scores


def random_walk_restart(G: nx.Graph, seed_genes: list[str],
                         restart_prob: float = RWR_RESTART_PROB,
                         max_iter: int = RWR_MAX_ITER,
                         tol: float = RWR_TOL) -> dict:
    """Random Walk with Restart (implementazione matriciale)."""
    print("[Fase2] Running Random Walk with Restart...")
    nodes = list(G.nodes())
    node_idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)

    # Matrice di adiacenza pesata
    A = nx.adjacency_matrix(G, nodelist=nodes, weight="weight").toarray()

    # Matrice di transizione (column-stochastic)
    col_sums = A.sum(axis=0)
    col_sums[col_sums == 0] = 1  # Evita divisione per zero
    W = A / col_sums  # Normalizza per colonna

    # Vettore iniziale (seed genes)
    p0 = np.zeros(n)
    seeds_in_network = [g for g in seed_genes if g in node_idx]
    for g in seeds_in_network:
        p0[node_idx[g]] = 1.0 / len(seeds_in_network)

    # Iterazione
    p = p0.copy()
    for iteration in range(max_iter):
        p_new = (1 - restart_prob) * (W @ p) + restart_prob * p0
        diff = np.linalg.norm(p_new - p, 1)
        p = p_new
        if diff < tol:
            print(f"  -> Converged after {iteration + 1} iterations")
            break
    else:
        print(f"  -> Max iterations ({max_iter}) reached")

    return dict(zip(nodes, p))


# ============================================================
# STEP 4: Druggability filtering
# ============================================================
def query_dgidb(genes: list[str]) -> pd.DataFrame:
    """Interroga DGIdb per informazioni sulla druggability."""
    print(f"[Fase2] Querying DGIdb for {len(genes)} genes...")

    # DGIdb API v2
    gene_str = ",".join(genes[:100])  # Limita a 100 per query
    url = f"https://dgidb.org/api/v2/interactions.json?genes={gene_str}"

    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        data = response.json()

        results = []
        for match in data.get("matchedTerms", []):
            gene = match["searchTerm"]
            categories = match.get("geneCategories", [])
            for interaction in match.get("interactions", []):
                results.append({
                    "gene": gene,
                    "drug_name": interaction.get("drugName"),
                    "interaction_type": interaction.get("interactionTypes"),
                    "categories": ", ".join(categories),
                })
        return pd.DataFrame(results)

    except Exception as e:
        print(f"  [WARN] DGIdb query failed: {e}")
        return pd.DataFrame()


def filter_druggable(ranking_df: pd.DataFrame,
                      dgidb_df: pd.DataFrame) -> pd.DataFrame:
    """Filtra target per druggability."""
    if dgidb_df.empty:
        return ranking_df

    druggable_genes = set(dgidb_df["gene"].unique())
    ranking_df["is_druggable"] = ranking_df.index.isin(druggable_genes)
    ranking_df["dgidb_categories"] = ranking_df.index.map(
        dgidb_df.groupby("gene")["categories"].first()
    )
    return ranking_df


# ============================================================
# STEP 5: LLM-based literature refinement
# ============================================================
def pubmed_text_mining(genes: list[str], top_n: int = 20):
    """Scarica abstract PubMed e analizzali con un LLM.

    Questo step e' opzionale e richiede accesso ad API LLM.
    """
    try:
        from pymed import PubMed
    except ImportError:
        print("[Fase2] pymed not installed. Skipping text mining.")
        return {}

    pubmed = PubMed(tool="MagicSolution", email="your@email.com")

    gene_evidence = {}
    for gene in genes[:top_n]:
        query = f'"{gene}" AND ("ischemia reperfusion" OR "acute kidney injury")'
        articles = list(pubmed.query(query, max_results=5))
        abstracts = [a.abstract for a in articles if a.abstract]

        if abstracts:
            # TODO: Invia ad API LLM per analisi
            # prompt = f"Analizza il ruolo di {gene} nel danno renale IR basandoti su: {abstracts}"
            # response = llm_client.generate(prompt)
            gene_evidence[gene] = {
                "n_articles": len(abstracts),
                "abstracts_found": True,
            }
        else:
            gene_evidence[gene] = {"n_articles": 0, "abstracts_found": False}

    return gene_evidence


# ============================================================
# MAIN PIPELINE
# ============================================================
def main():
    print("=" * 60)
    print("MAGIC SOLUTION - Fase 2: Network Medicine")
    print("=" * 60)

    # --- Carica seed genes ---
    seed_file = INPUT_DIR / "seed_genes.csv"
    if seed_file.exists():
        seed_df = pd.read_csv(seed_file, index_col=0)
        seed_genes = seed_df.index.tolist()
        print(f"\nLoaded {len(seed_genes)} seed genes from Fase 1")
    else:
        print(f"\n[WARN] {seed_file} not found. Using demo genes.")
        # Demo seed genes (pathway IR-AKI noti dalla letteratura)
        seed_genes = [
            "TNF", "IL6", "IL1B", "NFKB1", "CASP3", "BCL2", "HMGB1",
            "TLR4", "VEGFA", "HIF1A", "SOD2", "NOS2", "PTGS2",
            "MMP9", "ICAM1", "CCL2", "CXCL8", "IL10", "TGFB1",
        ]

    # --- Scarica interazioni PPI ---
    string_df = get_string_interactions(seed_genes)
    biogrid_df = get_biogrid_interactions(seed_genes)

    # --- Costruisci rete ---
    G = build_ppi_network(string_df, biogrid_df)
    nx.write_graphml(G, str(OUTPUT_DIR / "ppi_network.graphml"))

    # --- Network diffusion ---
    pr_scores = personalized_pagerank(G, seed_genes)
    rwr_scores = random_walk_restart(G, seed_genes)

    # --- Combina score ---
    combined = pd.DataFrame({
        "pagerank": pd.Series(pr_scores),
        "rwr": pd.Series(rwr_scores),
    }).fillna(0)

    # Rank combination
    combined["pr_rank"] = combined["pagerank"].rank(ascending=False)
    combined["rwr_rank"] = combined["rwr"].rank(ascending=False)
    combined["combined_rank"] = (combined["pr_rank"] + combined["rwr_rank"]) / 2
    combined = combined.sort_values("combined_rank")

    print(f"\n--- Top 20 targets by combined rank ---")
    print(combined.head(20)[["pagerank", "rwr", "combined_rank"]])

    # --- Druggability check ---
    top_genes = combined.head(50).index.tolist()
    dgidb_results = query_dgidb(top_genes)
    combined = filter_druggable(combined, dgidb_results)

    # --- Text mining (opzionale) ---
    # evidence = pubmed_text_mining(top_genes)

    # --- Salva risultati ---
    combined.to_csv(OUTPUT_DIR / "ranked_targets.csv")
    if not dgidb_results.empty:
        dgidb_results.to_csv(OUTPUT_DIR / "dgidb_interactions.csv")

    print(f"\n[Fase2] Completata. Output in {OUTPUT_DIR}")
    print("  -> ranked_targets.csv")
    print("  -> ppi_network.graphml")


if __name__ == "__main__":
    main()
