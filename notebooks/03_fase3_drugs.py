"""
MAGIC SOLUTION - Fase 3: Drug Repurposing
==========================================
Interroga database farmacologici per trovare farmaci approvati
che interagiscano con i target prioritizzati, e li filtra per
idoneita' alla formulazione in soluzione acquosa.

Requisiti:
    pip install chembl_webresource_client rdkit-pypi requests pandas lxml

Input:  ranked_targets.csv (dalla Fase 2)
Output: candidate_drugs.csv, drug_target_matrix.csv
"""

import pandas as pd
import numpy as np
import requests
import time
from pathlib import Path

# ============================================================
# CONFIGURAZIONE
# ============================================================
INPUT_DIR = Path("../output/fase2")
OUTPUT_DIR = Path("../output/fase3")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Numero massimo di target da interrogare
MAX_TARGETS = 30

# Soglie di potenza (pChEMBL >= 6 = IC50 <= 1 uM)
PCHEMBL_THRESHOLD = 6.0

# Filtri per soluzione acquosa
MAX_LOGP = 2.0          # Solubilit\u00E0 in acqua
MAX_MW = 500             # Peso molecolare
MAX_PHASE = 4            # Farmaco approvato (fase 4)
MIN_PHASE = 3            # Almeno fase 3 clinica


# ============================================================
# STEP 1: Query ChEMBL
# ============================================================
def query_chembl_for_target(gene_name: str) -> pd.DataFrame:
    """Cerca farmaci attivi su un target in ChEMBL."""
    from chembl_webresource_client.new_client import new_client

    target_client = new_client.target
    activity_client = new_client.activity

    results = []

    try:
        # Cerca il target
        targets = target_client.search(gene_name).filter(
            organism="Homo sapiens",
            target_type="SINGLE PROTEIN",
        )

        for t in targets[:3]:  # Max 3 target ChEMBL per gene
            chembl_id = t["target_chembl_id"]

            # Cerca attivita' con buona potenza
            activities = activity_client.filter(
                target_chembl_id=chembl_id,
                standard_type__in=["IC50", "Ki", "EC50", "Kd"],
                pchembl_value__gte=PCHEMBL_THRESHOLD,
            ).only([
                "molecule_chembl_id", "canonical_smiles",
                "standard_type", "standard_value", "standard_units",
                "pchembl_value", "target_chembl_id",
            ])

            for act in activities:
                results.append({
                    "gene": gene_name,
                    "target_chembl_id": chembl_id,
                    "molecule_chembl_id": act.get("molecule_chembl_id"),
                    "smiles": act.get("canonical_smiles"),
                    "activity_type": act.get("standard_type"),
                    "activity_value": act.get("standard_value"),
                    "activity_units": act.get("standard_units"),
                    "pchembl": act.get("pchembl_value"),
                })

    except Exception as e:
        print(f"  [WARN] ChEMBL query failed for {gene_name}: {e}")

    return pd.DataFrame(results)


def get_molecule_info(chembl_id: str) -> dict:
    """Recupera info dettagliate su una molecola da ChEMBL."""
    from chembl_webresource_client.new_client import new_client

    molecule_client = new_client.molecule

    try:
        mol = molecule_client.get(chembl_id)
        if mol:
            props = mol.get("molecule_properties", {}) or {}
            return {
                "pref_name": mol.get("pref_name"),
                "max_phase": mol.get("max_phase"),
                "molecule_type": mol.get("molecule_type"),
                "oral": mol.get("oral"),
                "alogp": props.get("alogp"),
                "mw_freebase": props.get("mw_freebase"),
                "hba": props.get("hba"),
                "hbd": props.get("hbd"),
                "psa": props.get("psa"),
                "ro5_violations": props.get("num_ro5_violations"),
            }
    except Exception:
        pass
    return {}


# ============================================================
# STEP 2: Query DGIdb
# ============================================================
def query_dgidb_drugs(gene: str) -> list[dict]:
    """Interroga DGIdb per interazioni farmaco-gene."""
    url = f"https://dgidb.org/api/v2/interactions.json?genes={gene}"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        results = []
        for match in data.get("matchedTerms", []):
            for interaction in match.get("interactions", []):
                results.append({
                    "gene": gene,
                    "drug_name": interaction.get("drugName"),
                    "interaction_type": interaction.get("interactionTypes"),
                    "source": ", ".join(interaction.get("sources", [])),
                    "pmids": ", ".join(interaction.get("pmids", [])),
                })
        return results
    except Exception:
        return []


# ============================================================
# STEP 3: Filtri per soluzione di perfusione
# ============================================================
def compute_physicochemical_properties(smiles: str) -> dict:
    """Calcola proprieta' fisico-chimiche con RDKit."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Crippen

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"valid": False}

        return {
            "valid": True,
            "MW": round(Descriptors.MolWt(mol), 2),
            "LogP": round(Crippen.MolLogP(mol), 2),
            "TPSA": round(Descriptors.TPSA(mol), 2),
            "HBD": Descriptors.NumHDonors(mol),
            "HBA": Descriptors.NumHAcceptors(mol),
            "RotBonds": Descriptors.NumRotatableBonds(mol),
            "Lipinski_violations": sum([
                Descriptors.MolWt(mol) > 500,
                Crippen.MolLogP(mol) > 5,
                Descriptors.NumHDonors(mol) > 5,
                Descriptors.NumHAcceptors(mol) > 10,
            ]),
        }
    except ImportError:
        print("[WARN] RDKit not available")
        return {"valid": False}


def filter_for_perfusion_solution(drugs_df: pd.DataFrame) -> pd.DataFrame:
    """Applica filtri specifici per soluzione di perfusione acquosa."""
    print("[Fase3] Applying perfusion solution filters...")

    # 1. Solubilit\u00E0 in acqua (LogP basso)
    mask_soluble = drugs_df["LogP"] < MAX_LOGP
    print(f"  -> After LogP < {MAX_LOGP}: {mask_soluble.sum()} molecules")

    # 2. Peso molecolare ragionevole
    mask_mw = drugs_df["MW"] < MAX_MW
    print(f"  -> After MW < {MAX_MW}: {(mask_soluble & mask_mw).sum()} molecules")

    # 3. Status regolatorio (farmaco approvato o fase 3+)
    mask_approved = drugs_df["max_phase"] >= MIN_PHASE
    print(f"  -> After phase >= {MIN_PHASE}: {(mask_soluble & mask_mw & mask_approved).sum()} molecules")

    filtered = drugs_df[mask_soluble & mask_mw & mask_approved]
    return filtered


# ============================================================
# STEP 4: Connectivity Map (Signature Reversal)
# ============================================================
def signature_reversal_cmap(up_genes: list[str], down_genes: list[str]):
    """Interroga CMap/LINCS per signature reversal.

    Trova farmaci il cui profilo trascrizionale inverte la firma del danno IR.
    Richiede: up_genes (sovra-espressi in IR) e down_genes (sotto-espressi).
    """
    # CMap/LINCS L1000 API (CLUE.io)
    # Nota: richiede registrazione su https://clue.io
    print("[Fase3] CMap signature reversal - placeholder")
    print("  -> Registrarsi su https://clue.io per accesso API")
    print("  -> Alternativa: usare cmapPy per analisi locale")

    # Esempio di query CLUE.io
    # url = "https://api.clue.io/api/sigs"
    # headers = {"user_key": "YOUR_API_KEY"}
    # query = {
    #     "where": {"pert_type": "trt_cp"},
    #     "fields": ["pert_iname", "score"]
    # }

    return pd.DataFrame()  # Placeholder


# ============================================================
# STEP 5: Selezione finale del cocktail
# ============================================================
def select_cocktail(candidates: pd.DataFrame,
                     n_drugs: int = 5) -> pd.DataFrame:
    """Seleziona il cocktail finale ottimizzando per:
    - Copertura di pathway diversi
    - Massima potenza (pChEMBL)
    - Status approvato (max_phase = 4)
    - Minimo overlap di target
    """
    print(f"[Fase3] Selecting top {n_drugs} candidates for cocktail...")

    if candidates.empty:
        return candidates

    # Ordina per: approvato > potenza > solubilit\u00E0
    candidates["score"] = (
        candidates["max_phase"].fillna(0) * 10 +
        candidates["pchembl"].fillna(0) +
        (MAX_LOGP - candidates["LogP"].fillna(0))
    )

    # Seleziona diversificando per gene target
    selected = []
    used_genes = set()
    for _, row in candidates.sort_values("score", ascending=False).iterrows():
        if row["gene"] not in used_genes:
            selected.append(row)
            used_genes.add(row["gene"])
            if len(selected) >= n_drugs:
                break

    return pd.DataFrame(selected)


# ============================================================
# MAIN PIPELINE
# ============================================================
def main():
    print("=" * 60)
    print("MAGIC SOLUTION - Fase 3: Drug Repurposing")
    print("=" * 60)

    # --- Carica target ranked ---
    target_file = INPUT_DIR / "ranked_targets.csv"
    if target_file.exists():
        targets = pd.read_csv(target_file, index_col=0)
        top_targets = targets.head(MAX_TARGETS).index.tolist()
        print(f"\nLoaded {len(top_targets)} top targets from Fase 2")
    else:
        print(f"\n[WARN] {target_file} not found. Using demo targets.")
        top_targets = ["TNF", "IL6", "CASP3", "HMGB1", "TLR4",
                        "HIF1A", "NOS2", "PTGS2", "MMP9", "NFKB1"]

    # --- Query ChEMBL per ogni target ---
    all_chembl_drugs = []
    for i, gene in enumerate(top_targets):
        print(f"\n[{i+1}/{len(top_targets)}] Querying ChEMBL for {gene}...")
        drugs = query_chembl_for_target(gene)
        if not drugs.empty:
            all_chembl_drugs.append(drugs)
            print(f"  -> Found {len(drugs)} activities")
        time.sleep(0.5)

    if all_chembl_drugs:
        chembl_df = pd.concat(all_chembl_drugs).drop_duplicates("molecule_chembl_id")
        print(f"\nTotal unique molecules from ChEMBL: {len(chembl_df)}")
    else:
        chembl_df = pd.DataFrame()
        print("\n[WARN] No ChEMBL results found")

    # --- Proprieta' fisico-chimiche ---
    if not chembl_df.empty and "smiles" in chembl_df.columns:
        print("\n[Fase3] Computing physicochemical properties...")
        props = chembl_df["smiles"].dropna().apply(compute_physicochemical_properties)
        props_df = pd.DataFrame(props.tolist(), index=props.index)
        chembl_df = pd.concat([chembl_df, props_df], axis=1)

        # Recupera info molecola (max_phase, nome, ecc.)
        print("[Fase3] Fetching molecule details...")
        for idx, row in chembl_df.iterrows():
            mol_info = get_molecule_info(row["molecule_chembl_id"])
            for k, v in mol_info.items():
                chembl_df.loc[idx, k] = v
            time.sleep(0.3)

    # --- Filtra per soluzione di perfusione ---
    if not chembl_df.empty:
        filtered = filter_for_perfusion_solution(chembl_df)
    else:
        filtered = pd.DataFrame()

    # --- Query DGIdb per conferma ---
    dgidb_results = []
    for gene in top_targets[:10]:
        dgidb_results.extend(query_dgidb_drugs(gene))
        time.sleep(0.3)
    dgidb_df = pd.DataFrame(dgidb_results)

    # --- Seleziona cocktail ---
    cocktail = select_cocktail(filtered)

    # --- Salva risultati ---
    if not chembl_df.empty:
        chembl_df.to_csv(OUTPUT_DIR / "all_chembl_hits.csv", index=False)
    if not filtered.empty:
        filtered.to_csv(OUTPUT_DIR / "filtered_candidates.csv", index=False)
    if not dgidb_df.empty:
        dgidb_df.to_csv(OUTPUT_DIR / "dgidb_interactions.csv", index=False)
    if not cocktail.empty:
        cocktail.to_csv(OUTPUT_DIR / "candidate_drugs.csv", index=False)

    print(f"\n[Fase3] Completata. Output in {OUTPUT_DIR}")
    print("  -> all_chembl_hits.csv")
    print("  -> filtered_candidates.csv")
    print("  -> candidate_drugs.csv (cocktail finale)")


if __name__ == "__main__":
    main()
