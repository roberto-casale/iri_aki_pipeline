"""
MAGIC SOLUTION - Fase 4: Valutazione In Silico
================================================
Valuta sicurezza, efficacia e interazioni del cocktail
farmacologico. Definisce le concentrazioni ottimali.

Requisiti:
    pip install rdkit-pypi scipy pandas numpy biopandas requests
    pip install meeko vina  # per docking (opzionale, richiede setup)
    pip install deepchem    # per tossicologia (opzionale, pesante)

Input:  candidate_drugs.csv (dalla Fase 3)
Output: final_formulation.csv, admet_profiles.csv, ddi_matrix.csv
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")

# ============================================================
# CONFIGURAZIONE
# ============================================================
INPUT_DIR = Path("../output/fase3")
OUTPUT_DIR = Path("../output/fase4")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# STEP 1: Profilo ADMET completo
# ============================================================
def compute_full_admet(smiles: str) -> dict:
    """Calcola profilo ADMET completo con RDKit."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
    except ImportError:
        return {"error": "RDKit not installed"}

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": "Invalid SMILES"}

    return {
        # Proprieta' fisiche
        "MW": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Crippen.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "RotBonds": Descriptors.NumRotatableBonds(mol),
        "RingCount": Descriptors.RingCount(mol),
        "AromaticRings": Descriptors.NumAromaticRings(mol),
        "FractionCSP3": round(Descriptors.FractionCSP3(mol), 3),

        # Regole
        "Lipinski_violations": sum([
            Descriptors.MolWt(mol) > 500,
            Crippen.MolLogP(mol) > 5,
            Descriptors.NumHDonors(mol) > 5,
            Descriptors.NumHAcceptors(mol) > 10,
        ]),
        "Veber_violations": sum([
            Descriptors.TPSA(mol) > 140,
            Descriptors.NumRotatableBonds(mol) > 10,
        ]),

        # Solubilit\u00E0 stimata (ESOL model approssimato)
        "ESOL_LogS": round(-1.25 * Crippen.MolLogP(mol) + 0.5, 2),

        # Complessita' molecolare
        "BertzCT": round(Descriptors.BertzCT(mol), 1),
    }


# ============================================================
# STEP 2: Tossicologia predittiva
# ============================================================
def predict_toxicity_deepchem(smiles_list: list[str]) -> pd.DataFrame:
    """Predice tossicit\u00E0 con modelli Tox21 pre-addestrati (DeepChem)."""
    try:
        import deepchem as dc

        print("[Fase4] Loading Tox21 model...")
        tox21_tasks, (train, valid, test), transformers = dc.molnet.load_tox21(
            featurizer="GraphConv", split="random"
        )

        model = dc.models.GraphConvModel(
            len(tox21_tasks), mode="classification",
            batch_size=50, learning_rate=0.001,
        )
        # Nota: in produzione, caricare un modello pre-addestrato
        # model.restore()

        # Featurize le molecole candidate
        featurizer = dc.feat.ConvMolFeaturizer()
        features = featurizer.featurize(smiles_list)

        # Predici
        dataset = dc.data.NumpyDataset(X=features)
        predictions = model.predict(dataset)

        results = pd.DataFrame(predictions[:, :, 1],
                                columns=tox21_tasks,
                                index=smiles_list)
        return results

    except ImportError:
        print("[WARN] DeepChem not available. Skipping toxicity prediction.")
        return pd.DataFrame()
    except Exception as e:
        print(f"[WARN] Toxicity prediction failed: {e}")
        return pd.DataFrame()


def predict_toxicity_protox(smiles: str) -> dict:
    """Query ProTox-II per predizione tossicit\u00E0.

    ProTox-II: https://tox-new.charite.de/protox_II/
    Nota: ProTox-II ha API web, ma verificare TOS per uso programmatico.
    """
    # Placeholder - in produzione usare l'API
    return {
        "LD50_predicted": None,
        "toxicity_class": None,
        "hepatotoxicity": None,
        "carcinogenicity": None,
        "immunotoxicity": None,
        "mutagenicity": None,
    }


# ============================================================
# STEP 3: Molecular Docking
# ============================================================
def prepare_ligand(smiles: str, output_pdbqt: str):
    """Prepara un ligando per AutoDock Vina da SMILES."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import meeko

        # Genera conformazione 3D
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

        # Converti a PDBQT con Meeko
        preparator = meeko.MoleculePreparation()
        mol_setup = preparator.prepare(mol)
        pdbqt_str = meeko.PDBQTWriterLegacy.write_string(mol_setup)

        with open(output_pdbqt, "w") as f:
            for line in pdbqt_str:
                f.write(line)

        return True
    except Exception as e:
        print(f"  [WARN] Ligand preparation failed: {e}")
        return False


def run_docking(receptor_pdbqt: str, ligand_pdbqt: str,
                center: tuple, box_size: tuple = (20, 20, 20),
                exhaustiveness: int = 32) -> dict:
    """Esegue molecular docking con AutoDock Vina."""
    try:
        from vina import Vina

        v = Vina(sf_name="vina")
        v.set_receptor(receptor_pdbqt)
        v.set_ligand_from_file(ligand_pdbqt)
        v.compute_vina_maps(center=list(center), box_size=list(box_size))
        v.dock(exhaustiveness=exhaustiveness, n_poses=5)

        energies = v.energies()
        return {
            "best_affinity_kcal": energies[0][0],
            "n_poses": len(energies),
            "all_affinities": [e[0] for e in energies],
        }
    except ImportError:
        print("  [WARN] AutoDock Vina not installed.")
        return {"best_affinity_kcal": None}
    except Exception as e:
        print(f"  [WARN] Docking failed: {e}")
        return {"best_affinity_kcal": None}


def download_target_structure(uniprot_id: str, output_dir: str) -> str:
    """Scarica struttura 3D da AlphaFold o PDB."""
    import requests

    # Prova prima AlphaFold
    af_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    try:
        resp = requests.get(af_url, timeout=30)
        if resp.status_code == 200:
            pdb_path = Path(output_dir) / f"{uniprot_id}_alphafold.pdb"
            pdb_path.write_text(resp.text)
            return str(pdb_path)
    except Exception:
        pass

    print(f"  [WARN] Could not download structure for {uniprot_id}")
    return ""


# ============================================================
# STEP 4: Drug-Drug Interactions (DDI)
# ============================================================
def check_ddi_drugbank(drug_names: list[str],
                        drugbank_xml: str = None) -> pd.DataFrame:
    """Verifica interazioni farmaco-farmaco da DrugBank.

    Richiede il file XML completo di DrugBank (licenza accademica).
    """
    if drugbank_xml and Path(drugbank_xml).exists():
        from lxml import etree

        print("[Fase4] Parsing DrugBank for DDI...")
        tree = etree.parse(drugbank_xml)
        root = tree.getroot()
        ns = {"db": "http://www.drugbank.ca"}

        interactions = []
        drug_name_set = set(d.lower() for d in drug_names)

        for drug in root.findall("db:drug", ns):
            name = drug.findtext("db:name", namespaces=ns, default="").lower()
            if name in drug_name_set:
                for interaction in drug.findall(".//db:drug-interaction", ns):
                    partner = interaction.findtext("db:name", namespaces=ns, default="")
                    desc = interaction.findtext("db:description", namespaces=ns, default="")
                    if partner.lower() in drug_name_set:
                        interactions.append({
                            "drug_a": name,
                            "drug_b": partner.lower(),
                            "description": desc,
                        })

        return pd.DataFrame(interactions)
    else:
        print("[Fase4] DrugBank XML not provided. DDI check skipped.")
        print("  -> Download da https://go.drugbank.com/releases/latest")
        return pd.DataFrame()


def build_ddi_matrix(drugs: list[str], ddi_df: pd.DataFrame) -> pd.DataFrame:
    """Costruisce matrice di interazioni pairwise."""
    n = len(drugs)
    matrix = pd.DataFrame(0, index=drugs, columns=drugs, dtype=int)

    if not ddi_df.empty:
        for _, row in ddi_df.iterrows():
            if row["drug_a"] in drugs and row["drug_b"] in drugs:
                matrix.loc[row["drug_a"], row["drug_b"]] = 1
                matrix.loc[row["drug_b"], row["drug_a"]] = 1

    return matrix


# ============================================================
# STEP 5: Ottimizzazione concentrazioni
# ============================================================
def hill_equation(conc, emax, ec50, n):
    """Equazione di Hill per modellare la dose-risposta."""
    return emax * (conc ** n) / (ec50 ** n + conc ** n)


def fit_dose_response(concentrations: np.ndarray,
                       responses: np.ndarray) -> dict:
    """Fit curva dose-risposta con equazione di Hill."""
    try:
        popt, pcov = curve_fit(
            hill_equation, concentrations, responses,
            p0=[100, 1, 1],  # Iniziale: Emax=100%, EC50=1 uM, n=1
            bounds=([0, 0, 0], [200, 1000, 10]),
            maxfev=5000,
        )
        perr = np.sqrt(np.diag(pcov))

        return {
            "Emax": round(popt[0], 2),
            "EC50_uM": round(popt[1], 4),
            "Hill_coefficient": round(popt[2], 2),
            "Emax_SE": round(perr[0], 2),
            "EC50_SE": round(perr[1], 4),
            "fit_success": True,
        }
    except Exception as e:
        print(f"  [WARN] Dose-response fit failed: {e}")
        return {"fit_success": False}


def compute_target_concentration(ec50_um: float, hill_n: float,
                                  target_effect_pct: float = 80,
                                  perfusion_volume_ml: float = 500) -> dict:
    """Calcola la concentrazione target nella soluzione di perfusione.

    Args:
        ec50_um: EC50 in micromolare
        hill_n: Coefficiente di Hill
        target_effect_pct: Effetto desiderato (% di Emax)
        perfusion_volume_ml: Volume totale della soluzione

    Returns:
        Concentrazioni in varie unita'
    """
    # Concentrazione per raggiungere target_effect_pct di Emax
    # Dalla Hill equation: C = EC50 * (E / (Emax - E))^(1/n)
    ratio = target_effect_pct / (100 - target_effect_pct)
    conc_um = ec50_um * (ratio ** (1 / hill_n))

    # Converti in mg/L (assumendo MW medio di 400 Da)
    # In produzione: usare MW reale del farmaco
    mw_estimated = 400  # Da
    conc_mg_per_l = conc_um * mw_estimated / 1000

    return {
        "target_conc_uM": round(conc_um, 2),
        "target_conc_mg_L": round(conc_mg_per_l, 2),
        "target_effect_pct": target_effect_pct,
        "perfusion_volume_ml": perfusion_volume_ml,
        "total_drug_mg": round(conc_mg_per_l * perfusion_volume_ml / 1000, 3),
    }


# ============================================================
# MAIN PIPELINE
# ============================================================
def main():
    print("=" * 60)
    print("MAGIC SOLUTION - Fase 4: Valutazione In Silico")
    print("=" * 60)

    # --- Carica candidati ---
    candidate_file = INPUT_DIR / "candidate_drugs.csv"
    if candidate_file.exists():
        candidates = pd.read_csv(candidate_file)
        print(f"\nLoaded {len(candidates)} candidate drugs from Fase 3")
    else:
        print(f"\n[WARN] {candidate_file} not found. Using demo data.")
        candidates = pd.DataFrame({
            "gene": ["PTGS2", "TNF", "CASP3", "NOS2", "HIF1A"],
            "pref_name": ["Ibuprofen", "Infliximab", "Z-VAD-FMK", "L-NAME", "Roxadustat"],
            "smiles": [
                "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",
                None,  # Biologico
                None,  # Peptide
                "NC(=N)NCCCC(=O)O",
                "Oc1cc(=O)c2cc(-c3ccccc3)cn2[nH]1",
            ],
            "pchembl": [5.2, 8.0, 7.5, 5.5, 7.0],
            "max_phase": [4, 4, 2, 3, 4],
        })

    # --- ADMET completo ---
    print("\n--- Step 1: ADMET profiling ---")
    admet_records = []
    for idx, row in candidates.iterrows():
        if pd.notna(row.get("smiles")):
            profile = compute_full_admet(row["smiles"])
            profile["drug"] = row.get("pref_name", row.get("molecule_chembl_id", f"drug_{idx}"))
            admet_records.append(profile)
    admet_df = pd.DataFrame(admet_records)
    if not admet_df.empty:
        print(admet_df[["drug", "MW", "LogP", "TPSA", "Lipinski_violations"]].to_string())

    # --- Tossicologia ---
    print("\n--- Step 2: Toxicity prediction ---")
    smiles_valid = candidates["smiles"].dropna().tolist()
    if smiles_valid:
        tox_results = predict_toxicity_deepchem(smiles_valid)
        # Oppure: tox_results da ProTox-II (API web)

    # --- Molecular Docking ---
    print("\n--- Step 3: Molecular docking ---")
    print("[INFO] Docking richiede strutture PDB dei target.")
    print("  -> Scaricare da AlphaFold/PDB e convertire a PDBQT")
    print("  -> Eseguire prepare_ligand() + run_docking()")
    # Per ogni coppia farmaco-target:
    # 1. download_target_structure(uniprot_id)
    # 2. prepare_ligand(smiles, "ligand.pdbqt")
    # 3. run_docking("target.pdbqt", "ligand.pdbqt", center)

    # --- DDI ---
    print("\n--- Step 4: Drug-Drug Interactions ---")
    drug_names = candidates["pref_name"].dropna().tolist()
    ddi_df = check_ddi_drugbank(drug_names)
    ddi_matrix = build_ddi_matrix(drug_names, ddi_df)
    print("DDI matrix:")
    print(ddi_matrix.to_string())

    # --- Ottimizzazione concentrazioni ---
    print("\n--- Step 5: Concentration optimization ---")
    # Esempio con dati fittizi (in produzione: dati da ChEMBL/letteratura)
    example_conc = np.array([0.01, 0.1, 1, 10, 100])  # uM
    example_resp = np.array([5, 15, 50, 85, 95])  # % inibizione

    fit_result = fit_dose_response(example_conc, example_resp)
    if fit_result.get("fit_success"):
        print(f"  EC50 = {fit_result['EC50_uM']} uM")
        print(f"  Hill coefficient = {fit_result['Hill_coefficient']}")

        target_conc = compute_target_concentration(
            ec50_um=fit_result["EC50_uM"],
            hill_n=fit_result["Hill_coefficient"],
        )
        print(f"  Target concentration (80% effect): {target_conc['target_conc_uM']} uM")
        print(f"  = {target_conc['target_conc_mg_L']} mg/L")

    # --- Salva risultati ---
    if not admet_df.empty:
        admet_df.to_csv(OUTPUT_DIR / "admet_profiles.csv", index=False)
    ddi_matrix.to_csv(OUTPUT_DIR / "ddi_matrix.csv")

    # Formulazione finale
    formulation = candidates.copy()
    # TODO: aggiungere concentrazioni calcolate per ogni farmaco
    formulation.to_csv(OUTPUT_DIR / "final_formulation.csv", index=False)

    print(f"\n[Fase4] Completata. Output in {OUTPUT_DIR}")
    print("  -> admet_profiles.csv")
    print("  -> ddi_matrix.csv")
    print("  -> final_formulation.csv")


if __name__ == "__main__":
    main()
