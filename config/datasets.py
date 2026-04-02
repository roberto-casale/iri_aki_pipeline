"""
Configurazione dei dataset GEO per la pipeline Magic Solution.
Ogni entry contiene i metadati necessari per download e processing.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class Platform(Enum):
    ILLUMINA_HT12 = "illumina_ht12"       # Illumina HumanHT-12 V4.0
    AFFYMETRIX_GENE10ST = "affy_gene10st"  # Affymetrix Human Gene 1.0 ST
    AFFYMETRIX_GENE20ST = "affy_gene20st"  # Affymetrix Human Gene 2.0 ST
    AFFYMETRIX_U133PLUS2 = "affy_u133p2"   # Affymetrix HG-U133 Plus 2.0
    RNASEQ = "rnaseq"


class Role(Enum):
    DISCOVERY = "discovery"
    VALIDATION = "validation"
    SUPPORT = "support"


@dataclass
class DatasetConfig:
    accession: str
    title: str
    platform: Platform
    gpl: str
    n_samples: int
    organism: str
    role: Role
    description: str
    # Condizioni campione
    conditions: dict = field(default_factory=dict)
    # File supplementari disponibili su GEO
    supplementary_files: list = field(default_factory=list)
    # Come scaricare i dati processati
    download_method: str = "geoparse"  # "geoparse", "supplementary", "series_matrix"
    # Colonna nei metadati che indica la condizione IR/ctrl
    condition_column: Optional[str] = None
    # Note specifiche per il processing
    notes: str = ""


# ===========================================================
# DATASET DEFINITIONS
# ===========================================================

DATASETS = {
    # -------------------------------------------------------
    # DISCOVERY: Dataset principale (microarray, n=554)
    # -------------------------------------------------------
    "GSE43974": DatasetConfig(
        accession="GSE43974",
        title="Pathways for intervention to optimize donor organ quality",
        platform=Platform.ILLUMINA_HT12,
        gpl="GPL10558",
        n_samples=554,
        organism="Homo sapiens",
        role=Role.DISCOVERY,
        description=(
            "Biopsie renali da donatori viventi, DBD e DCD a 3 timepoint: "
            "T1=pre-retrieval, T2=post cold ischemia, T3=60min post reperfusion. "
            "DGF annotato per campioni post-reperfusion."
        ),
        conditions={
            "donor_types": ["BD", "DCD", "Living"],
            "timepoints": ["T1", "T2", "T3"],
            "n_BD_T1": 81, "n_DCD_T1": 37,
            "n_BD_T2_T3_paired": 67, "n_BD_T2_unpaired": 43, "n_BD_T3_unpaired": 38,
            "n_DCD_T2_T3_paired": 29, "n_DCD_T2_unpaired": 24, "n_DCD_T3_unpaired": 35,
            "n_Living_T1": 37, "n_Living_T3_paired": 34,
        },
        supplementary_files=[
            "GSE43974_non-normalized_data.txt.gz",
            "GSE43974_paired_analyses_DBD_DCD_Living.xls.gz",
        ],
        download_method="series_matrix",
        condition_column="title",  # Contiene info su donor type e timepoint
        notes=(
            "Dataset di riferimento nella letteratura IRI/DGF. "
            "Confronto chiave: T1 (ctrl) vs T3 (post-reperfusion IRI). "
            "I sample names seguono il pattern: Kidney_donor_{BD|DCD|Living}_T{1|2|3}_{n}. "
            "Series Matrix file contiene dati gia' normalizzati (quantile). "
            "Il file non-normalized_data.txt.gz contiene raw intensities."
        ),
    ),

    # -------------------------------------------------------
    # VALIDATION: RNA-seq datasets
    # -------------------------------------------------------
    "GSE90861": DatasetConfig(
        accession="GSE90861",
        title="A Molecular signature for Delayed Graft Function [RNA-Seq]",
        platform=Platform.RNASEQ,
        gpl="GPL18573",
        n_samples=46,
        organism="Homo sapiens",
        role=Role.VALIDATION,
        description=(
            "23 coppie di biopsie renali (pre-perfusione b vs post-perfusione b1). "
            "12 DGF + 11 IGF. Design paired."
        ),
        conditions={
            "n_DGF": 12, "n_IGF": 11,
            "biopsy_types": ["pre-perfusion (b)", "post-perfusion (b1)"],
        },
        supplementary_files=["GSE90861_RAW.tar"],
        download_method="supplementary",
        condition_column="title",  # Contiene IGF/DGF e b/b1
        notes=(
            "RAW.tar contiene file .txt per campione con counts. "
            "Sample names: {IGF|DGF}_{id}{b|b1}. "
            "b = pre-perfusion, b1 = post-perfusion."
        ),
    ),

    "GSE126805": DatasetConfig(
        accession="GSE126805",
        title="Transcriptional trajectories of human kidney injury progression",
        platform=Platform.RNASEQ,
        gpl="GPL21290",
        n_samples=163,
        organism="Homo sapiens",
        role=Role.VALIDATION,
        description=(
            "Biopsie protocollari da 42 trapianti renali a 4 timepoint: "
            "pre-reperfusion (baseline_pre), post-reperfusion (baseline_post), "
            "3 mesi, 12 mesi."
        ),
        conditions={
            "timepoints": ["baseline_pre", "baseline_post", "3months", "1year"],
            "n_patients": 42,
        },
        supplementary_files=[
            "GSE126805_all.gene_counts.txt.gz",
            "GSE126805_all.gene_RPKM.txt.gz",
        ],
        download_method="supplementary",
        condition_column="title",  # Contiene patient_id e timepoint
        notes=(
            "File counts gia' disponibile come matrice (gene x sample). "
            "Confronto chiave: baseline_pre vs baseline_post (= effetto IRI). "
            "I timepoint 3m e 12m servono per traiettorie di recupero."
        ),
    ),

    "GSE142077": DatasetConfig(
        accession="GSE142077",
        title="RNA-Seq identifies condition-specific biological signatures of IRI",
        platform=Platform.RNASEQ,
        gpl="GPL16791",
        n_samples=15,
        organism="Homo sapiens",
        role=Role.VALIDATION,
        description=(
            "5 pazienti maschi, ciascuno con 3 condizioni: "
            "pre-ischemia (P), ischemia 15min (I), reperfusione 10min (R10). "
            "Design paired."
        ),
        conditions={
            "n_patients": 5,
            "conditions": ["P (pre-ischemia)", "I (ischemia)", "R10 (reperfusion)"],
        },
        supplementary_files=["GSE142077_RAW.tar"],
        download_method="geoparse",
        condition_column="title",  # P-IRI-h-{n}, I-IRI-h-{n}, R10-IRI-h-{n}
        notes=(
            "Design paired eccellente. Campioni da nefrectomia per tumore. "
            "Sample naming: {P|I|R10}-IRI-h-{patient_id}. "
            "RAW.tar contiene count files per campione."
        ),
    ),

    "GSE293480": DatasetConfig(
        accession="GSE293480",
        title="Transcriptomic profiling during NMP of human kidneys [bulk RNA-seq]",
        platform=Platform.RNASEQ,
        gpl="GPL24676",
        n_samples=7,
        organism="Homo sapiens",
        role=Role.VALIDATION,
        description=(
            "Bulk RNA-seq da biopsie renali durante NMP. "
            "Confronto DGF prolungato vs IGF post-trapianto."
        ),
        conditions={
            "groups": ["PDGF (prolonged DGF)", "IGF"],
        },
        supplementary_files=["GSE293480_rawCounts.csv.gz"],
        download_method="supplementary",
        condition_column="title",  # PDGF{n} o IGF{n}
        notes=(
            "Matrice counts gia' disponibile come CSV. "
            "Molto rilevante per Aim 3 (kidney transplant NMP). "
            "Sample size piccolo (7), usare come supporto, non discovery."
        ),
    ),

    # -------------------------------------------------------
    # SUPPORT: Affymetrix microarray datasets
    # -------------------------------------------------------
    "GSE54888": DatasetConfig(
        accession="GSE54888",
        title="Transcriptome analysis in preimplantation kidney biopsies",
        platform=Platform.AFFYMETRIX_GENE10ST,
        gpl="GPL6244",
        n_samples=54,
        organism="Homo sapiens",
        role=Role.SUPPORT,
        description=(
            "54 biopsie pre-impianto da donatori deceduti. "
            "27 DGF + 27 non-DGF. 13 con DGF prolungato (>14gg)."
        ),
        conditions={
            "n_DGF": 27, "n_nonDGF": 27, "n_pDGF": 13,
        },
        supplementary_files=["GSE54888_RAW.tar"],
        download_method="geoparse",  # Processed data in Sample table
        condition_column="characteristics_ch1",
        notes=(
            "Dati processati gia' nella Sample table di GEO. "
            "CEL files disponibili per re-processing con affy/oligo. "
            "Biopsie PRE-impianto: cattura firme precoci pre-IRI."
        ),
    ),

    "GSE37838": DatasetConfig(
        accession="GSE37838",
        title="Molecular assessment of implantation biopsies",
        platform=Platform.AFFYMETRIX_U133PLUS2,
        gpl="GPL570",
        n_samples=78,
        organism="Homo sapiens",
        role=Role.SUPPORT,
        description=(
            "78 biopsie da 53 donatori deceduti, prese durante reperfusione. "
            "Creatinina >265 umol/L a day 7 o necessita' di dialisi = poor function."
        ),
        conditions={
            "n_donors": 53,
            "endpoint": "serum creatinine >265 umol/L at day 7 or dialysis",
        },
        supplementary_files=["GSE37838_RAW.tar"],
        download_method="geoparse",
        condition_column="characteristics_ch1",
        notes=(
            "Biopsie DURANTE reperfusione (=T3 di GSE43974). "
            "Dati processati nella Sample table. "
            "Affymetrix U133 Plus 2.0 — piattaforma molto comune."
        ),
    ),

    "GSE53769": DatasetConfig(
        accession="GSE53769",
        title="Molecular regulation of acute kidney injury (mRNA)",
        platform=Platform.AFFYMETRIX_GENE20ST,
        gpl="GPL16686",
        n_samples=36,
        organism="Homo sapiens",
        role=Role.SUPPORT,
        description=(
            "18 biopsie zero-hour + 18 post-trapianto da 18 allotrapianti. "
            "8 AKI (necrosi tubulare acuta) + 10 controlli (biopsie protocollari)."
        ),
        conditions={
            "n_AKI": 8, "n_control": 10,
            "biopsy_types": ["zero-hour", "post-transplant"],
        },
        supplementary_files=["GSE53769_RAW.tar"],
        download_method="geoparse",
        condition_column="characteristics_ch1",
        notes=(
            "Design paired: zero-hour vs post-Tx per ogni allotrapianto. "
            "AKI = necrosi tubulare acuta istologica senza rigetto."
        ),
    ),
}


# ===========================================================
# HELPER FUNCTIONS
# ===========================================================

def get_discovery_datasets():
    return {k: v for k, v in DATASETS.items() if v.role == Role.DISCOVERY}

def get_validation_datasets():
    return {k: v for k, v in DATASETS.items() if v.role == Role.VALIDATION}

def get_support_datasets():
    return {k: v for k, v in DATASETS.items() if v.role == Role.SUPPORT}

def get_rnaseq_datasets():
    return {k: v for k, v in DATASETS.items() if v.platform == Platform.RNASEQ}

def get_microarray_datasets():
    return {k: v for k, v in DATASETS.items() if v.platform != Platform.RNASEQ}

def get_all_human_iri():
    """Tutti i dataset umani di IR-AKI ordinati per ruolo."""
    order = {Role.DISCOVERY: 0, Role.VALIDATION: 1, Role.SUPPORT: 2}
    return dict(sorted(DATASETS.items(), key=lambda x: order[x[1].role]))


def print_dataset_summary():
    """Stampa un riassunto di tutti i dataset."""
    print(f"{'Accession':<12} {'Piattaforma':<20} {'Campioni':>8} {'Ruolo':<12} {'Titolo'}")
    print("-" * 90)
    for acc, ds in get_all_human_iri().items():
        print(f"{acc:<12} {ds.platform.value:<20} {ds.n_samples:>8} {ds.role.value:<12} {ds.title[:40]}")
    total = sum(ds.n_samples for ds in DATASETS.values())
    print(f"\n{'Totale campioni:':>42} {total}")


if __name__ == "__main__":
    print_dataset_summary()
