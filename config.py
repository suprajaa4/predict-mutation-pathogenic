from dataclasses import dataclass

@dataclass(frozen=True)
class Paths:
    # ---- Raw data ----
    raw_clinvar: str = "data/raw/clinvar_variant_summary.txt.gz"

    # ---- Processed core variants ----
    processed_variants: str = "data/processed/variants.parquet"

    # ---- VEP-derived features (Grantham, SIFT, PolyPhen) ----
    vep_features: str = "data/external/vep_features.parquet"

    # ---- Final ML table ----
    ml_table: str = "data/processed/ml_table.parquet"

    # ---- Outputs ----
    results_dir: str = "results"
    figures_dir: str = "results/figures"
    models_dir: str = "models"
