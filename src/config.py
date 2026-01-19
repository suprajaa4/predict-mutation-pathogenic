from dataclasses import dataclass

@dataclass(frozen=True)
class Paths:

    raw_clinvar: str = "data/raw/clinvar_variant_summary.txt.gz"

    processed_variants: str = "data/processed/variants.parquet"

    vep_features: str = "data/external/vep_features.parquet"
    ml_table: str = "data/processed/ml_table.parquet"

    results_dir: str = "results"
    figures_dir: str = "results/figures"
    models_dir: str = "models"

