import pandas as pd
from pathlib import Path
from .config import Paths

def to_ucsc(c):
    c = str(c)
    if c.startswith("chr"):
        return c
    if c in {"MT","M"}:
        return "chrM"
    return "chr" + c

def main():
    p = Paths()
    variants = pd.read_parquet(p.processed_variants)
    vep = pd.read_parquet("data/external/vep_features.parquet")

    variants["chrom"] = variants["chrom"].map(to_ucsc)
    vep["chrom"] = vep["chrom"].map(to_ucsc)

    df = variants.merge(
        vep[["chrom","pos","ref","alt","ref_aa","alt_aa","aa_pos","grantham","sift_score","polyphen_score"]],
        on=["chrom","pos","ref","alt"],
        how="inner"
    )

    # Impute numeric
    for col in ["grantham","aa_pos","sift_score","polyphen_score"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
        df[col] = df[col].fillna(df[col].median())

    Path(p.ml_table).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(p.ml_table, index=False)
    print("Saved:", p.ml_table, "rows:", len(df))

if __name__ == "__main__":
    main()
