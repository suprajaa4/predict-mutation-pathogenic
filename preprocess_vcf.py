import gzip
import re
from pathlib import Path

import pandas as pd
from .config import Paths

POS = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic/Pathogenic"}
NEG = {"Benign", "Likely_benign", "Benign/Likely_benign", "Likely_benign/Benign"}

def parse_info(info: str) -> dict:
    d = {}
    for part in info.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        else:
            d[part] = True
    return d

def normalize_clnsig(raw: str) -> str | None:
    if raw is None:
        return None
    # CLNSIG often has comma-separated values for multiple submissions.
    # We'll treat any clear benign-ish as NEG and any clear pathogenic-ish as POS,
    # and drop conflicting/uncertain.
    s = raw

    # Drop clearly noisy categories
    bad = ["Conflicting_interpretations", "Uncertain_significance"]
    if any(b in s for b in bad):
        return None

    # If it contains both benign and pathogenic tokens → ambiguous
    has_pos = ("Pathogenic" in s) or ("Likely_pathogenic" in s)
    has_neg = ("Benign" in s) or ("Likely_benign" in s)
    if has_pos and has_neg:
        return None

    if has_pos:
        return "POS"
    if has_neg:
        return "NEG"
    return None

def main():
    paths = Paths()
    in_vcf = "data/raw/clinvar.grch37.vcf.gz"
    out_path = paths.processed_variants
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)

    rows = []
    with gzip.open(in_vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, pos, vid, ref, alt, qual, flt, info = line.rstrip("\n").split("\t")[:8]

            # SNVs only (single ALT, single nucleotide ref+alt)
            if "," in alt:
                continue
            if len(ref) != 1 or len(alt) != 1:
                continue

            info_d = parse_info(info)
            clnsig = info_d.get("CLNSIG")
            bucket = normalize_clnsig(clnsig)
            if bucket is None:
                continue

            label = 1 if bucket == "POS" else 0

            # VCF positions are 1-based already
            rows.append({
                "chrom": str(chrom),
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "label": label,
                "clnsig_raw": clnsig
            })

    df = pd.DataFrame(rows)
    df.to_parquet(out_path, index=False)

    print("Saved:", out_path)
    print("Rows:", len(df))
    print("Positive rate:", df["label"].mean())

if __name__ == "__main__":
    main()
