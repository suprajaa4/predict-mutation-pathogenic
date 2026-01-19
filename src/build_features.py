import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

from .config import Paths


# Grantham distance (full matrix)

_AA_ORDER = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
_GRANTHAM = [
    [0,112,111,126,195, 91,107, 60, 86, 94, 96,106, 84,113, 27, 99, 58,148,112, 64],
    [112,0, 86, 96,180, 43, 54,125, 29, 97,102, 26, 91, 97,103,110, 71,101, 77, 96],
    [111,86, 0, 23,139, 46, 42, 80, 68,149,153, 94,142,158, 91, 46, 65,174,143,133],
    [126,96,23, 0,154, 61, 45, 94, 81,168,172,101,160,177,108, 65, 85,181,160,152],
    [195,180,139,154,0,154,170,159,174,198,198,202,196,205,169,112,149,215,194,192],
    [91, 43,46, 61,154, 0, 29, 87, 24,109,113, 53,101,116, 76, 68, 42,130, 99, 96],
    [107,54,42, 45,170, 29, 0, 98, 40,134,138, 56,126,140, 93, 80, 65,152,122,121],
    [60,125,80, 94,159, 87, 98, 0, 98,135,138,127,127,153, 42, 56, 59,184,147,109],
    [86, 29,68, 81,174, 24, 40, 98, 0, 94, 99, 32, 87,100, 77, 89, 47,115, 83, 84],
    [94, 97,149,168,198,109,134,135, 94, 0,  5,102, 10, 21, 95,142, 89, 61, 33, 29],
    [96,102,153,172,198,113,138,138, 99, 5,  0,107, 15, 22, 98,145, 92, 61, 36, 32],
    [106,26,94,101,202, 53, 56,127, 32,102,107, 0, 95,102,103,121, 78,110, 85, 97],
    [84, 91,142,160,196,101,126,127, 87, 10, 15, 95, 0, 28, 87,135, 81, 67, 36, 21],
    [113,97,158,177,205,116,140,153,100, 21, 22,102, 28, 0,114,155,103, 40, 22, 50],
    [27,103,91,108,169, 76, 93, 42, 77, 95, 98,103, 87,114, 0, 74, 38,147,110, 68],
    [99,110,46, 65,112, 68, 80, 56, 89,142,145,121,135,155, 74, 0, 58,177,144,124],
    [58, 71,65, 85,149, 42, 65, 59, 47, 89, 92, 78, 81,103, 38, 58, 0,128, 92, 69],
    [148,101,174,181,215,130,152,184,115, 61, 61,110, 67, 40,147,177,128, 0, 37, 88],
    [112,77,143,160,194, 99,122,147, 83, 33, 36, 85, 36, 22,110,144, 92, 37, 0, 55],
    [64, 96,133,152,192, 96,121,109, 84, 29, 32, 97, 21, 50, 68,124, 69, 88, 55, 0],
]
_AA_IDX = {aa: i for i, aa in enumerate(_AA_ORDER)}

def grantham_distance(ref_aa, alt_aa) -> float:
    if ref_aa is None or alt_aa is None:
        return np.nan
    r = str(ref_aa).strip().upper()
    a = str(alt_aa).strip().upper()
    if len(r) != 1 or len(a) != 1:
        return np.nan
    if r not in _AA_IDX or a not in _AA_IDX:
        return np.nan
    return float(_GRANTHAM[_AA_IDX[r]][_AA_IDX[a]])



def normalize_chrom_for_vep(chrom: str) -> str:
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    if c == "M":
        c = "MT"
    return c

def to_vep_vcf_line(chrom: str, pos_1based: int, ref: str, alt: str) -> str:
    # VEP region endpoint accepts VCF-like: chrom pos id ref alt qual filter info
    return f"{chrom} {int(pos_1based)} . {ref} {alt} . . ."

def vep_url(assembly: str, pick: bool, sift: bool, polyphen: bool) -> str:
    base = "https://rest.ensembl.org"
    if assembly.lower() in {"grch37", "hg19"}:
        base = "https://grch37.rest.ensembl.org"

    params = []
    if pick:
        params.append("pick=1")
    if sift:
        params.append("sift=1")
    if polyphen:
        params.append("polyphen=1")

    q = ("?" + "&".join(params)) if params else ""
    return f"{base}/vep/homo_sapiens/region{q}"

def post_vep(url: str, variant_lines: list[str], timeout_s: int = 60) -> requests.Response:
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    payload = {"variants": variant_lines}
    return requests.post(url, headers=headers, data=json.dumps(payload), timeout=timeout_s)

def call_vep_batch(url: str, variant_lines: list[str], max_retries: int = 6) -> list[dict]:
    # retry 429; fallback for 400 handled outside
    for attempt in range(max_retries):
        r = post_vep(url, variant_lines)
        if r.status_code == 429:
            time.sleep(1.0 * (2 ** attempt))
            continue
        if r.status_code >= 400:
            # let caller decide how to handle (e.g., fallback per variant)
            r.raise_for_status()
        return r.json()
    raise RuntimeError("VEP failed after retries (rate limiting?)")

def call_vep_safe(url: str, variant_lines: list[str]) -> list[dict]:
    """
    If the batch fails with 400 (often due to one malformed record),
    fall back to per-variant calls and skip the bad ones.
    """
    try:
        return call_vep_batch(url, variant_lines)
    except requests.HTTPError as e:
        # If it's not 400, re-raise
        resp = getattr(e, "response", None)
        if resp is None or resp.status_code != 400:
            raise

        # Per-variant fallback
        out = []
        for v in variant_lines:
            try:
                out.extend(call_vep_batch(url, [v]))
            except requests.HTTPError as e2:
                resp2 = getattr(e2, "response", None)
                code = resp2.status_code if resp2 is not None else "?"
                msg = (resp2.text[:200] if resp2 is not None else str(e2))
                print(f"Skipping bad variant (HTTP {code}): {v} | {msg}")
                continue
        return out

def extract_missense_tc(rec: dict) -> dict | None:
    """
    Find a transcript consequence that is missense and includes amino_acids + protein_start.
    """
    tcs = rec.get("transcript_consequences") or []
    for tc in tcs:
        if "missense_variant" not in (tc.get("consequence_terms") or []):
            continue

        aa = tc.get("amino_acids")          # e.g., "K/R"
        aa_pos = tc.get("protein_start")    # e.g., 150
        if not aa or "/" not in aa or aa_pos is None:
            continue

        ref_aa, alt_aa = aa.split("/", 1)
        return {
            "ref_aa": ref_aa,
            "alt_aa": alt_aa,
            "aa_pos": int(aa_pos),
            "sift_pred": tc.get("sift_prediction"),
            "sift_score": tc.get("sift_score"),
            "polyphen_pred": tc.get("polyphen_prediction"),
            "polyphen_score": tc.get("polyphen_score"),
        }
    return None


def main():
    p = Paths()
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", default=p.processed_variants)
    ap.add_argument("--out", default="data/external/vep_features.parquet")
    ap.add_argument("--batch_size", type=int, default=200)
    ap.add_argument("--max_rows", type=int, default=None)
    ap.add_argument("--assembly", choices=["grch38", "grch37"], default="grch38")
    ap.add_argument("--no_pick", action="store_true", help="Disable VEP pick=1 (use all transcripts; slower)")
    ap.add_argument("--no_sift", action="store_true", help="Disable requesting SIFT")
    ap.add_argument("--no_polyphen", action="store_true", help="Disable requesting PolyPhen")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()


    df = pd.read_parquet(args.variants)

    need = {"chrom", "pos", "ref", "alt"}
    if not need.issubset(df.columns):
        raise ValueError(f"Expected columns {sorted(need)} in {args.variants}, got {list(df.columns)}")



    VALID = {"A", "C", "G", "T"}
    df["ref"] = df["ref"].astype(str).str.upper().str.strip()
    df["alt"] = df["alt"].astype(str).str.upper().str.strip()
    df = df[df["ref"].isin(VALID) & df["alt"].isin(VALID)].copy()

    if len(df) == 0:
        raise ValueError(
        "After filtering to A/C/G/T SNVs, no variants remain.\n"
        "This likely means your preprocess output has ref/alt like '-' '.' or multi-allelic values.\n"
        "Inspect a few rows of variants.parquet to confirm.")


    if args.max_rows:
        df = df.head(args.max_rows).copy()

    df["chrom_vep"] = df["chrom"].map(normalize_chrom_for_vep)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["chrom_vep", "pos", "ref", "alt"]).copy()

    url = vep_url(
        assembly=args.assembly,
        pick=(not args.no_pick),
        sift=(not args.no_sift),
        polyphen=(not args.no_polyphen),
    )

    vep_lines = [
        to_vep_vcf_line(c, int(pos), ref, alt)
        for c, pos, ref, alt in df[["chrom_vep", "pos", "ref", "alt"]].itertuples(index=False)
    ]

    out_rows = []
    kept_missense = 0
    skipped_nonmissense = 0

    for start in tqdm(range(0, len(vep_lines), args.batch_size)):
        batch = vep_lines[start:start + args.batch_size]
        recs = call_vep_safe(url, batch)

        if args.debug and start == 0 and recs:
            r0 = recs[0]
            print("\nDEBUG first VEP record keys:", list(r0.keys()))
            print("DEBUG input:", r0.get("input"))
            print("DEBUG most_severe_consequence:", r0.get("most_severe_consequence"))
            print("DEBUG allele_string:", r0.get("allele_string"))
            print("DEBUG n transcript_consequences:", len(r0.get("transcript_consequences") or []))

        for rec in recs:
            miss = extract_missense_tc(rec)
            if miss is None:
                skipped_nonmissense += 1
                continue

            allele = rec.get("allele_string") or ""
            ref2 = allele.split("/")[0] if "/" in allele else None
            alt2 = allele.split("/")[1] if "/" in allele else None

            g = grantham_distance(miss["ref_aa"], miss["alt_aa"])

            out_rows.append({
                "chrom": str(rec.get("seq_region_name")),
                "pos": int(rec.get("start")) if rec.get("start") is not None else np.nan,
                "ref": ref2,
                "alt": alt2,
                "ref_aa": miss["ref_aa"],
                "alt_aa": miss["alt_aa"],
                "aa_pos": miss["aa_pos"],
                "grantham": g,
                "sift_pred": miss.get("sift_pred"),
                "sift_score": miss.get("sift_score"),
                "polyphen_pred": miss.get("polyphen_pred"),
                "polyphen_score": miss.get("polyphen_score"),
                "hgvsg": rec.get("hgvsg"),
            })
            kept_missense += 1

    out_df = pd.DataFrame(out_rows)

    # Ensure columns exist even if empty
    expected_cols = [
        "chrom","pos","ref","alt","ref_aa","alt_aa","aa_pos","grantham",
        "sift_pred","sift_score","polyphen_pred","polyphen_score","hgvsg"
    ]
    for c in expected_cols:
        if c not in out_df.columns:
            out_df[c] = pd.NA

    for c in ["pos","aa_pos","grantham","sift_score","polyphen_score"]:
        out_df[c] = pd.to_numeric(out_df[c], errors="coerce")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(args.out, index=False)

    print(f"\nSaved: {args.out}")
    print(f"Input SNVs sent to VEP (after A/C/G/T filter): {len(df):,}")
    print(f"Missense rows kept: {kept_missense:,}")
    print(f"Non-missense skipped (or missing amino_acids): {skipped_nonmissense:,}")

    # If empty, suggest GRCh37 test
    if kept_missense == 0:
        print("\nNOTE: 0 missense variants found.")
        print("Try: --assembly grch37 (your coordinates may be closer to GRCh37), OR increase --max_rows and rerun.")
    else:
        print(out_df[["grantham","sift_score","polyphen_score"]].describe(include="all"))


if __name__ == "__main__":
    main()

