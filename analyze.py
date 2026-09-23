import argparse
import numpy as np
import pandas as pd

SOURCE_URL = "https://raw.githubusercontent.com/ab604/lung-neoantigen-supplement/main/Supplementary-Data-S9-nsclc-tested-neoantigens-2024-05-11.csv"

def analyze(source: str):
    df = pd.read_csv(source)
    x = df[df["table_name"].str.contains("HLA-", na=False)].copy()
    for c in ["median_mt_score", "median_fold_change", "rank_percent"]:
        x[c] = pd.to_numeric(x[c], errors="coerce")
    x = x.dropna(subset=["median_mt_score", "median_fold_change", "rank_percent"])

    x["rank_binding"] = x.groupby("accel_id")["median_mt_score"].rank(method="min", ascending=True)
    x["rank_agretopicity"] = x.groupby("accel_id")["median_fold_change"].rank(method="min", ascending=False)
    x["rank_composite"] = x.groupby("accel_id")["rank_percent"].rank(method="min", ascending=True)

    pairs = [("binding", "agretopicity"), ("binding", "composite"), ("agretopicity", "composite")]
    rows = []
    for patient, g in x.groupby("accel_id"):
        top = {
            "binding": set(g.nsmallest(min(3, len(g)), "rank_binding").index),
            "agretopicity": set(g.nsmallest(min(3, len(g)), "rank_agretopicity").index),
            "composite": set(g.nsmallest(min(3, len(g)), "rank_composite").index),
        }
        for a, b in pairs:
            union = top[a] | top[b]
            inter = top[a] & top[b]
            rows.append({
                "patient": patient,
                "criterion_a": a,
                "criterion_b": b,
                "jaccard": len(inter) / len(union) if union else np.nan,
            })

    overlap = pd.DataFrame(rows)
    summary = overlap.groupby(["criterion_a", "criterion_b"])["jaccard"].mean().reset_index()
    return overlap, summary

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=SOURCE_URL, help="Local CSV path or URL to Supplementary Data S9")
    ap.add_argument("--out", default="top3_overlap_by_patient.csv")
    args = ap.parse_args()
    overlap, summary = analyze(args.csv)
    overlap.to_csv(args.out, index=False)
    print(summary.to_string(index=False))
