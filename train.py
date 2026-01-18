import argparse, json
from pathlib import Path
import pandas as pd
import joblib

from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score

from xgboost import XGBClassifier

from .plots import ensure_dir, plot_roc, plot_pr, plot_calibration
from .config import Paths

def eval_probs(y_true, y_prob):
    return {
        "roc_auc": float(roc_auc_score(y_true, y_prob)),
        "pr_auc": float(average_precision_score(y_true, y_prob)),
    }

def main():
    p = Paths()
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=p.ml_table)
    args = ap.parse_args()

    df = pd.read_parquet(args.data)
    CONS = ["grantham","aa_pos"]
    FULL = ["grantham","aa_pos","sift_score","polyphen_score"]

    for c in set(CONS + FULL):
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(df[c].median())

    y = df["label"].astype(int)
    train_df, test_df, y_train, y_test = train_test_split(
        df, y, test_size=0.2, random_state=42, stratify=y
    )

    ensure_dir(p.figures_dir)
    ensure_dir(p.models_dir)
    Path(p.results_dir).mkdir(parents=True, exist_ok=True)

    metrics = {"n_rows": int(len(df)), "pos_rate": float(y.mean()), "models": {}, "features": {"cons": CONS, "full": FULL}}

    # LogReg baseline (conservation only)
    lr_cons = Pipeline([("scaler", StandardScaler()), ("clf", LogisticRegression(max_iter=2000))])
    lr_cons.fit(train_df[CONS], y_train)
    pr1 = lr_cons.predict_proba(test_df[CONS])[:, 1]
    metrics["models"]["logreg_cons"] = eval_probs(y_test, pr1)
    plot_roc(y_test, pr1, f"{p.figures_dir}/roc_logreg_cons.png")
    plot_pr(y_test, pr1, f"{p.figures_dir}/pr_logreg_cons.png")
    plot_calibration(y_test, pr1, f"{p.figures_dir}/cal_logreg_cons.png")
    joblib.dump(lr_cons, f"{p.models_dir}/logreg_cons.joblib")

    # XGB full
    xgb = XGBClassifier(
        n_estimators=400, max_depth=4, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        eval_metric="logloss", random_state=42
    )
    xgb.fit(train_df[FULL], y_train)
    pr2 = xgb.predict_proba(test_df[FULL])[:, 1]
    metrics["models"]["xgb_full"] = eval_probs(y_test, pr2)
    plot_roc(y_test, pr2, f"{p.figures_dir}/roc_xgb_full.png")
    plot_pr(y_test, pr2, f"{p.figures_dir}/pr_xgb_full.png")
    plot_calibration(y_test, pr2, f"{p.figures_dir}/cal_xgb_full.png")
    joblib.dump(xgb, f"{p.models_dir}/xgb_full.joblib")

    # Save importance + predictions
    pd.DataFrame({"feature": FULL, "importance": xgb.feature_importances_})\
      .sort_values("importance", ascending=False)\
      .to_csv(f"{p.results_dir}/xgb_feature_importance.csv", index=False)

    pred = test_df[["chrom","pos","ref","alt","label"]].copy()
    pred["p_logreg_cons"] = pr1
    pred["p_xgb_full"] = pr2
    pred.to_parquet(f"{p.results_dir}/test_predictions.parquet", index=False)

    with open(f"{p.results_dir}/metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)
    print(json.dumps(metrics, indent=2))

if __name__ == "__main__":
    main()
