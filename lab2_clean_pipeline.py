
import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.metrics import (
    matthews_corrcoef, accuracy_score, precision_score, recall_score,
    confusion_matrix, roc_auc_score, average_precision_score,
    RocCurveDisplay, PrecisionRecallDisplay, make_scorer
)
import matplotlib.pyplot as plt

# ===== CONFIG =====
PROJECT_DIR = "."
PROC_DIR = f"{PROJECT_DIR}/data/processed"
RESULTS_DIR = f"{PROJECT_DIR}/results"
FIG_DIR = f"{RESULTS_DIR}/figures"
TABLE_DIR = f"{RESULTS_DIR}/tables"

POS_REP_FASTA = f"{PROC_DIR}/pos_rep.fasta"
NEG_A_REP_FASTA = f"{PROC_DIR}/negA_rep.fasta"
NEG_B_REP_FASTA = f"{PROC_DIR}/negB_rep.fasta"

TRAIN_CSV = f"{PROC_DIR}/train_df.csv"
TEST_CSV  = f"{PROC_DIR}/test_df.csv"

RANDOM_SEED = 42
N_FOLDS = 5
NTERM_FEATURE_LEN = 22

def ensure_dirs():
    for d in [PROC_DIR, RESULTS_DIR, FIG_DIR, TABLE_DIR]:
        Path(d).mkdir(parents=True, exist_ok=True)

def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "seq" not in df.columns and "sequence" in df.columns:
        df = df.rename(columns={"sequence":"seq"})
    if "accession" not in df.columns:
        for cand in ["Accession","Entry","entry","id","ID"]:
            if cand in df.columns:
                df = df.rename(columns={cand:"accession"})
                break
    if "label" not in df.columns:
        for cand in ["Label","y","Y","target"]:
            if cand in df.columns:
                df = df.rename(columns={cand:"label"})
                break

    if "subset" not in df.columns:
        if "Transmembrane" in df.columns:
            df["subset"] = "unknown"
            df.loc[df["label"]==1, "subset"] = "pos"
            df.loc[(df["label"]==0) & (df["Transmembrane"].notna()), "subset"] = "negB"
            df.loc[(df["label"]==0) & (df["Transmembrane"].isna()), "subset"] = "negA"
        else:
            df["subset"] = np.where(df["label"]==1, "pos", "neg")

    missing = [c for c in ["accession","seq","label"] if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns {missing}. Got columns={list(df.columns)}")
    return df

# --- Features ---
AA = list("ACDEFGHIKLMNPQRSTVWY")
AA_IDX = {a:i for i,a in enumerate(AA)}
KD = {
 "A": 1.8, "C": 2.5, "D":-3.5, "E":-3.5, "F": 2.8, "G":-0.4, "H":-3.2, "I": 4.5, "K":-3.9,
 "L": 3.8, "M": 1.9, "N":-3.5, "P":-1.6, "Q":-3.5, "R":-4.5, "S":-0.8, "T":-0.7, "V": 4.2, "W":-0.9, "Y":-1.3
}
def nterm(seq, L=NTERM_FEATURE_LEN):
    return seq[:L] if len(seq) >= L else seq + "X"*(L-len(seq))
def aac_features(seq):
    s = nterm(seq)
    vec = np.zeros(len(AA), dtype=float); valid = 0
    for ch in s:
        if ch in AA_IDX:
            vec[AA_IDX[ch]] += 1; valid += 1
    if valid > 0: vec /= valid
    return vec
def kd_mean(seq):
    s = nterm(seq)
    vals = [KD.get(ch, np.nan) for ch in s]
    vals = [v for v in vals if not np.isnan(v)]
    return float(np.mean(vals)) if vals else 0.0
def featurize(df):
    X_aac = np.vstack(df["seq"].apply(aac_features).values)
    X_kd = df["seq"].apply(kd_mean).values.reshape(-1,1)
    return np.hstack([X_aac, X_kd])

def train_svm_grid(X_train, y_train):
    pipe = Pipeline([("scaler", StandardScaler()),
                     ("svc", SVC(kernel="rbf"))])
    param_grid = {
        "svc__C": [0.1, 1, 10, 100],
        "svc__gamma": ["scale", "auto", 0.1, 0.01, 0.001]
    }
    cv = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_SEED)
    scorer = make_scorer(matthews_corrcoef)
    grid = GridSearchCV(pipe, param_grid=param_grid, scoring=scorer, cv=cv, n_jobs=-1, refit=True)
    grid.fit(X_train, y_train)
    return grid

def evaluate(y_true, y_pred, y_score):
    m = {
        "accuracy": accuracy_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0),
        "mcc": matthews_corrcoef(y_true, y_pred),
        "roc_auc": roc_auc_score(y_true, y_score),
        "pr_auc": average_precision_score(y_true, y_score),
        "cm": confusion_matrix(y_true, y_pred),
    }
    return m

def attach_predictions(df, y_pred, y_score):
    out = df.copy()
    out["y_pred"] = y_pred
    out["y_score"] = y_score
    out["is_fp"] = (out["label"]==0) & (out["y_pred"]==1)
    out["is_fn"] = (out["label"]==1) & (out["y_pred"]==0)
    return out

def main():
    ensure_dirs()

    # Load split created by you
    train_df = normalize_columns(pd.read_csv(TRAIN_CSV))
    test_df  = normalize_columns(pd.read_csv(TEST_CSV))

    X_train = featurize(train_df); y_train = train_df["label"].values
    X_test  = featurize(test_df);  y_test  = test_df["label"].values

    grid = train_svm_grid(X_train, y_train)
    best_model = grid.best_estimator_
    print("Best params:", grid.best_params_, "CV MCC:", grid.best_score_)

    y_pred = best_model.predict(X_test)
    y_score = best_model.decision_function(X_test)

    m = evaluate(y_test, y_pred, y_score)
    print(m)

    pd.DataFrame([{
        "model":"svm_rbf_grid",
        "accuracy": m["accuracy"],
        "precision": m["precision"],
        "recall": m["recall"],
        "mcc": m["mcc"],
        "roc_auc": m["roc_auc"],
        "pr_auc": m["pr_auc"],
        "best_params": str(grid.best_params_),
        "cv_mcc": grid.best_score_
    }]).to_csv(f"{TABLE_DIR}/svm_metrics.csv", index=False)

    RocCurveDisplay.from_predictions(y_test, y_score)
    plt.savefig(f"{FIG_DIR}/svm_roc.png", dpi=200, bbox_inches="tight"); plt.close()

    PrecisionRecallDisplay.from_predictions(y_test, y_score)
    plt.savefig(f"{FIG_DIR}/svm_pr.png", dpi=200, bbox_inches="tight"); plt.close()

    cm = m["cm"]
    plt.figure()
    plt.imshow(cm)
    plt.title("SVM Confusion Matrix")
    plt.xlabel("Predicted"); plt.ylabel("True")
    plt.xticks([0,1],["No SP","SP"]); plt.yticks([0,1],["No SP","SP"])
    for (i,j), v in np.ndenumerate(cm):
        plt.text(j, i, str(v), ha="center", va="center")
    plt.savefig(f"{FIG_DIR}/svm_cm.png", dpi=200, bbox_inches="tight"); plt.close()

    pred_df = attach_predictions(test_df, y_pred, y_score)
    pred_df[pred_df["is_fp"]][["accession","subset","y_score"]].to_csv(f"{TABLE_DIR}/false_positives.csv", index=False)
    pred_df[pred_df["is_fn"]][["accession","subset","y_score"]].to_csv(f"{TABLE_DIR}/false_negatives.csv", index=False)

    print("Outputs written to results/")

if __name__ == "__main__":
    main()
