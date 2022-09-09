import logging

from joblib import dump
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (roc_curve, auc, precision_recall_curve,
                             average_precision_score)
from sklearn.model_selection import RepeatedStratifiedKFold


logger = logging.getLogger("qadabra")
logger.setLevel(logging.INFO)
fh = logging.FileHandler(snakemake.log[0], mode="w")
formatter = logging.Formatter(
    f"[%(asctime)s - {snakemake.rule}] :: %(message)s"
)
fh.setFormatter(formatter)
logger.addHandler(fh)

logging.captureWarnings(True)
logging.getLogger("py.warnings").addHandler(fh)

covariate = snakemake.params["factor_name"]
target = snakemake.params["target"]

logger.info(f"Covariate: {covariate}")
logger.info(f"Target: {target}")

metadata = pd.read_table(snakemake.input["metadata"], sep="\t",
                         index_col=0).dropna(subset=[covariate])
log_ratios = pd.read_table(snakemake.input["log_ratios"], sep="\t",
                           index_col=0).join(metadata, how="inner")
X = log_ratios["log_ratio"].values.reshape(-1, 1)
y = (log_ratios[covariate] == target).astype(int).values

mean_fpr = np.linspace(0, 1, 100)
mean_rec = np.linspace(0, 1, 100)

n_splits = snakemake.config["ml_params"]["n_splits"]
n_repeats = snakemake.config["ml_params"]["n_repeats"]
logger.info(
    f"Using repeated stratified k-fold CV with {n_splits} splits & "
    f"{n_repeats} repeats"
)
cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats,
                             random_state=1)

model = LogisticRegression(penalty="none")
folds = [(train, test) for train, test in cv.split(X, y)]
model_data = {
    "folds": {
        f"fold_{i+1}": {"train_idx": train, "test_idx": test}
        for i, (train, test) in enumerate(folds)
    }
}
model_data["fprs"] = mean_fpr
model_data["recalls"] = mean_rec

model_data["log_ratios"] = X
model_data["truth"] = y

# https://stackoverflow.com/a/67968945
tprs = []
precs = []
roc_aucs = []
pr_aucs = []
for i, (train, test) in enumerate(folds):
    prediction = model.fit(X[train], y[train]).predict_proba(X[test])

    fpr, tpr, _ = roc_curve(y[test], prediction[:, 1])
    tpr_interp = np.interp(mean_fpr, fpr, tpr)
    tprs.append(tpr_interp)
    roc_auc = auc(mean_fpr, tpr_interp)
    roc_aucs.append(roc_auc)

    prec, rec, _ = precision_recall_curve(y[test], prediction[:, 1])
    prec = prec[::-1]
    rec = rec[::-1]
    prec_interp = np.interp(mean_rec, rec, prec)
    precs.append(prec_interp)
    pr_auc = auc(mean_rec, prec_interp)
    pr_aucs.append(pr_auc)

    logger.info(
        f"Fold {i+1}/{len(folds)}: ROC AUC = {roc_auc:.2f}, "
        f"PR AUC = {pr_auc:.2f}"
    )

model_data["tprs"] = tprs
model_data["precisions"] = precs
model_data["roc_aucs"] = roc_aucs
model_data["pr_aucs"] = pr_aucs

mean_roc_auc = np.mean(roc_aucs)
mean_pr_auc = np.mean(pr_aucs)

logger.info(f"Mean ROC AUC = {mean_roc_auc:.2f}")
logger.info(f"Mean PR AUC = {mean_pr_auc:.2f}")

dump(model_data, snakemake.output[0], compress=True)
logger.info(f"Saving model to {snakemake.output[0]}")
