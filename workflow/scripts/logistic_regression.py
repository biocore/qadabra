from joblib import dump
import numpy as np
import pandas as pd
from scipy import interp
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import RepeatedKFold

from utils import get_logger


logger = get_logger(snakemake.log[0], snakemake.rule)

covariate = snakemake.config["model"]["covariate"]
target = snakemake.config["model"]["target"]

logger.info(f"Covariate: {covariate}")
logger.info(f"Target: {target}")

metadata = pd.read_table(snakemake.input["metadata"], sep="\t",
                         index_col=0).dropna(subset=[covariate])
log_ratios = pd.read_table(snakemake.input["log_ratios"], sep="\t",
                           index_col=0).join(metadata, how="inner")
X = log_ratios["log_ratio"].values.reshape(-1, 1)
y = (log_ratios[covariate] == target).astype(int).values

mean_fpr = np.linspace(0, 1, 100)

n_splits = snakemake.config["ml_params"]["n_splits"]
n_repeats = snakemake.config["ml_params"]["n_repeats"]
logger.info(
    f"Using repeated k-fold CV with {n_splits} splits & {n_repeats} repeats"
)
cv = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=1)

model = LogisticRegression(penalty="none")
folds = [(train, test) for train, test in cv.split(X, y)]
model_data = {
    "folds": {
        f"fold_{i+1}": {"train_idx": train, "test_idx": test}
        for i, (train, test) in enumerate(folds)
    }
}
model_data["fprs"] = mean_fpr

model_data["log_ratios"] = X
model_data["truth"] = y

tprs = []
aucs = []
for i, (train, test) in enumerate(folds):
    prediction = model.fit(X[train], y[train]).predict_proba(X[test])
    fpr, tpr, _ = roc_curve(y[test], prediction[:, 1])
    tprs.append(np.interp(mean_fpr, fpr, tpr))
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)
    logger.info(f"Fold {i}/{len(folds)}: AUC = ({roc_auc:.2f})")

model_data["tprs"] = tprs
model_data["aucs"] = aucs

dump(model_data, snakemake.output[0], compress=True)
logger.info(f"Saving model to {snakemake.output[0]}")
