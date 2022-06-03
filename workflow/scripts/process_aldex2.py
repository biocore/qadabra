import pandas as pd

from utils import get_logger


logger = get_logger(snakemake.log[0], snakemake.rule)
logger.info("Loading differentials...")
diffs = pd.read_table(snakemake.input[0], sep="\t", index_col=0)

covariate = snakemake.config["model"]["covariate"]
reference = snakemake.config["model"]["reference"]
target = snakemake.config["model"]["target"]
col = f"model.{covariate}{target} Estimate"
logger.info(f"Using {col}")

diffs = diffs[col]
diffs.name = "aldex2"
diffs.index.name = "feature_id"
diffs.to_csv(snakemake.output[0], sep="\t", index=True)
logger.info(f"Saved to {snakemake.output[0]}")
