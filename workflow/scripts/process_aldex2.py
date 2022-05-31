import pandas as pd


diffs = pd.read_table(snakemake.input[0], sep="\t", index_col=0)

covariate = snakemake.config["model"]["covariate"]
reference = snakemake.config["model"]["reference"]
target = snakemake.config["model"]["target"]
col = f"model.{covariate}{target} Estimate"

diffs = diffs[col]
diffs.name = "aldex2"
diffs.index.name = "feature_id"
diffs.to_csv(snakemake.output[0], sep="\t", index=True)
