import pandas as pd


diffs = pd.read_table(snakemake.input[0], sep="\t", index_col=0)

covariate = snakemake.config["model"]["covariate"]
target = snakemake.config["model"]["target"]
col = f"{covariate}{target}"

diffs = diffs[col]
diffs.name = "ancombc"
diffs.index.name = "feature_id"
diffs.to_csv(snakemake.output[0], sep="\t", index=True)
