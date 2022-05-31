import pandas as pd


diffs = pd.read_table(snakemake.input[0], sep="\t", index_col=0)
diffs = diffs["logFC"]
diffs.name = "edger"
diffs.index.name = "feature_id"
diffs.to_csv(snakemake.output[0], sep="\t", index=True)
