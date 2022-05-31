import pandas as pd


diffs = pd.read_table(snakemake.input[0], sep="\t", index_col=0)
diffs = diffs["log2FoldChange"]
diffs.name = "deseq2"
diffs.index.name = "feature_id"
diffs.to_csv(snakemake.output[0], sep="\t", index=True)
