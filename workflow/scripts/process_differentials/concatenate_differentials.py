import pandas as pd


diffs_files = [
    pd.read_table(x, sep="\t", index_col=0) for x in snakemake.input
]
concat_diffs = pd.concat(diffs_files, axis=1)
concat_diffs.to_csv(snakemake.output[0], sep="\t", index=True)
