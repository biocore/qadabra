import pandas as pd

from utils import get_logger


logger = get_logger(snakemake.log[0], snakemake.rule)
logger.info("Loading differentials...")
diffs = pd.read_table(snakemake.input[0], sep="\t", index_col=0)
col = "log2FoldChange"
logger.info(f"Using {col}")

diffs = diffs[col]
diffs.name = "deseq2"
diffs.index.name = "feature_id"
diffs.to_csv(snakemake.output[0], sep="\t", index=True)
logger.info(f"Saved to {snakemake.output[0]}")
