import pandas as pd

from utils import get_logger


logger = get_logger(snakemake.log[0], snakemake.rule)
logger.info(f"Loading {snakemake.wildcards['tool']} differentials...")
col = snakemake.params["col"]
logger.info(f"Using '{col}' as column")
diffs = pd.read_table(snakemake.input[0], sep="\t", index_col=0)[col]
diffs.name = snakemake.wildcards["tool"]
diffs.index.name = "feature_id"
diffs.to_csv(snakemake.output[0], sep="\t", index=True)
logger.info(f"Saved to {snakemake.output[0]}")
