import pandas as pd

from utils import get_logger


logger = get_logger(snakemake.log[0], snakemake.rule)

logger.info("Loading differentials...")
diffs_files = [
    pd.read_table(x, sep="\t", index_col=0) for x in snakemake.input
]
concat_diffs = pd.concat(diffs_files, axis=1)
concat_diffs.to_csv(snakemake.output[0], sep="\t", index=True)
logger.info(f"Saved to {snakemake.output[0]}")
