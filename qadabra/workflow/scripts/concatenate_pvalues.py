import logging
import pandas as pd

logger = logging.getLogger("qadabra")
logger.setLevel(logging.INFO)
fh = logging.FileHandler(snakemake.log[0], mode="w")
formatter = logging.Formatter(
    f"[%(asctime)s - {snakemake.rule}] :: %(message)s"
)
fh.setFormatter(formatter)
logger.addHandler(fh)

logging.captureWarnings(True)
logging.getLogger("py.warnings").addHandler(fh)

logger.info("Loading p-values...")
pvalue_files = [
    pd.read_table(x, sep="\t", index_col=0) for x in snakemake.input
]
concat_pvalues = pd.concat(pvalue_files, axis=1)

concat_pvalues.to_csv(snakemake.output[0], sep="\t", index=True)
logger.info(f"Saved to {snakemake.output[0]}")
