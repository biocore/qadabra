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

logger.info("Loading differentials...")
diffs_files = [
    pd.read_table(x, sep="\t", index_col=0) for x in snakemake.input
]
concat_diffs = pd.concat(diffs_files, axis=1)

for col in concat_diffs.columns:
    new_col_name = col + " differentials"
    concat_diffs.rename(columns={col: new_col_name}, inplace=True)

concat_diffs.to_csv(snakemake.output[0], sep="\t", index=True)
logger.info(f"Saved to {snakemake.output[0]}")
