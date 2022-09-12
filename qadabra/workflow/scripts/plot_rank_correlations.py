import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


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

plt.style.use(snakemake.config["stylesheet"])

logger.info("Loading differentials...")
diffs_df = pd.read_table(snakemake.input[0], sep="\t", index_col=0)

corr = diffs_df.corr("kendall")
mask = np.zeros_like(corr)
mask[np.triu_indices_from(mask)] = True

fig, ax = plt.subplots(1, 1)
sns.heatmap(
    corr,
    annot=True,
    mask=mask,
    square=True,
    ax=ax
)
ax.set_title("Kendall Rank Correlations")
plt.savefig(snakemake.output[0])
logger.info(f"Saved to {snakemake.output[0]}")
