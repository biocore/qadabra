import logging

import matplotlib.pyplot as plt
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

plt.style.use(snakemake.params[0])

logger.info("Loading differentials...")
diffs_df = pd.read_table(snakemake.input[0], sep="\t", index_col=0)

fig, ax = plt.subplots(1, 1)
sns.heatmap(
    diffs_df.corr("spearman"),
    annot=True,
    ax=ax
)
ax.set_title("Spearman Rank Correlations")
plt.savefig(snakemake.output[0])
logger.info(f"Saved to {snakemake.output[0]}")
