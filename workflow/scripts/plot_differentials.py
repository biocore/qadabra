import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from utils import get_logger

logger = get_logger(snakemake.log[0], snakemake.rule)
plt.style.use(snakemake.params[0])

logger.info("Loading differentials...")
diffs = pd.read_table(snakemake.input[0], sep="\t", index_col=0).squeeze()
values = diffs.sort_values().values

logger.info("Plotting differentials...")
fig, ax = plt.subplots(1, 1)
plt.bar(
    np.arange(len(values)),
    values,
    width=1
)
tool_name = snakemake.wildcards["tool"]
ax.set_title(f"{tool_name} Differentials")
ax.set_xlabel("Rank")
ax.grid("both")
plt.savefig(snakemake.output[0])
logger.info(f"Saved to {snakemake.output[0]}")
