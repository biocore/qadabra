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
