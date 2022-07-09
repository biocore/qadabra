import logging

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
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

logger.info("Loading PCA results...")
feat_df = pd.read_table(snakemake.input["features"], sep="\t", index_col=0)
tool_df = pd.read_table(snakemake.input["tools"], sep="\t", index_col=0)
prop_exp = pd.read_table(snakemake.input["prop_exp"], sep="\t", index_col=0)
prop_exp = prop_exp.squeeze()

pc_cols = tool_df.columns.tolist()
tool_list = tool_df.index.tolist()

arrow_palette = dict(zip(
    tool_list,
    sns.color_palette("colorblind", len(tool_list)).as_hex()
))

fig, ax = plt.subplots(1, 1)

logger.info("Plotting scatterplot...")
sns.scatterplot(
    data=feat_df,
    x="PC1",
    y="PC2",
    color="lightgray",
    edgecolor="black",
    linewidth=0.1,
    ax=ax
)

xlim = ax.get_xlim()
ylim = ax.get_ylim()

xrange = xlim[1] - xlim[0]
yrange = ylim[1] - ylim[0]

logger.info("Plotting arrows...")
for i, row in tool_df.iterrows():
    tool_name = row.name
    color = arrow_palette[tool_name]
    ax.arrow(
        x=0, y=0,
        dx=row["PC1"], dy=row["PC2"],
        facecolor=color,
        linewidth=0.5,
        head_width=xrange*0.005*5,
        width=xrange*0.005,
        edgecolor="black",
    )


patches = []
for i, row in tool_df.iterrows():
    tool_name = row.name
    color = arrow_palette[tool_name]
    ax.arrow(
        x=0, y=0,
        dx=row["PC1"], dy=row["PC2"],
        facecolor=color,
        linewidth=0.5,
        head_width=xrange*0.005*5,
        width=xrange*0.005,
        edgecolor="black",
    )
    tool_patch = Patch(
        facecolor=color,
        edgecolor="black",
        label=tool_name,
        linewidth=0.5,
    )
    patches.append(tool_patch)

ax.legend(
    handles=patches,
    bbox_to_anchor=[1, 1],
    loc="upper left",
    frameon=False
)

prop_exp_labels = [
    f"{i} ({x*100:.2f}%)"
    for i, x in prop_exp.iteritems()
]
ax.set_xlabel(prop_exp_labels[0])
ax.set_ylabel(prop_exp_labels[1])

ax.set_title("PCA")
plt.savefig(snakemake.output[0])
logger.info(f"Saved to {snakemake.output[0]}")
