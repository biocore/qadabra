from collections import defaultdict
import logging

import joblib
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
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

tools = snakemake.config["tools"] + ["pca_pc1"]
palette = dict(zip(
    tools, sns.color_palette("colorblind", len(tools))
))

pctile = snakemake.wildcards["pctile"]

models = defaultdict(dict)
for tool, model_loc in zip(tools, snakemake.input):
    logger.info(f"Loading {model_loc}...")
    models[tool] = joblib.load(model_loc)

fig, ax = plt.subplots(1, 1)
ax.set_aspect("equal")

best_tool = None
best_auc = 0
leg_lines = []
for tool in tools:
    this_model_data = models[tool]
    mean_fprs = this_model_data["fprs"]
    mean_tprs = np.mean(this_model_data["tprs"], axis=0)

    roc_aucs = this_model_data["roc_aucs"]
    mean_auc = np.mean(roc_aucs)
    std_auc = np.std(roc_aucs)

    if mean_auc > best_auc:
        best_auc = mean_auc
        best_tool = tool

    color = palette[tool]
    line = Line2D([0], [0], color=color, lw=2,
                  label=f"{tool} ({mean_auc:.2f} $\pm$ {std_auc:.2f})")
    leg_lines.append(line)

    ax.plot(mean_fprs, mean_tprs, lw=2, color=color)
    logger.info(f"{tool} Mean AUC = {mean_auc:.2f} +- {std_auc:.2f}")

leg = ax.legend(handles=leg_lines, loc="lower right", frameon=False)
for text in leg.get_texts():
    if best_tool in text._text:
        text.set_weight("bold")

logger.info(f"Best model: {best_tool} ({best_auc:.2f})")

ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="black")
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title(f"{pctile}% Features")

plt.savefig(snakemake.output[0])
logger.info(f"Saved to {snakemake.output[0]}")
