import logging

import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import UpSet, from_contents


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

tool_names = snakemake.config["tools"]
pctile = snakemake.wildcards["pctile"]
logger.info("Loading percentile features...")
feat_df_dict = {
    tool: pd.read_table(f, sep="\t", index_col=0)
    for tool, f in zip(tool_names, snakemake.input)
}

num_list = {
    tool: df.query("location == 'numerator'").index.tolist()
    for tool, df in feat_df_dict.items()
}

denom_list = {
    tool: df.query("location == 'denominator'").index.tolist()
    for tool, df in feat_df_dict.items()
}

num_plt = UpSet(from_contents(num_list), subset_size="count",
                show_counts=True).plot()
plt.title(f"Numerator - {pctile}%")
plt.savefig(snakemake.output["numerator"])
logger.info(f"Saved to {snakemake.output['numerator']}")

denom_plt = UpSet(from_contents(denom_list), subset_size="count",
                  show_counts=True).plot()
plt.title(f"Denominator - {pctile}%")
plt.savefig(snakemake.output["denominator"])
logger.info(f"Saved to {snakemake.output['denominator']}")
