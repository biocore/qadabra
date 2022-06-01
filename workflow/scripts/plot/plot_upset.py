import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import UpSet, from_contents


plt.style.use(snakemake.params[0])

tool_names = snakemake.config["tools"]
pctile = snakemake.wildcards["pctile"]
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

denom_plt = UpSet(from_contents(denom_list), subset_size="count",
                  show_counts=True).plot()
plt.title(f"Denominator - {pctile}%")
plt.savefig(snakemake.output["denominator"])
