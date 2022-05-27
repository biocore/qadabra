import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


plt.style.use(snakemake.params[0])

diffs_df = pd.read_table(snakemake.input[0], sep="\t", index_col=0)

fig, ax = plt.subplots(1, 1)
sns.heatmap(
    diffs_df.corr("spearman"),
    annot=True,
    ax=ax
)
ax.set_title("Spearman Rank Correlations")
plt.savefig(snakemake.output[0])
