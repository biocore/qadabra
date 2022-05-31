import pandas as pd


diffs = pd.read_table(snakemake.input[0], sep="\t", index_col=0)
tool = snakemake.wildcards["tool"]
pctile = int(snakemake.wildcards["pctile"])

sorted_diffs = diffs.sort_values(by=tool, ascending=False)

n = int(diffs.shape[0]*pctile/100)
top_feats = sorted_diffs.head(n).index.tolist()
bot_feats = sorted_diffs.tail(n).index.tolist()
top_feats = pd.DataFrame(
    dict(feature_id=top_feats, location="numerator")
)
bot_feats = pd.DataFrame(
    dict(feature_id=bot_feats, location="denominator")
)
df = pd.concat([top_feats, bot_feats])
df["pctile"] = pctile
df["num_feats"] = n

df.to_csv(snakemake.output[0], sep="\t", index=False)
