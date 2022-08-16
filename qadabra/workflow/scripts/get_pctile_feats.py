import logging

import pandas as pd


logger = logging.getLogger("qadabra")
logger.setLevel(logging.INFO)
fh = logging.FileHandler(snakemake.log[0], mode="w")
formatter = logging.Formatter(
    f"[%(asctime)s - {snakemake.rule}] :: %(message)s"
)
fh.setFormatter(formatter)
logger.addHandler(fh)

logger.info("Loading differentials...")
diffs = pd.read_table(snakemake.input[0], sep="\t", index_col=0)
if snakemake.wildcards.get("tool") is not None:
    tool = snakemake.wildcards["tool"]
else:
    tool = "PC1"
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
logger.info(f"Saved to {snakemake.output[0]}")
