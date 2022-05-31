import biom
import numpy as np
import pandas as pd


table = biom.load_table(snakemake.input["table"])
table = table.to_dataframe(dense=True).T

def log_ratio(table, top_feats, bot_feats):
    num_sum = table.loc[:, top_feats].sum(axis=1)
    denom_sum = table.loc[:, bot_feats].sum(axis=1)
    lr_df = pd.concat([num_sum, denom_sum], axis=1)
    lr_df.columns = ["num", "denom"]
    lr_df = lr_df.dropna(how="all")
    lr_df = lr_df + 1
    lr_df["log_ratio"] = np.log(lr_df["num"]/lr_df["denom"]).to_frame()
    return lr_df


all_lr_dfs = []
for f in snakemake.input["feat_lists"]:
    feat_df = pd.read_table(f, sep="\t")
    top_feats = feat_df.query("location == 'numerator'")["feature_id"]
    bot_feats = feat_df.query("location == 'denominator'")["feature_id"]
    lr_df = log_ratio(table, top_feats, bot_feats).reset_index()
    lr_df["pctile"] = feat_df["pctile"].unique().item()
    lr_df["num_feats"] = feat_df["num_feats"].unique().item()
    all_lr_dfs.append(lr_df)

total_lr_df = pd.concat(all_lr_dfs)
total_lr_df["tool"] = snakemake.wildcards["tool"]
total_lr_df = total_lr_df.rename(columns={"index": "sample_name"})
total_lr_df.to_csv(snakemake.output[0], sep="\t", index=False)
