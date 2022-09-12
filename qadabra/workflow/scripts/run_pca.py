import logging

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


logger = logging.getLogger("qadabra")
logger.setLevel(logging.INFO)
fh = logging.FileHandler(snakemake.log[0], mode="w")
formatter = logging.Formatter(
    f"[%(asctime)s - {snakemake.rule}] :: %(message)s"
)
fh.setFormatter(formatter)
logger.addHandler(fh)

logger.info("Loading differentials...")
diff_df = pd.read_table(snakemake.input[0], sep="\t", index_col=0)
rank_df = (
    diff_df.rank(ascending=False)
    .rename(columns={f"{x}": f"{x}_rank" for x in diff_df.columns})
    .join(diff_df)
)

tool_list = diff_df.columns.tolist()
tool_rank_list = [x + "_rank" for x in tool_list]

logger.info("Scaling values...")
scaled_diff_values = StandardScaler().fit_transform(diff_df.values)

pc_cols = [f"PC{x+1}" for x in range(len(tool_list))]
logger.info("Running PCA...")
pca = PCA(n_components=len(tool_list)).fit(scaled_diff_values)
feat_df = pd.DataFrame(
    pca.transform(scaled_diff_values),
    index=diff_df.index,
    columns=pc_cols
).join(rank_df)
max_vals = feat_df[pc_cols].max()

tool_df = pd.DataFrame(
    np.dot(pca.components_, np.diag(max_vals)),
    index=diff_df.columns,
    columns=pc_cols
)
tool_df.index.name = "tool"

feat_df.to_csv(snakemake.output["features"], sep="\t", index=True)
tool_df.to_csv(snakemake.output["tools"], sep="\t", index=True)

logger.info(f"Saved feature PCs to {snakemake.output['features']}")
logger.info(f"Saved tool PCs to {snakemake.output['tools']}")

prop_exp = pca.explained_variance_ratio_
prop_exp = pd.Series(prop_exp, index=pc_cols)
prop_exp.name = "proportion_var_explained"
prop_exp.to_csv(snakemake.output["prop_exp"], sep="\t", index=True)
logger.info(f"Saved explained variance to {snakemake.output['prop_exp']}")
