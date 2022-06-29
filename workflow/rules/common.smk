all_diff_files = expand(
    "results/{tool}/differentials.processed.tsv", tool=config["tools"]
)
all_diff_files.extend(["results/concatenated_differentials.tsv", "results/qurro"])

all_ml = expand(
    "results/{tool}/ml/regression/model_data.pctile_{pctile}.joblib",
    tool=config["tools"],
    pctile=config["log_ratio_feat_pcts"],
)

all_viz_files = expand("figures/{tool}_differentials.svg", tool=config["tools"])
all_viz_files.extend(expand("figures/{viz}.svg", viz=["spearman_heatmap"]))
all_viz_files.extend(
    expand(
        "figures/upset/upset.pctile_{pctile}.{location}.svg",
        pctile=config["log_ratio_feat_pcts"],
        location=["numerator", "denominator"],
    )
)
all_viz_files.extend(
    expand(
        "figures/roc/roc.pctile_{pctile}.svg",
        pctile=config["log_ratio_feat_pcts"],
    )
)
all_viz_files.extend(
    expand(
        "figures/{fig}.html",
        fig=["pca", "rank_comparison"]
    )
)

all_input = all_viz_files.copy()
all_input.append("results/qurro")

covariate = config["model"]["covariate"]
reference = config["model"]["reference"]
target = config["model"]["target"]


def build_songbird_formula(wildcards):
    return f"C({covariate}, Treatment('{reference}'))"


diffab_tool_columns = {
    "edger": "logFC",
    "deseq2": "log2FoldChange",
    "ancombc": f"{covariate}{target}",
    "aldex2": f"model.{covariate}{target} Estimate",
    "songbird": f"C({covariate}, Treatment('{reference}'))[T.{target}]",
}
