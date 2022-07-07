all_diff_files = expand(
    "results/tools/{tool}/differentials.processed.tsv", tool=config["tools"]
)
all_diff_files.extend(
    [
        "results/concatenated_differentials.tsv",
        "results/qurro",
        "results/differentials_table.html",
    ]
)

all_ml = expand(
    "results/ml/{tool}/regression/model_data.pctile_{pctile}.joblib",
    tool=config["tools"] + ["pca_pc1"],
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
all_viz_files.extend(expand("figures/{fig}.html", fig=["pca", "rank_comparisons"]))

all_input = all_viz_files.copy()
all_input.extend(["results/qurro", "results/differentials_table.html"])

covariate = config["model"]["covariate"]
reference = config["model"]["reference"]
target = config["model"]["target"]
confounders = config["model"]["confounders"]


songbird_formula = f"C({covariate}, Treatment('{reference}'))"
songbird_formula += " + ".join(confounders)


diffab_tool_columns = {
    "edger": "logFC",
    "deseq2": "log2FoldChange",
    "ancombc": f"{covariate}{target}",
    "aldex2": f"model.{covariate}{target} Estimate",
    "songbird": f"C({covariate}, Treatment('{reference}'))[T.{target}]",
}
