import numpy as np
import pandas as pd


datasets = pd.read_table("config/datasets.tsv", sep="\t", index_col=0)
names = datasets.index

def get_dataset_cfg(wildcards, keys):
    d = datasets.loc[wildcards.dataset, keys].to_dict()
    if "confounders" in keys:
        if not np.isnan(d["confounders"]):
            d["confounders"] = d["confounders"].split(";")
        else:
            d["confounders"] = []
    return d

def get_songbird_formula(wildcards):
    d = datasets.loc[wildcards.dataset].to_dict()

    covariate = d["factor_name"]
    reference = d["reference_level"]
    formula = f"C({covariate}, Treatment('{reference}'))"
    if not np.isnan(d["confounders"]):
        confounders = d["confounders"].split(";")
        formula = f"{formula} + {' + '.join(confounders)}"
    return formula

def get_diffab_tool_columns(wildcards):
    d = datasets.loc[wildcards.dataset].to_dict()
    covariate = d["factor_name"]
    target = d["target_level"]
    reference = d["reference_level"]

    columns = {
        "edger": f"{covariate}{target}",
        "deseq2": "log2FoldChange",
        "ancombc": "coefs",
        "aldex2": f"model.{covariate}{target} Estimate",
        "songbird": f"C({covariate}, Treatment('{reference}'))[T.{target}]",
        "maaslin2": "coef",
        "metagenomeseq": f"{covariate}{target}",
        "corncob": "coefs",
    }
    return columns[wildcards.tool]


def get_pvalue_tool_columns(wildcards):
    d = datasets.loc[wildcards.dataset].to_dict()
    covariate = d["factor_name"]
    target = d["target_level"]
    reference = d["reference_level"]

    columns = {
        "edger": "PValue",
        "deseq2": "pvalue",
        "ancombc": "pvals",
        "aldex2": f"model.{covariate}{target} Pr(>|t|)",
        "maaslin2": "pval",
        "metagenomeseq": "pvalues",
        "corncob": "fit.p",
    }
    return columns[wildcards.tool]
      
all_differentials = expand(
    "results/{dataset}/{out}",
    dataset=names,
    out=["concatenated_differentials.tsv", "qurro", "differentials_table.html"]
)

all_pvalues = expand(
    "results/{dataset}/{out}",
    dataset=names,
    out=["concatenated_pvalues.tsv", "pvalues_table.html"]
)

all_results = expand(
    "results/{dataset}/{out}",
    dataset=names,
    out=["qadabra_all_result.tsv"]
)

pvalue_volcanoes = expand(
    "figures/{dataset}/{tool}_pvalue_volcanoes.html",
    dataset=names,
    tool=config["ptools"]
)

all_ml = expand(
    "results/{dataset}/ml/{tool}/regression/model_data.pctile_{pctile}.joblib",
    dataset=names,
    tool=config["tools"] + ["pca_pc1"],
    pctile=config["log_ratio_feat_pcts"],
)

all_diff_viz = expand(
    "figures/{dataset}/{tool}_differentials.html",
    dataset=names,
    tool=config["tools"]
)

all_viz_files = expand(
    "figures/{dataset}/{viz}",
    dataset=names,
    viz=["kendall_diff_heatmap.svg", "kendall_pvalue_heatmap.svg", "differential_pw_comparisons.html", "pvalue_pw_comparisons.html", "pca.svg"]
)

all_viz_files.extend(expand(
    "figures/{dataset}/upset/upset.pctile_{pctile}.{location}.svg",
    dataset=names,
    pctile=config["log_ratio_feat_pcts"],
    location=["numerator", "denominator"],
))
# all_viz_files.extend(expand(
#     "figures/{dataset}/{curve}/{curve}.pctile_{pctile}.svg",
#     dataset=names,
#     pctile=config["log_ratio_feat_pcts"],
#     curve=["pr", "roc"],
# ))

all_input = all_differentials + all_pvalues + pvalue_volcanoes + all_viz_files + all_ml + all_diff_viz + all_results

for dataset in datasets.iterrows():
    if not pd.isna(dataset[1]['tree']):
        empress_output = expand(
            "results/{dataset}/{out}",
            dataset=names,
            out=["empress"]
            )
        all_input = all_input + empress_output