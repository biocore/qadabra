import numpy as np
import pandas as pd


datasets = pd.read_table("config/datasets.tsv", sep="\t", index_col=0)
names = datasets.index

def get_dataset_cfg(wildcards, keys):
    d = datasets.loc[wildcards.dataset, keys].to_dict()
    if d.get("confounders"):
        d["confounders"] = d["confounders"].split(";")
    return d

def get_songbird_formula(wildcards):
    d = datasets.loc[wildcards.dataset].to_dict()

    covariate = d["factor_name"]
    reference = d["reference_level"]
    formula = f"C({covariate}, Treatment('{reference}'))"
    if d.get("confounders"):
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
        "ancombc": f"{covariate}{target}",
        "aldex2": f"model.{covariate}{target} Estimate",
        "songbird": f"C({covariate}, Treatment('{reference}'))[T.{target}]",
        "maaslin2": "coef",
        "metagenomeseq": f"{covariate}{target}",
        "corncob": "coefs",
    }
    return columns[wildcards.tool]
