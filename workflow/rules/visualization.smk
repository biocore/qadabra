import re


stylesheet = "config/qadabra.mplstyle"


rule plot_differentials:
    input:
        "results/{tool}/differentials.processed.tsv"
    output:
        "figures/{tool}_differentials.pdf"
    params:
        stylesheet
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot/plot_differentials.py"


rule plot_rank_correlation:
    input:
        "results/concatenated_differentials.tsv"
    output:
        "figures/spearman_heatmap.pdf"
    params:
        stylesheet
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot/plot_rank_correlations.py"


rule interactive:
    input:
        "results/concatenated_differentials.tsv"
    output:
        "results/qadabra.html"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/interactive_app.py"


rule upset:
    input:
        expand(
            "results/{tool}/ml/pctile_feats/pctile_{{pctile}}.tsv",
            tool=config["tools"],
        )
    output:
        numerator="figures/upset/upset_pctile_{pctile}.numerator.pdf",
        denominator="figures/upset/upset_pctile_{pctile}.denominator.pdf"
    params:
        stylesheet
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot/plot_upset.py"


rule plot_roc:
    input:
        expand(
            "results/{tool}/ml/regression/model_data.pctile_{{pctile}}.joblib",
            tool=config["tools"]
        )
    output:
        "figures/roc/roc.pctile_{pctile}.pdf"
    params:
        stylesheet
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot/plot_roc.py"
