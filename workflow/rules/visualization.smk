import re


stylesheet = "config/qadabra.mplstyle"


rule plot_differentials:
    input:
        "results/{tool}/differentials.processed.tsv"
    output:
        report(
            "figures/{tool}_differentials.svg",
            caption="../report/plot_differentials.rst",
            category="Visualization",
            subcategory="Differentials",
            labels={"tool": "{tool}"}
        )
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
        report(
            "figures/spearman_heatmap.svg",
            caption="../report/plot_rank_correlation.rst",
            category="Visualization",
        )
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
        numerator=report(
            "figures/upset/upset_pctile_{pctile}.numerator.svg",
            caption="../report/plot_upset.rst",
            category="Visualization",
            subcategory="UpSet",
            labels={"percentile": "{pctile}", "type": "numerator"}
        ),
        denominator=report(
            "figures/upset/upset_pctile_{pctile}.denominator.svg",
            caption="../report/plot_upset.rst",
            category="Visualization",
            subcategory="UpSet",
            labels={"percentile": "{pctile}", "type": "denominator"}
        )
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
        report(
            "figures/roc/roc.pctile_{pctile}.svg",
            caption="../report/plot_roc.rst",
            category="Visualization",
            subcategory="ROC",
            labels={"percentile": "{pctile}"}
        )
    params:
        stylesheet
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot/plot_roc.py"
