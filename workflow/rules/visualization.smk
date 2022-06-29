import re


stylesheet = "config/qadabra.mplstyle"


rule plot_differentials:
    input:
        "results/{tool}/differentials.processed.tsv",
    output:
        report(
            "figures/{tool}_differentials.svg",
            caption="../report/plot_differentials.rst",
            category="Visualization",
            subcategory="Differentials",
            labels={"tool": "{tool}"},
        ),
    log:
        "log/plot_differentials.{tool}.log",
    params:
        stylesheet,
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_differentials.py"


rule plot_rank_correlation:
    input:
        "results/concatenated_differentials.tsv",
    output:
        report(
            "figures/spearman_heatmap.svg",
            caption="../report/plot_rank_correlation.rst",
            category="Visualization",
        ),
    log:
        "log/plot_rank_correlation.log",
    params:
        stylesheet,
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_rank_correlations.py"


rule interactive:
    input:
        "results/concatenated_differentials.tsv",
    output:
        report("results/qadabra.html", category="Visualization"),
    log:
        "log/interactive.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/interactive_app.py"


rule upset:
    input:
        expand(
            "results/{tool}/ml/pctile_feats/pctile_{{pctile}}.tsv",
            tool=config["tools"],
        ),
    output:
        numerator=report(
            "figures/upset/upset.pctile_{pctile}.numerator.svg",
            caption="../report/plot_upset.rst",
            category="Visualization",
            subcategory="UpSet",
            labels={"percentile": "{pctile}", "type": "numerator"},
        ),
        denominator=report(
            "figures/upset/upset.pctile_{pctile}.denominator.svg",
            caption="../report/plot_upset.rst",
            category="Visualization",
            subcategory="UpSet",
            labels={"percentile": "{pctile}", "type": "denominator"},
        ),
    log:
        "log/upset.pctile_{pctile}.log",
    params:
        stylesheet,
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_upset.py"


rule plot_roc:
    input:
        expand(
            "results/{tool}/ml/regression/model_data.pctile_{{pctile}}.joblib",
            tool=config["tools"],
        ),
    output:
        report(
            "figures/roc/roc.pctile_{pctile}.svg",
            caption="../report/plot_roc.rst",
            category="Visualization",
            subcategory="ROC",
            labels={"percentile": "{pctile}"},
        ),
    log:
        "log/plot_roc.{pctile}.log",
    params:
        stylesheet,
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_roc.py"


rule plot_pca:
    input:
        features="results/pca/pca_features.tsv",
        tools="results/pca/pca_tools.tsv",
        prop_exp="results/pca/proportion_explained.tsv"
    output:
        report("figures/pca.html", category="Visualization")
    log:
        "log/plot_pca.log"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_pca.py"
