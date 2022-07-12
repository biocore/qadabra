import re


stylesheet = config["stylesheet"]


rule plot_differentials:
    input:
        "results/tools/{tool}/differentials.processed.tsv",
    output:
        report(
            "figures/{tool}_differentials.svg",
            caption="../report/plot_differentials.rst",
            category="Differentials",
            subcategory="Rank Plots",
            labels={"tool": "{tool}"},
        ),
    log:
        "log/plot_differentials.{tool}.log",
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
            category="Differentials",
            subcategory="Comparison",
        ),
    log:
        "log/plot_rank_correlation.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_rank_correlations.py"


rule plot_upset:
    input:
        expand(
            "results/ml/{tool}/pctile_feats/pctile_{{pctile}}.tsv",
            tool=config["tools"],
        ),
    output:
        numerator=report(
            "figures/upset/upset.pctile_{pctile}.numerator.svg",
            caption="../report/plot_upset.rst",
            category="UpSet",
            labels={"percentile": "{pctile}", "type": "numerator"},
        ),
        denominator=report(
            "figures/upset/upset.pctile_{pctile}.denominator.svg",
            caption="../report/plot_upset.rst",
            category="UpSet",
            labels={"percentile": "{pctile}", "type": "denominator"},
        ),
    wildcard_constraints:
        tool="(?!pca_pc1)\w*",
    log:
        "log/upset.pctile_{pctile}.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_upset.py"


rule plot_roc:
    input:
        expand(
            "results/ml/{tool}/regression/model_data.pctile_{{pctile}}.joblib",
            tool=config["tools"] + ["pca_pc1"],
        ),
    output:
        report(
            "figures/roc/roc.pctile_{pctile}.svg",
            caption="../report/plot_roc.rst",
            category="Logistic Regression",
            subcategory="ROC",
            labels={"percentile": "{pctile}"},
        ),
    log:
        "log/plot_roc.{pctile}.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_roc.py"


rule plot_pr:
    input:
        expand(
            "results/ml/{tool}/regression/model_data.pctile_{{pctile}}.joblib",
            tool=config["tools"] + ["pca_pc1"],
        ),
    output:
        report(
            "figures/pr/pr.pctile_{pctile}.svg",
            caption="../report/plot_prc.rst",
            category="Logistic Regression",
            subcategory="PR",
            labels={"percentile": "{pctile}"},
        ),
    log:
        "log/plot_prc.{pctile}.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_prc.py"


rule plot_pca:
    input:
        features="results/pca/pca_features.tsv",
        tools="results/pca/pca_tools.tsv",
        prop_exp="results/pca/proportion_explained.tsv",
    output:
        report(
            "figures/pca.svg",
            caption="../report/plot_pca.rst",
            category="Differentials",
            subcategory="Comparison",
        ),
    log:
        "log/plot_pca.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_pca.py"


rule plot_rank_comparison:
    input:
        "results/concatenated_differentials.tsv",
    output:
        report(
            "figures/rank_comparisons.html",
            caption="../report/plot_rank_comparison.rst",
            category="Differentials",
            subcategory="Comparison",
        ),
    log:
        "log/plot_rank_comparison.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_rank_comparison.py"


rule qurro:
    input:
        ranks="results/concatenated_differentials.tsv",
        table=config["table"],
        metadata=config["metadata"],
    output:
        report(
            directory("results/qurro"),
            caption="../report/qurro.rst",
            htmlindex="index.html",
            category="Differentials",
            subcategory="Comparison",
        ),
    log:
        "log/qurro.log",
    conda:
        "../envs/qadabra-songbird.yaml"
    shell:
        """
        qurro \
            --ranks {input.ranks} \
            --table {input.table} \
            --sample-metadata {input.metadata} \
            --output-dir {output} > {log} 2>&1
        """


rule create_table:
    input:
        "results/concatenated_differentials.tsv",
    output:
        report(
            "results/differentials_table.html",
            caption="../report/create_table.rst",
            category="Differentials",
            subcategory="Comparison",
        ),
    log:
        "log/create_table.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/create_table.py"


if config["tree"]:

    rule empress:
        input:
            tree=config["tree"],
            diff="results/concatenated_differentials.tsv",
        output:
            report(
                directory("results/empress"),
                caption="../report/empress.rst",
                htmlindex="empress.html",
                category="Tree",
            ),
        log:
            "log/empress.log",
        conda:
            "../envs/qadabra-songbird.yaml"
        shell:
            """
            empress tree-plot \
                --tree {input.tree} \
                --feature-metadata {input.diff} \
                --output-dir {output} > {log} 2>&1
            """
