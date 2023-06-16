stylesheet = config["stylesheet"]


rule plot_differentials:
    input:
        "results/{dataset}/tools/{tool}/differentials.processed.tsv",
    output:
        report(
            "figures/{dataset}/{tool}_differentials.html",
            caption="../report/plot_differentials.rst",
            category="Differentials",
            subcategory="Rank Plots",
            labels={"tool": "{tool}"},
        ),
    log:
        "log/{dataset}/plot_differentials.{tool}.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_differentials.py"


rule plot_diff_kendall_heatmap:
    input:
        "results/{dataset}/concatenated_differentials.tsv",
    output:
        report(
            "figures/{dataset}/kendall_diff_heatmap.svg",
            caption="../report/plot_diff_kendall_heatmap.rst",
            category="Differentials",
            subcategory="Comparison",
        ),
    log:
        "log/{dataset}/plot_diff_kendall_heatmap.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_diff_kendall_heatmap.py"


rule plot_pvalue_kendall_heatmap:
    input:
        "results/{dataset}/concatenated_pvalues.tsv",
    output:
        report(
            "figures/{dataset}/kendall_pvalue_heatmap.svg",
            caption="../report/plot_pvalue_kendall_heatmap.rst",
            category="P-values",
            subcategory="Comparison",
        ),
    log:
        "log/{dataset}/plot_pvalue_kendall_heatmap.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_pvalue_kendall_heatmap.py"


rule plot_upset:
    input:
        expand(
            "results/{{dataset}}/ml/{tool}/pctile_feats/pctile_{{pctile}}.tsv",
            tool=config["tools"],
        ),
    output:
        numerator=report(
            "figures/{dataset}/upset/upset.pctile_{pctile}.numerator.svg",
            caption="../report/plot_upset.rst",
            category="UpSet",
            labels={"percentile": "{pctile}", "type": "numerator"},
        ),
        denominator=report(
            "figures/{dataset}/upset/upset.pctile_{pctile}.denominator.svg",
            caption="../report/plot_upset.rst",
            category="UpSet",
            labels={"percentile": "{pctile}", "type": "denominator"},
        ),
    wildcard_constraints:
        tool="(?!pca_pc1)\w*",
    log:
        "log/{dataset}/upset.pctile_{pctile}.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_upset.py"


rule plot_roc:
    input:
        expand(
            "results/{{dataset}}/ml/{tool}/regression/model_data.pctile_{{pctile}}.joblib",
            tool=config["tools"] + ["pca_pc1"],
        ),
    output:
        report(
            "figures/{dataset}/roc/roc.pctile_{pctile}.svg",
            caption="../report/plot_roc.rst",
            category="Logistic Regression",
            subcategory="ROC",
            labels={"percentile": "{pctile}"},
        ),
    log:
        "log/{dataset}/plot_roc.{pctile}.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_roc.py"


rule plot_pr:
    input:
        expand(
            "results/{{dataset}}/ml/{tool}/regression/model_data.pctile_{{pctile}}.joblib",
            tool=config["tools"] + ["pca_pc1"],
        ),
    output:
        report(
            "figures/{dataset}/pr/pr.pctile_{pctile}.svg",
            caption="../report/plot_prc.rst",
            category="Logistic Regression",
            subcategory="PR",
            labels={"percentile": "{pctile}"},
        ),
    log:
        "log/{dataset}/plot_prc.{pctile}.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_prc.py"


rule plot_pca:
    input:
        features="results/{dataset}/pca/pca_features.tsv",
        tools="results/{dataset}/pca/pca_tools.tsv",
        prop_exp="results/{dataset}/pca/proportion_explained.tsv",
    output:
        report(
            "figures/{dataset}/pca.svg",
            caption="../report/plot_pca.rst",
            category="Differentials",
            subcategory="Comparison",
        ),
    log:
        "log/{dataset}/plot_pca.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_pca.py"


rule plot_rank_comparison:
    input:
        "results/{dataset}/concatenated_differentials.tsv",
    output:
        report(
            "figures/{dataset}/rank_comparisons.html",
            caption="../report/plot_rank_comparison.rst",
            category="Differentials",
            subcategory="Comparison",
        ),
    log:
        "log/{dataset}/plot_rank_comparison.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_rank_comparison.py"


rule plot_pvalue_pw_comparisons:
    input:
        "results/{dataset}/concatenated_pvalues.tsv",
    output:
        report(
            "figures/{dataset}/pvalue_pw_comparisons.html",
            caption="../report/plot_pvalue_pw_comparisons.rst",
            category="P-values",
            subcategory="Comparison",
        ),
    log:
        "log/{dataset}/plot_pvalue_pw_comparisons.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_pvalue_pw_comparisons.py"


rule plot_pvalue_volcanoes:
    input:
        "results/{dataset}/tools/{tool}/differentials.tsv",
    output:
        report(
            "figures/{dataset}/{tool}_pvalue_volcanoes.html",
            caption="../report/plot_pvalue_volcanoes.rst",
            category="P-values",
            subcategory="Volcano Plots",
            labels={"tool": "{tool}"},
        ),
    log:
        "log/{dataset}/plot_pvalue_volcanoes.{tool}.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot_pvalue_volcanoes.py"

        
rule qurro:
    input:
        unpack(lambda wc: get_dataset_cfg(wc, ["table", "metadata"])),
        ranks="results/{dataset}/concatenated_differentials.tsv",
    output:
        report(
            directory("results/{dataset}/qurro"),
            caption="../report/qurro.rst",
            htmlindex="index.html",
            category="Differentials",
            subcategory="Comparison",
        ),
    log:
        "log/{dataset}/qurro.log",
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


rule create_diff_table:
    input:
        "results/{dataset}/concatenated_differentials.tsv",
    output:
        report(
            "results/{dataset}/differentials_table.html",
            caption="../report/create_diff_table.rst",
            category="Differentials",
            subcategory="Comparison",
        ),
    log:
        "log/{dataset}/create_diff_table.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/create_diff_table.py"

rule create_pval_table:
    input:
        "results/{dataset}/concatenated_pvalues.tsv",
    output:
        report(
            "results/{dataset}/pvalues_table.html",
            caption="../report/create_pval_table.rst",
            category="P-values",
            subcategory="Comparison",
        ),
    log:
        "log/{dataset}/create_pval_table.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/create_pval_table.py"