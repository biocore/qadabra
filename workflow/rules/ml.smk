rule pctile_feats:
    input:
        "results/concatenated_differentials.tsv"
    output:
        report(
            "results/{tool}/ml/pctile_feats/pctile_{pctile}.tsv",
            category="Machine Learning",
            subcategory="Feature Rankings",
            labels={"tool": "{tool}", "percentile": "{pctile}"}
        )
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/ml/get_pctile_feats.py"


rule log_ratios:
    input:
        table=config["table"],
        feats="results/{tool}/ml/pctile_feats/pctile_{pctile}.tsv",
    output:
        report(
            "results/{tool}/ml/log_ratios/log_ratios.pctile_{pctile}.tsv",
            category="Machine Learning",
            subcategory="Sample Log-Ratios",
            labels={"tool": "{tool}", "percentile": "{pctile}"}
        )
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/ml/get_log_ratios.py"


rule logistic_regression:
    input:
        log_ratios="results/{tool}/ml/log_ratios/log_ratios.pctile_{pctile}.tsv",
        metadata=config["metadata"]
    output:
        "results/{tool}/ml/regression/model_data.pctile_{pctile}.joblib"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/ml/logistic_regression.py"
