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
    log:
        "log/pctile_feats.{tool}.pctile_{pctile}.log"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/get_pctile_feats.py"


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
    log:
        "log/log_ratios.{tool}.pctile_{pctile}.log"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/get_log_ratios.py"


rule logistic_regression:
    input:
        log_ratios="results/{tool}/ml/log_ratios/log_ratios.pctile_{pctile}.tsv",
        metadata=config["metadata"]
    output:
        "results/{tool}/ml/regression/model_data.pctile_{pctile}.joblib"
    log:
        "log/logistic_regression.{tool}.pctile_{pctile}.log"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/logistic_regression.py"