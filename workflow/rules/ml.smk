rule run_pca:
    input:
        "results/concatenated_differentials.tsv",
    output:
        features="results/pca/pca_features.tsv",
        tools="results/pca/pca_tools.tsv",
        prop_exp="results/pca/proportion_explained.tsv"
    log:
        "log/run_pca.log"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/run_pca.py"


rule tool_pctile_feats:
    input:
        "results/concatenated_differentials.tsv",
    output:
        "results/ml/{tool}/pctile_feats/pctile_{pctile}.tsv"
    log:
        "log/pctile_feats.{tool}.pctile_{pctile}.log",
    wildcard_constraints:
        tool="(?!pca_pc1)\w+"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/get_pctile_feats.py"


rule pca_pctile_feats:
    input:
        rules.run_pca.output.features
    output:
        "results/ml/pca_pc1/pctile_feats/pctile_{pctile}.tsv",
    log:
        "log/pctile_feats.pca_pc1.pctile_{pctile}.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/get_pctile_feats.py"


rule log_ratios:
    input:
        table=config["table"],
        feats="results/ml/{tool}/pctile_feats/pctile_{pctile}.tsv",
    output:
        "results/ml/{tool}/log_ratios/log_ratios.pctile_{pctile}.tsv",
    log:
        "log/log_ratios.{tool}.pctile_{pctile}.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/get_log_ratios.py"


rule logistic_regression:
    input:
        log_ratios="results/ml/{tool}/log_ratios/log_ratios.pctile_{pctile}.tsv",
        metadata=config["metadata"],
    output:
        "results/ml/{tool}/regression/model_data.pctile_{pctile}.joblib",
    log:
        "log/logistic_regression.{tool}.pctile_{pctile}.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/logistic_regression.py"
