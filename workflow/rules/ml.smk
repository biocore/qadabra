rule run_pca:
    input:
        "results/{dataset}/concatenated_differentials.tsv",
    output:
        features="results/{dataset}/pca/pca_features.tsv",
        tools="results/{dataset}/pca/pca_tools.tsv",
        prop_exp="results/{dataset}/pca/proportion_explained.tsv",
    log:
        "log/{dataset}/run_pca.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/run_pca.py"


rule tool_pctile_feats:
    input:
        "results/{dataset}/tools/{tool}/differentials.processed.tsv",
    output:
        "results/{dataset}/ml/{tool}/pctile_feats/pctile_{pctile}.tsv",
    log:
        "log/{dataset}/pctile_feats.{tool}.pctile_{pctile}.log",
    wildcard_constraints:
        tool="(?!pca_pc1)\w*",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/get_pctile_feats.py"


rule pca_pctile_feats:
    input:
        rules.run_pca.output.features,
    output:
        "results/{dataset}/ml/pca_pc1/pctile_feats/pctile_{pctile}.tsv",
    log:
        "log/{dataset}/pctile_feats.pca_pc1.pctile_{pctile}.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/get_pctile_feats.py"


rule log_ratios:
    input:
        table=lambda wc: get_dataset_cfg(wc, ["table"])["table"],
        feats="results/{dataset}/ml/{tool}/pctile_feats/pctile_{pctile}.tsv",
    output:
        "results/{dataset}/ml/{tool}/log_ratios/log_ratios.pctile_{pctile}.tsv",
    log:
        "log/{dataset}/log_ratios.{tool}.pctile_{pctile}.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/get_log_ratios.py"


rule logistic_regression:
    input:
        log_ratios="results/{dataset}/ml/{tool}/log_ratios/log_ratios.pctile_{pctile}.tsv",
        metadata=lambda wc: get_dataset_cfg(wc, ["metadata"])["metadata"]
    output:
        "results/{dataset}/ml/{tool}/regression/model_data.pctile_{pctile}.joblib",
    log:
        "log/{dataset}/logistic_regression.{tool}.pctile_{pctile}.log",
    params:
        # I have no idea why this way works but the other way doesn't
        factor_name=lambda wc: get_dataset_cfg(wc, ["factor_name"])["factor_name"],
        target=lambda wc: get_dataset_cfg(wc, ["target_level"])["target_level"],
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/logistic_regression.py"
