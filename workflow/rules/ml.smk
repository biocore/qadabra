rule pctile_feats:
    input:
        "results/concatenated_differentials.tsv"
    output:
        "results/{tool}/ml_feats/pctile_{pctile}.tsv",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/ml/get_pctile_feats.py"


rule log_ratios:
    input:
        table=config["table"],
        feat_lists=expand(
            "results/{{tool}}/ml_feats/pctile_{pctile}.tsv",
            pctile=config["log_ratio_feat_pcts"]
        )
    output:
        "results/{tool}/ml_feats/log_ratios.tsv"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/ml/get_log_ratios.py"
