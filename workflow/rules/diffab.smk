import os

da_args = ["table", "metadata"]
da_params = ["factor_name", "target_level", "reference_level", "confounders"]


rule deseq2:
    input:
        unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
        "results/{dataset}/tools/deseq2/differentials.tsv",
        "results/{dataset}/tools/deseq2/results.rds",
    log:
        "log/{dataset}/deseq2.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params)
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/deseq2.R"


rule ancombc:
    input:
        unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
         "results/{dataset}/tools/ancombc/differentials.tsv",
         "results/{dataset}/tools/ancombc/results.rds",
    log:
        "log/{dataset}/ancombc.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params)
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/ancombc.R"


rule aldex2:
    input:
        unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
        "results/{dataset}/tools/aldex2/differentials.tsv",
        "results/{dataset}/tools/aldex2/results.rds",
    log:
        "log/{dataset}/aldex2.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params)
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/aldex2.R"


rule edger:
    input:
        unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
        "results/{dataset}/tools/edger/differentials.tsv",
        "results/{dataset}/tools/edger/results.rds",
    log:
        "log/{dataset}/edger.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params)
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/edger.R"


rule songbird:
    input:
       unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
        "results/{dataset}/tools/songbird/differentials.tsv",
    log:
        "log/{dataset}/songbird.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params),
        epochs=config["songbird_params"]["epochs"],
        diff_prior=config["songbird_params"]["differential_prior"],
        formula=get_songbird_formula,
        outdir=lambda wc, output: os.path.dirname(output[0]),
    conda:
        "../envs/qadabra-songbird.yaml"
    shell:
        """
        songbird multinomial \
            --input-biom {input.table} \
            --metadata-file {input.metadata} \
            --formula "{params.formula}" \
            --epochs {params.epochs} \
            --differential-prior {params.diff_prior} \
            --summary-interval 1 \
            --min-feature-count 0 \
            --min-sample-count 0 \
            --random-seed 1 \
            --summary-dir {params.outdir} > {log} 2>&1
        """


rule maaslin2:
    input:
       unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
        diff_file="results/{dataset}/tools/maaslin2/differentials.tsv",
        out_dir=directory("results/{dataset}/tools/maaslin2/output"),
    log:
        "log/{dataset}/maaslin2.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params)
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/maaslin2.R"


rule metagenomeseq:
    input:
       unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
        "results/{dataset}/tools/metagenomeseq/differentials.tsv",
        "results/{dataset}/tools/metagenomeseq/results.rds",
    log:
        "log/{dataset}/metagenomeseq.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params)
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/metagenomeseq.R"


rule corncob:
    input:
       unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
        "results/{dataset}/tools/corncob/differentials.tsv",
        "results/{dataset}/tools/corncob/results.rds",
    log:
        "log/{dataset}/corncob.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params)
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/corncob.R"


rule process_differentials:
    input:
        "results/{dataset}/tools/{tool}/differentials.tsv",
    output:
        "results/{dataset}/tools/{tool}/differentials.processed.tsv",
    log:
        "log/{dataset}/process_differentials.{tool}.log",
    params:
        col=get_diffab_tool_columns
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/process_differentials.py"


rule combine_differentials:
    input:
        expand(
            "results/{{dataset}}/tools/{tool}/differentials.processed.tsv",
            tool=config["tools"]
        ),
    output:
        "results/{dataset}/concatenated_differentials.tsv",
    log:
        "log/{dataset}/combine_differentials.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/concatenate_differentials.py"
