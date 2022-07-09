rule deseq2:
    input:
        table=config["table"],
        metadata=config["metadata"],
    output:
        "results/tools/deseq2/differentials.tsv",
        "results/tools/deseq2/results.rds",
    log:
        "log/deseq2.log",
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/deseq2.R"


rule ancombc:
    input:
        table=config["table"],
        metadata=config["metadata"],
    output:
        "results/tools/ancombc/differentials.tsv",
        "results/tools/ancombc/results.rds",
    log:
        "log/ancombc.log",
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/ancombc.R"


rule aldex2:
    input:
        table=config["table"],
        metadata=config["metadata"],
    output:
        "results/tools/aldex2/differentials.tsv",
        "results/tools/aldex2/results.rds",
    log:
        "log/aldex2.log",
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/aldex2.R"


rule edger:
    input:
        table=config["table"],
        metadata=config["metadata"],
    output:
        "results/tools/edger/differentials.tsv",
        "results/tools/edger/results.rds",
    log:
        "log/edger.log",
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/edger.R"


rule songbird:
    input:
        table=config["table"],
        metadata=config["metadata"],
    output:
        "results/tools/songbird/differentials.tsv",
    log:
        "log/songbird.log",
    params:
        epochs=config["songbird_params"]["epochs"],
        diff_prior=config["songbird_params"]["differential_prior"],
        formula=songbird_formula,
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
            --summary-dir results/tools/songbird > {log} 2>&1
        """


rule maaslin2:
    input:
        table=config["table"],
        metadata=config["metadata"],
    output:
        diff_file="results/tools/maaslin2/differentials.tsv",
        out_dir=directory("results/tools/maaslin2/output"),
    log:
        "log/maaslin2.log",
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/maaslin2.R"


rule metagenomeseq:
    input:
        table=config["table"],
        metadata=config["metadata"],
    output:
        "results/tools/metagenomeseq/differentials.tsv",
        "results/tools/metagenomeseq/results.rds",
    log:
        "log/metagenomeseq.log",
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/metagenomeseq.R"


rule corncob:
    input:
        table=config["table"],
        metadata=config["metadata"],
    output:
        "results/tools/corncob/differentials.tsv",
        "results/tools/corncob/results.rds",
    log:
        "log/corncob.log",
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/corncob.R"


rule process_differentials:
    input:
        "results/tools/{tool}/differentials.tsv",
    output:
        "results/tools/{tool}/differentials.processed.tsv",
    log:
        "log/process_differentials.{tool}.log",
    params:
        col=lambda wildcards: diffab_tool_columns[wildcards.tool],
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/process_differentials.py"


rule combine_differentials:
    input:
        expand("results/tools/{tool}/differentials.processed.tsv", tool=config["tools"]),
    output:
        "results/concatenated_differentials.tsv",
    log:
        "log/combine_differentials.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/concatenate_differentials.py"
