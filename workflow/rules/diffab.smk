rule deseq2:
    input:
        table=config["table"],
        metadata=config["metadata"]
    output:
        "results/deseq2/differentials.tsv",
        "results/deseq2/results.rds"
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/deseq2.R"


rule ancombc:
    input:
        table=config["table"],
        metadata=config["metadata"]
    output:
        "results/ancombc/differentials.tsv",
        "results/ancombc/results.rds"
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/ancombc.R"


rule aldex2:
    input:
        table=config["table"],
        metadata=config["metadata"]
    output:
        "results/aldex2/differentials.tsv",
        "results/aldex2/results.rds"
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/aldex2.R"


rule songbird:
    input:
        table=config["table"],
        metadata=config["metadata"]
    output:
        "results/songbird/differentials.tsv"
    conda:
        "../envs/qadabra-songbird.yaml"
    shell:
        """
        songbird multinomial \
            --input-biom {input.table} \
            --metadata-file {input.metadata} \
            --formula "C({config[model][covariate]}, Treatment('{config[model][reference]}'))" \
            --epochs {config[songbird_params][epochs]} \
            --differential-prior {config[songbird_params][differential_prior]} \
            --summary-interval 1 \
            --min-feature-count 0 \
            --min-sample-count 0 \
            --summary-dir results/songbird
        """


rule process_differentials:
    input:
        "results/{tool}/differentials.tsv"
    output:
        "results/{tool}/differentials.processed.tsv"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/process_differentials/process_{wildcards.tool}.py"


rule combine_differentials:
    input:
        expand("results/{tool}/differentials.processed.tsv", tool=config["tools"])
    output:
        "results/concatenated_differentials.tsv"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/process_differentials/concatenate_differentials.py"
