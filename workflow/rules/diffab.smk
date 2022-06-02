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


rule edger:
    input:
        table=config["table"],
        metadata=config["metadata"]
    output:
        "results/edger/differentials.tsv",
        "results/edger/results.rds"
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/edger.R"


rule songbird:
    input:
        table=config["table"],
        metadata=config["metadata"]
    output:
        "results/songbird/differentials.tsv",
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


rule process_deseq2:
    input:
        "results/deseq2/differentials.tsv"
    output:
        "results/deseq2/differentials.processed.tsv"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/process_differentials/process_deseq2.py"


rule process_ancombc:
    input:
        "results/ancombc/differentials.tsv"
    output:
        "results/ancombc/differentials.processed.tsv"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/process_differentials/process_ancombc.py"


rule process_aldex2:
    input:
        "results/aldex2/differentials.tsv"
    output:
        "results/aldex2/differentials.processed.tsv"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/process_differentials/process_aldex2.py"


rule process_edger:
    input:
        "results/edger/differentials.tsv"
    output:
        "results/edger/differentials.processed.tsv"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/process_differentials/process_edger.py"


rule process_songbird:
    input:
        "results/songbird/differentials.tsv"
    output:
        "results/songbird/differentials.processed.tsv"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/process_differentials/process_songbird.py"


rule combine_differentials:
    input:
        expand("results/{tool}/differentials.processed.tsv", tool=config["tools"])
    output:
        "results/concatenated_differentials.tsv",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/process_differentials/concatenate_differentials.py"


rule qurro:
    input:
        "results/concatenated_differentials.tsv"
    output:
        report(
            directory("results/qurro"),
            htmlindex="index.html",
            category="Differential Abundance"
        )
    conda:
        "../envs/qadabra-songbird.yaml"
    shell:
        """
        qurro \
            --ranks {input} \
            --table {config[table]} \
            --sample-metadata {config[metadata]} \
            --output-dir results/qurro
        """
