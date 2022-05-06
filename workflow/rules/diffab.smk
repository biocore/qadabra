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
        "../scripts/deseq2.R"


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
        "../scripts/ancombc.R"

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
            --summary-dir results/songbird
        """
