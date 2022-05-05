rule deseq2:
    input:
        table="/home/grahman/projects/qadabra/workflow/data/table.biom",
        metadata="/home/grahman/projects/qadabra/workflow/data/metadata.tsv"
    output:
        "results/deseq2.tsv"
    conda:
        "../envs/qadabra-da-R.yml"
    script:
        "../scripts/deseq2.R"
