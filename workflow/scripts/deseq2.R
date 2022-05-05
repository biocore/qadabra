library(biomformat)
library(DESeq2)

table <- biomformat::read_biom(snakemake@input[["table"]])
table <- as.matrix(biomformat::biom_data(table))

metadata <- read.table(snakemake@input[["metadata"]], sep="\t", header=T,
                       row.names=1)

covariate <- snakemake@config[["covariate"]]
target <- snakemake@config[["target"]]
reference <- snakemake@config[["reference"]]

samples <- colnames(table)
metadata <- subset(metadata, rownames(metadata) %in% samples)
metadata[[covariate]] <- as.factor(metadata[[covariate]])
metadata[[covariate]] <- relevel(metadata[[covariate]], reference)
sample_order <- row.names(metadata)
table <- table[, sample_order]

design.formula <- as.formula(paste0("~", snakemake@config[["covariate"]]))
dds <- DESeq2::DESeqDataSetFromMatrix(
    countData=table,
    colData=metadata,
    design=design.formula
)
dds_results <- DESeq2::DESeq(dds, sfType="poscounts")

results <- DESeq2::results(
    dds_results,
    format="DataFrame",
    tidy=TRUE,
    cooksCutoff=FALSE,
    contrast=c(covariate, target, reference)
)
row.names(results) <- rownames(table)
write.table(results, file=snakemake@output[[1]], sep="\t")
