library(biomformat)
library(DESeq2)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

table <- biomformat::read_biom(snakemake@input[["table"]])
table <- as.matrix(biomformat::biom_data(table))

metadata <- read.table(snakemake@input[["metadata"]], sep="\t", header=T,
                       row.names=1)

covariate <- snakemake@config[["model"]][["covariate"]]
target <- snakemake@config[["model"]][["target"]]
reference <- snakemake@config[["model"]][["reference"]]
confounders <- snakemake@config[["model"]][["confounders"]]

samples <- colnames(table)
metadata <- subset(metadata, rownames(metadata) %in% samples)
metadata[[covariate]] <- as.factor(metadata[[covariate]])
metadata[[covariate]] <- relevel(metadata[[covariate]], reference)
sample_order <- row.names(metadata)
table <- table[, sample_order]

confounders <- paste(confounders, collapse=" + ")
design.formula <- as.formula(paste0("~", covariate, " + ", confounders))
dds <- DESeq2::DESeqDataSetFromMatrix(
    countData=table,
    colData=metadata,
    design=design.formula
)
dds.results <- DESeq2::DESeq(dds, sfType="poscounts")
saveRDS(dds.results, snakemake@output[[2]])

results <- DESeq2::results(
    dds.results,
    format="DataFrame",
    tidy=TRUE,
    cooksCutoff=FALSE,
    contrast=c(covariate, target, reference)
)
row.names(results) <- rownames(table)
write.table(results, file=snakemake@output[[1]], sep="\t")
