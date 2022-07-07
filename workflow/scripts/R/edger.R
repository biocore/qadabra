library(biomformat)
library(edgeR)

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
mm <- model.matrix(design.formula, metadata)

d <- edgeR::DGEList(counts=table, group=metadata[[covariate]])
d <- edgeR::calcNormFactors(d)
d <- edgeR::estimateDisp(d, design=mm)
fit <- edgeR::glmQLFit(d, mm, robust=T)
saveRDS(fit, snakemake@output[[2]])

write.table(fit$coefficients, file=snakemake@output[[1]], sep="\t")
