library(biomformat)
library(edgeR)

table <- biomformat::read_biom(snakemake@input[["table"]])
table <- as.matrix(biomformat::biom_data(table))

metadata <- read.table(snakemake@input[["metadata"]], sep="\t", header=T,
                       row.names=1)

covariate <- snakemake@config[["model"]][["covariate"]]
target <- snakemake@config[["model"]][["target"]]
reference <- snakemake@config[["model"]][["reference"]]

samples <- colnames(table)
metadata <- subset(metadata, rownames(metadata) %in% samples)
metadata[[covariate]] <- as.factor(metadata[[covariate]])
metadata[[covariate]] <- relevel(metadata[[covariate]], reference)
sample_order <- row.names(metadata)
table <- table[, sample_order]

d <- edgeR::DGEList(counts=table, group=metadata[[covariate]])
d <- edgeR::calcNormFactors(d)
d <- edgeR::estimateCommonDisp(d)
d <- edgeR::estimateTagwiseDisp(d)
et <- edgeR::exactTest(d, pair=c(reference, target))
saveRDS(et, snakemake@output[[2]])

write.table(et$table, file=snakemake@output[[1]], sep="\t")
