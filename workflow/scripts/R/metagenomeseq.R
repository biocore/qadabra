library(biomformat)
library(Biobase)
library(metagenomeSeq)


log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

table <- biomformat::read_biom(snakemake@input[["table"]])
table <- as.matrix(biomformat::biom_data(table))
table <- as.data.frame(table)

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

pheno <- Biobase::AnnotatedDataFrame(metadata)

experiment <- metagenomeSeq::newMRexperiment(
    counts=table,
    phenoData=pheno,
)

pctile <- metagenomeSeq::cumNormStat(experiment, pFlag=T)
experiment_norm <- metagenomeSeq::cumNorm(experiment, p=pctile)
pd <- Biobase::pData(experiment_norm)

design.formula <- paste0("~", covariate)
if (length(confounders) != 0) {
    for (c in confounders) {
        metadata[[c]] <- as.factor(metadata[[c]])
    }
    confounders_form = paste(confounders, collapse=" + ")
    design.formula <- paste0(design.formula, " + ", confounders_form)
}
design.formula <- as.formula(design.formula)

mm <- model.matrix(design.formula, data=pd)
fit <- metagenomeSeq::fitZig(experiment_norm, mm)
saveRDS(fit, snakemake@output[[2]])

results <- metagenomeSeq::MRcoefs(
    obj=fit,
    number=dim(table)[1]
)
write.table(results, file=snakemake@output[[1]], sep="\t")
