library(biomformat)
library(Biobase)
library(metagenomeSeq)


log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

print("Loading table...")
table <- biomformat::read_biom(snakemake@input[["table"]])
table <- as.matrix(biomformat::biom_data(table))
table <- as.data.frame(table)

print("Loading metadata...")
metadata <- read.table(snakemake@input[["metadata"]], sep="\t", header=T,
                       row.names=1)

covariate <- snakemake@params[[1]][["factor_name"]]
target <- snakemake@params[[1]][["target_level"]]
reference <- snakemake@params[[1]][["reference_level"]]
confounders <- snakemake@params[[1]][["confounders"]]

print("Harmonizing table and metadata samples...")
samples <- colnames(table)
metadata <- subset(metadata, rownames(metadata) %in% samples)
metadata[[covariate]] <- as.factor(metadata[[covariate]])
metadata[[covariate]] <- relevel(metadata[[covariate]], reference)
sample_order <- row.names(metadata)
table <- table[, sample_order]
# Append F_ to features to avoid R renaming
row.names(table) <- paste0("F_", row.names(table))

pheno <- Biobase::AnnotatedDataFrame(metadata)

experiment <- metagenomeSeq::newMRexperiment(
    counts=table,
    phenoData=pheno,
)

pctile <- metagenomeSeq::cumNormStat(experiment, pFlag=T)
experiment_norm <- metagenomeSeq::cumNorm(experiment, p=pctile)
pd <- Biobase::pData(experiment_norm)

print("Creating design formula...")
design.formula <- paste0("~", covariate)
if (length(confounders) != 0) {
    for (c in confounders) {
        metadata[[c]] <- as.factor(metadata[[c]])
    }
    confounders_form = paste(confounders, collapse=" + ")
    design.formula <- paste0(design.formula, " + ", confounders_form)
}
design.formula <- as.formula(design.formula)
print(design.formula)

mm <- model.matrix(design.formula, data=pd)
print("Running metagenomeSeq...")
fit <- metagenomeSeq::fitZig(experiment_norm, mm)

saveRDS(fit, snakemake@output[[2]])
print("Saved RDS!")

results <- metagenomeSeq::MRcoefs(
    obj=fit,
    number=dim(table)[1]
)
row.names(results) <- gsub("^F_", "", row.names(results))
write.table(results, file=snakemake@output[[1]], sep="\t")
print("Saved differentials!")
