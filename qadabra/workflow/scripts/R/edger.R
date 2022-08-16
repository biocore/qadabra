library(biomformat)
library(edgeR)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

print("Loading table...")
table <- biomformat::read_biom(snakemake@input[["table"]])
table <- as.matrix(biomformat::biom_data(table))

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

print("Creating design formula...")
design.formula <- paste0("~", covariate)
if (length(confounders) != 0) {
    confounders_form = paste(confounders, collapse=" + ")
    design.formula <- paste0(design.formula, " + ", confounders_form)
}
design.formula <- as.formula(design.formula)
print(design.formula)
mm <- model.matrix(design.formula, metadata)

print("Running edgeR...")
d <- edgeR::DGEList(counts=table, group=metadata[[covariate]])
d <- edgeR::calcNormFactors(d)
d <- edgeR::estimateDisp(d, design=mm)
fit <- edgeR::glmQLFit(d, mm, robust=T)
saveRDS(fit, snakemake@output[[2]])
print("Saved RDS!")

results <- fit$coefficients
row.names(results) <- gsub("^F_", "", row.names(results))

write.table(results, file=snakemake@output[[1]], sep="\t")
print("Saved differentials!")
