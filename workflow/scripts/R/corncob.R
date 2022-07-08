library(biomformat)
library(corncob)
library(phyloseq)

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

taxa <- phyloseq::otu_table(table, taxa_are_rows=T)
meta <- phyloseq::sample_data(metadata)
physeq <- phyloseq::phyloseq(taxa, meta)

design.formula <- paste0("~", covariate)
if (length(confounders) != 0) {
    confounders_form = paste(confounders, collapse=" + ")
    design.formula <- paste0(design.formula, " + ", confounders_form)
}
design.formula <- as.formula(design.formula)

fit <- corncob::differentialTest(
    formula=design.formula,
    formula_null=~1,
    phi.formula=design.formula,
    phi.formula_null=~1,
    test="Wald",
    data=physeq,
    boot=F,
    full_output=T
)
saveRDS(fit, snakemake@output[[2]])

taxa_names <- names(fit$p)
colname <- paste0("mu.", covariate, target)

all_models <- fit$all_models
coefs <- c()
for (mod in all_models) {
    coef <- mod$coefficients[c(colname), "Estimate"]
    coefs <- c(coefs, coef)
}
coefs <- as.data.frame(coefs)
row.names(coefs) <- taxa_names

write.table(coefs, snakemake@output[[1]], sep="\t")
