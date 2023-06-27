library(biomformat)
library(Biobase)
library(metagenomeSeq)

# Set logging information
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Load the input table
print("Loading table...")
table <- biomformat::read_biom(snakemake@input[["table"]])
table <- as.matrix(biomformat::biom_data(table))
table <- as.data.frame(table)

# Load the metadata
print("Loading metadata...")
metadata <- read.table(snakemake@input[["metadata"]], sep="\t", header=T,
                       row.names=1)

covariate <- snakemake@params[[1]][["factor_name"]]
target <- snakemake@params[[1]][["target_level"]]
reference <- snakemake@params[[1]][["reference_level"]]
confounders <- snakemake@params[[1]][["confounders"]]

# Harmonize table and metadata samples
print("Harmonizing table and metadata samples...")
samples <- colnames(table)
metadata <- subset(metadata, rownames(metadata) %in% samples)
metadata[[covariate]] <- as.factor(metadata[[covariate]])
metadata[[covariate]] <- relevel(metadata[[covariate]], reference)
sample_order <- row.names(metadata)
table <- table[, sample_order]
row.names(table) <- paste0("F_", row.names(table)) # Append F_ to features to avoid R renaming

# Create AnnotatedDataFrame from metadata
pheno <- Biobase::AnnotatedDataFrame(metadata)

# Create MRexperiment object
experiment <- metagenomeSeq::newMRexperiment(
    counts=table,
    phenoData=pheno,
)

# Perform cumulative sum scaling (CSS) normalization
pctile <- metagenomeSeq::cumNormStat(experiment, pFlag=T)
experiment_norm <- metagenomeSeq::cumNorm(experiment, p=pctile)
pd <- Biobase::pData(experiment_norm)

# Create design formula for model matrix
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

mm <- model.matrix(design.formula, data=pd)

# Run MetagenomeSeq
print("Running metagenomeSeq...")
fit <- metagenomeSeq::fitZig(experiment_norm, mm)
saveRDS(fit, snakemake@output[[2]])
print("Saved RDS!")

# Create results table and modify results names
results <- metagenomeSeq::MRcoefs(
    obj=fit,
    number=dim(table)[1],
    adjustMethod="BH",
)
row.names(results) <- gsub("^F_", "", row.names(results))

# Save results to output file
write.table(results, file=snakemake@output[[1]], sep="\t")
print("Saved differentials!")
