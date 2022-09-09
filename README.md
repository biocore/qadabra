![Main CI](https://github.com/gibsramen/qadabra/actions/workflows/main.yml/badge.svg)

# Qadabra

**Q**uantitative **A**nalysis of **D**ifferential **Ab**undance **Ra**nks

Qadabra is a Snakemake workflow for comparing the results of differential abundance tools.
Importantly, Qadabra focuses on feature *ranks* rather than FDR corrected p-values.

## Installation

```
pip install qadabra
```

Qadabra requires the following dependencies:

* snakemake
* click
* biom-format
* pandas
* numpy
* cython
* iow

## Usage

### Creating the workflow structure

Qadabra can be used on multiple datasets at once.
First, we want to create the worfklow structure to perfrom differential abundance with all tools.

```
qadabra create-workflow --workflow-dest my_qadabra
```

This command will initialize the workflow but we still need to point to our dataset(s) of interest.

### Adding a dataset

We can add datasets one-by-one with the `add-dataset` command.

```
qadabra add-dataset \
    --workflow-dest my_qadabra \
    --table data/table.biom \
    --metadata data/metadata.tsv \
    --name my_dataset_1 \
    --factor-name case_control \
    --target-level case \
    --reference-level control \
    --verbose
```

Let's walkthrough the arguments provided here:

* `workflow-dest`: The location of the workflow that we created earlier
* `table`: Feature table (features by samples) in BIOM format
* `metadata`: Sample metadata in TSV format
* `name`: Name to give this dataset
* `factor-name`: Metadata column to use for differential abundance
* `target-level`: The value in the chosen factor to use as the target
* `reference-level`: The reference level to which we want to compare our target
* `verbose`: Flag to show all preprocessing performed by Qadabra

You can use `qadabra add-dataset --help` for more details.
To add another dataset, just run this command again with the new dataset information.

### Running the workflow

The previous commands will create a subdirectory, `my_qadabra` in which the workflow structure is contained.
Navigate into this directory; you should see two folders: `config` and `workflow`.
If you open the `config/config.yaml` file, you can see a number of options with which to run Qadabra.
You can modify these as you like.
For example, if you want to only run DESeq2, ANCOM-BC, and Songbird, you can delete the other entries in the `tools` heading.

From the command line, execute `snakemake --use-conda <other options>` to start the workflow.
Please read the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for how to run Snakemake best on your system.

When this process is completed, you should have directories `figures`, `results`, and `log`.
Each of these directories will have a separate folder for each dataset you added.

### Generating a report

You can also generate a report of the workflow with the following command:

```
snakemake --report report.zip
```

This will create a zipped directory containing the report.
Unzip this file and open the `report.html` file to view the report in your browser.

## Additional workflow options

### Worfklow subset

In some cases you may not want to run the full workflow and may only be interested in just running the different tools.
You can use `snakemake all_differentials --use-conda <other options>` to eschew the machine learning and visualization parts of the workflow.

### Phylogenetic visualization

Qadabra allows users to visualize the differentials on a phylogenetic tree using [EMPress](https://journals.asm.org/doi/10.1128/mSystems.01216-20).
With EMPress, you can annotate the tree with the differentials as barplots.
This can be useful for determining phylogenetic signal in differential abundance.

### Incorporating confounders

You can also specify additional confounders to incorporate into your DA model.
When adding a dataset, use `--confounder <column name>` to add a confounder into your model.
You can add multiple confounders by adding more `--confounder <column name>` arguments to `add-dataset`.

## Workflow Overview

Qadabra runs several differential abundance tools on the same dataset.
The features are ranked according to their association with the given metadata covariate.
The top and bottom features are then used to create log-ratios according to [Morton 2019](https://doi.org/10.1038/s41467-019-10656-5) and [Fedarko 2020](https://github.com/biocore/qurro).
These log-ratios are used as predictors in logistic regression models to predict the class given the log-ratio.

### Output

Qadabra generates many results files including many intermediate files that can be explored further.

#### Results

Each tool's output is stored in a separate subdirectory.
For the R tools, an RDS object with the tool's R data is saved.
The raw outputs are processed and concatenated into a file called `concatenated_differentials.tsv`.
A Qurro visualization of all the tool ranks is generated at `results/<dataset>/qurro/index.html`.
An interactive table with all the tool outputs is at `results/<dataset>/differentials_table.html`.

For each tool, the ranked features are used for machine learning models.
The `config.yaml` file enumerates the percentile of feats to use for log-ratios.
For example, at the 5% percentile, the top 5% of features and the bottom 5% of features associated with `covariate` are used to compute a log-ratio for each sample.
This log-ratio is used in repeated K-fold cross-validation to determine how well this log-ratio can predict class membership using logistic regression.
The `ml` subdirectory of each tool contains the features used, sample log-ratios, and compressed model objects.

#### Figures

The differential rank plots of each tool are plotted as `<tool_name>_differentials.svg`.
A heatmap of the pairwise Kendall rank correlation among all pairs of tools is available as well.
We also generated interactive plots to help compare the ranks of different features from the tools.
`figures/pca.svg` generates a PCA plot of all the features, showing the concordance and discordance of results as well as the contribution of the tools.
You can use the `figures/rank_comparisons.html` webpage to dynamically explore the relationship between pairs of tools.
The `upset` subdirectory contains [UpSet](https://doi.org/10.1109%2FTVCG.2014.2346248) plots comparing the features from each tool.
Finally, the `roc` and `pr` subdirectories contain ROC and PR (respectively) plots of all tools at each percentile of features.
