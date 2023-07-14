![Main CI](https://github.com/gibsramen/qadabra/actions/workflows/main.yml/badge.svg)

# Qadabra: **Q**uantitative **A**nalysis of **D**ifferential **Ab**undance **Ra**nks

##### (Pronounced *ka-da-bra*)

Qadabra is a Snakemake workflow for running and comparing several differential abundance methods (tools) on the same microbiome dataset.

Importantly, Qadabra focuses on both FDR corrected p-values *and* [feature ranks](https://www.nature.com/articles/s41467-019-10656-5) and generates visualizations of differential abundance results.

[![Schematic](images/Qadabra_schematic.svg)](images/Qadabra_schematic.pdf)

## Installation

We reccommend installing [mamba](https://anaconda.org/conda-forge/mamba) to manage your Qadabra environment. Once mamba is installed, create and activate your Qadabra environment:
```
mamba create -n qadabra_env python=3.9
mamba activate qadabra_env
```
Install Qadabra and it's dependencies using [pip](https://pypi.org/project/pip/):
```
pip install qadabra snakemake click biom-format pandas numpy cython iow
```
## Usage

### 1. Creating the workflow structure

Qadabra can be used on multiple datasets at once.
First, we want to create the workflow structure to perfrom differential abundance with all methods:

```
qadabra create-workflow --workflow-dest my_qadabra
```

This command will initialize the workflow, but we still need to point to our dataset(s) of interest.

### 2. Adding a dataset

We can add datasets one-by-one with the `add-dataset` command:

```
qadabra add-dataset \
    --workflow-dest my_qadabra \
    --table data/table.biom \
    --metadata data/metadata.tsv \
    --tree data/my_tree.nwk \
    --name my_dataset_1 \
    --factor-name case_control \
    --target-level case \
    --reference-level control \
    --confounder confounding_variable(s) \
    --verbose
```

Let's walkthrough the arguments provided here, which represent the inputs to Qadabra:

* `workflow-dest`: The location of the workflow that we created earlier
* `table`: Feature table (features by samples) in [BIOM](https://biom-format.org/) format
* `metadata`: Sample metadata in TSV format
* `tree`: Phylogenetic tree in .nwk or other tree format (optional)
* `name`: Name to give this dataset
* `factor-name`: Metadata column to use for differential abundance
* `target-level`: The value in the chosen factor to use as the target
* `reference-level`: The reference level to which we want to compare our target
* `confounder`: Any confounding variable metadata columns (optional)
* `verbose`: Flag to show all preprocessing performed by Qadabra

Your dataset should now be added as a line in `my_qadabra/config/datasets.tsv`. 

You can use `qadabra add-dataset --help` for more details. 
To add another dataset, just run this command again with the new dataset information.

### 3. Running the workflow

The previous commands will create a subdirectory, `my_qadabra` in which the workflow structure is contained.
From the command line, execute the following to start the workflow:
```
snakemake --use-conda --cores <number of cores preferred> <other options>
```
Please read the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for how to run Snakemake best on your system.

When this process is completed, you should have directories `figures`, `results`, and `log`.
Each of these directories will have a separate folder for each dataset you added.

### 4. Generating a report

After Qadabra has finished running, you can generate a Snakemake report of the workflow with the following command:

```
snakemake --report report.zip
```

This will create a zipped directory containing the report.
Unzip this file and open the `report.html` file to view the report containing results and visualizations in your browser.

## Additional workflow options

### Workflow subsetting

In some cases you may not want to run the full workflow and may only be interested in just running certain methods.

If you navigate into your `my_qadabra` directory, you should see two folders: `config` and `workflow`. If you open the `config/config.yaml` file, you can see a number of options with which to run Qadabra. You can modify these as you like to eschew certain parts of the workflow.
For example, if you want to only run DESeq2, ANCOM-BC, and Songbird, you can delete the other entries in the `methods` heading.

### Incorporating confounders

You can also specify additional confounders to incorporate into your DA model.
When adding a dataset, use `--confounder <column name>` to add a confounder into your model.
You can add multiple confounders by adding more `--confounder <column name>` arguments to `add-dataset`.

### Phylogenetic visualization

Qadabra allows users to visualize the differentials and p-values on an interactive phylogenetic tree using [EMPress](https://journals.asm.org/doi/10.1128/mSystems.01216-20).
With EMPress, you can annotate the tree with the differentials as barplots.
This can be useful for determining phylogenetic signal in differential abundance. See the [EMPress GitHub](https://github.com/biocore/empress) page for more information and tutorial.

## Outputs overview

Qadabra generates many results files and intermediate files that can be explored further.

#### Results files

The differential abundance results from Qadabra are outputted in terms of FDR corrected p-values and feature ranks. 
These results can be found in the `results/<dataset_name>/` directory. Let's walkthrough the Qadabra results files: 
* `concatenated_differentials.tsv`: TSV table containing the differentials from each method.
* `concatenated_pvalues.tsv`: TSV table containing the FDR corrected p-values from each method.
* `differentials_table.html`: HTML table displaying concatenated_differentials.tsv.
* `pvalues_table.html`: HTML table displaying concatenated_pvalues.tsv.
* `qadabra_all_result.tsv`: TSV table containing differentials, FDR corrected p-values, and the number of methods passing significance threshold of 0.05 for each feature. (This table is used as the metadata for EMPress if a phylogenetic tree input is present.)

Each method's individual outputs are stored in a separate subdirectory under the `results/<dataset_name>/methods/<method>` subdirectories. 
* `differentials.tsv`: This file contains the differential abundance results as outputted by each individiual method.
* `differentials.processed.tsv`: This file extracts just the differentials column from differentials.tsv.
* `pvalues.processed.tsv`: This file extracts just the p-value column from differentials.tsv.
* `results.rds`: For the R methods (all except Songbird), an RDS object with the method's R data is saved.

A [Qurro](https://github.com/biocore/qurro) visualization of all the method ranks is generated at `results/<dataset_name>/qurro/index.html`.

For each method, the ranked features are used for machine learning models.
The `results/<dataset_name>/ml` subdirectory of each method contains the features used, sample log-ratios, and compressed model objects.

Results from the PCA analysis can be found under `results/<dataset_name>/pca`.

#### Figures

The generated Snakemake report contains the following folder structure:
* `Differentials`
    - `Comparison`
        - `differentials_table.html`: HTML table displaying concatenated_differentials.tsv.
        - `kendall_diff_heatmap.svg`: Heatmap showing degree of concordance of differentials between methods.
        - `pca.svg`: PCA plot showing method-specific effects on the ranking of features.
        - `qurro/index.html`: Interactively explore feature ranks with Qurro.
        - `differential_pw_comparisons.html`: Interactive plot displaying pairwise correlations of differential between any two methods.
    - `UpSet plots`: [UpSet](https://doi.org/10.1109%2FTVCG.2014.2346248) plots comparing the features from each method for top and bottom 20%, 15%, 10%, and 5% of features.
    - `Rank plots`: Differential rank plots of each method. Features with a positive log ratio are more associated with `target-level`. Featues with a negative log ratio are more associated with `reference-level`.
* `P-values`
    - `Comparison`
        -  `pvalues_table.html`: HTML table displaying concatenated_pvalues.tsv.
        - `kendall_pvalue_heatmap.svg`: Heatmap showing degree of concordance of p-values between methods.
        - `pvalue_pw_comparisons.html`: Interactive plot displaying pairwise correlations of p-values between any two p-value producing methods.
    - `Volcano plots`: [Volcano plots](https://en.wikipedia.org/wiki/Volcano_plot_(statistics)) for each p-value producing method.
* `EMPress plot`: An EMPress.html file to interactively explore differential abundance results with respect to phylogenetic relationships (if a tree was provided).


## Tutorial
Coming soon: An in-depth [tutorial](tutorial.md) on the full Qadabra workflow and subsequent analysis interpretations using a microbiome dataset.

## FAQs
Coming soon: An [FAQs](FAQs.md) page of commonly asked question on the statistics and code pertaining to Qadabra.

## Citation
The manuscript for Qadabra is currently in progress. Please cite this GitHub page if Qadabra is used for your analysis. This project is licensed under the MIT License. See the [license](LICENSE) file for details.

## Issues and contributing

If you encounter any problems with Qadabra, please open a New Issue in the GitHub Issues table. Contributions are welcome and greatly appreciated. If you have any improvements or bug fixes, please follow these steps:

1. Fork the repository.
2. Clone the repository to your local computer: `git clone <link-to-forked-repo>`
3. Create a new branch: `git checkout -b feature/your-feature`
4. Make your changes and commit them: `git commit -m 'Add your feature'`
5. Push the branch to your forked repository: `git push origin feature/your-feature`
6. Open a pull request detailing your changes.

Please ensure that your code adheres to the existing code style and that you include appropriate tests.