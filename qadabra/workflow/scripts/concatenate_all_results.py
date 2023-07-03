import logging
import pandas as pd

logger = logging.getLogger("qadabra")
logger.setLevel(logging.INFO)
fh = logging.FileHandler(snakemake.log[0], mode="w")
formatter = logging.Formatter(
    f"[%(asctime)s - {snakemake.rule}] :: %(message)s"
)
fh.setFormatter(formatter)
logger.addHandler(fh)

logging.captureWarnings(True)
logging.getLogger("py.warnings").addHandler(fh)

logger.info("Loading differentials...")
concatenated_diff_file = pd.read_csv(snakemake.input["concatenated_differentials"], sep="\t", index_col=0)

for col in concatenated_diff_file.columns:
    new_col_name = col + " differentials"
    concatenated_diff_file.rename(columns={col: new_col_name}, inplace=True)
    
logger.info("Loading p-values...")
concatenated_pvalue_file = pd.read_csv(snakemake.input["concatenated_pvalues"], sep="\t", index_col=0)

for col in concatenated_pvalue_file.columns:
    new_col_name = col + " significant features"
    concatenated_pvalue_file.rename(columns={col: new_col_name}, inplace=True)


qadabra_all_result = pd.concat([concatenated_diff_file, concatenated_pvalue_file], axis=1)

for index, row in qadabra_all_result.iterrows():
    count = 0
    for col in qadabra_all_result.columns:
        if col.endswith(" significant features"):
            if row[col] < 0.05:
                count += 1
                qadabra_all_result.at[index, col] = "p<0.05"
            else:
                qadabra_all_result.at[index, col] = "ns"
    qadabra_all_result.at[index, "# of tools p > 0.05"] = count


qadabra_all_result.to_csv(snakemake.output[0], sep="\t", index=True)

