TMPDIR := $(shell mktemp -d)
TABLE_FILE := $(shell realpath qadabra/test_data/table.biom)
MD_FILE := $(shell realpath qadabra/test_data/metadata.tsv)

snaketest:
	@set -e;
	echo $(TMPDIR); \
	qadabra create-workflow --workflow-dest $(TMPDIR); \
	qadabra add-dataset --workflow-dest $(TMPDIR) --table $(TABLE_FILE) --metadata $(MD_FILE) --name "ampharos" --factor-name anemia --target-level anemic --reference-level normal --verbose; \
	cd $(TMPDIR); \
	snakemake --use-conda --cores 4
