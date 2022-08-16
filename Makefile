TMPDIR := $(shell mktemp -d)
TABLE_FILE := $(shell realpath qadabra/test_data/table.biom)
MD_FILE := $(shell realpath qadabra/test_data/metadata.tsv)

create_rulegraph:
	snakemake -f --rulegraph | dot -Tpng > imgs/rule_graph.png

snaketest:
	@cd $(TMPDIR); \
	qadabra create-workflow --workflow-dest . ;\
	qadabra add-dataset --table $(TABLE_FILE) --metadata $(MD_FILE) --name "ampharos" --factor-name anemia --target-level anemic --reference-level normal --verbose ; \
	snakemake --use-conda --cores 2
