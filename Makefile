create_rulegraph:
	snakemake -f --rulegraph | dot -Tpng > imgs/rule_graph.png
