#import "lib.typ": *
#let s = sec("functional")
#let f = s.facts

= Functional profile (read level)
#f.backend profiled #f.n_samples samples (#f.ontologies) in #f.unit. A gene can carry several terms, so the terms of an ontology do not add up to the annotated reads; every ontology is a table (`function_read_<ontology>.tsv`). This profile is reported apart from the contig-level one#refs(s, figs: ("function_cog_heatmap",), tabs: ("per_sample",)).

#fig(s, "function_cog_heatmap")
#datatable(s, "per_sample")
#datatable(s, "ko")
