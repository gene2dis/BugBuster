#import "lib.typ": *
#let s = sec("contig_functional")
#let f = s.facts

= Functional profile (contig level)
eggNOG-mapper and dbCAN annotated the predicted genes of the assemblies; featureCounts counted them per sample (#f.n_samples samples; #f.ontologies) in #f.unit. TPM is relative to the reads on predicted genes, and a gene with several terms counts for each. Not comparable row by row with the read-level profile (#figref(s, "function_contig_cog_heatmap"), #tabref(s, "per_sample")).

#fig(s, "function_contig_cog_heatmap")
#datatable(s, "per_sample")
#datatable(s, "cog")
#datatable(s, "ko")
