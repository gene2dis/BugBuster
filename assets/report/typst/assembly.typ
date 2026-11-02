#import "lib.typ": *
#let s = sec("assembly")

= Assembly
#s.facts.n_assemblies MEGAHIT assembl#if s.facts.n_assemblies == 1 [y] else [ies] (#tabref(s, "stats")). Only contigs kept by the length filter go on to binning, contig taxonomy, ARGs and functions, and reads that are not on contigs are outside every contig-based table (#tabref(s, "read_use")).

#datatable(s, "stats")
#datatable(s, "read_use")
