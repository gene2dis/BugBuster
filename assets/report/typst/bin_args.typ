#import "lib.typ": *
#let s = sec("bin_args")

= ARGs of MAGs
The #s.facts.n_genes ARG genes of the MAGs were clustered at 90-100 % identity; #s.facts.n_shared clusters at 90 % are found in more than one MAG (#tabref(s, "sharing")). A cluster in several MAGs is the same sequence in different genomes: consistent with transfer, not evidence of it alone.

#datatable(s, "sharing")
#datatable(s, "long")
