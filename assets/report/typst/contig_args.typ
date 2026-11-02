#import "lib.typ": *
#let s = sec("contig_args")

= ARGs on contigs
DeepARG found #s.facts.n_genes ARG genes on the contigs of #s.facts.n_assemblies assembl#if s.facts.n_assemblies == 1 [y] else [ies]. Each is shown with its contig's taxon (Metabuli), bin and geNomad class: physical linkage, unlike the read-level candidate hosts. These are gene counts on contigs, not read abundances (#figref(s, "arg_contig_hosts"), #tabref(s, "hosts")).

#fig(s, "arg_contig_hosts")
#datatable(s, "hosts")
#datatable(s, "mobile")
#datatable(s, "long")
