#import "lib.typ": *
#let s = sec("compare")

= Comparison of ARG methods
#s.facts.methods on the same samples, gene by gene, with allele numbers removed (KPC-2 and KPC-97 are both KPC). K-MARVEL reports calls; RGI bwt reports every gene with a mapped read, so an RGI gene counts when reads cover at least #s.facts.rgi_cov % of its best reference sequence, and "RGI below cutoff" marks K-MARVEL genes RGI saw with partial coverage, which is not absence. ARGs-OAP uses SARG, whose names differ from CARD's, and is not compared here (#figref(s, "arg_method_upset"), #figref(s, "arg_method_by_category"), #tabref(s, "overlap")).

#fig(s, "arg_method_upset")
#fig(s, "arg_method_by_category")
#datatable(s, "overlap")
#datatable(s, "comparison")
