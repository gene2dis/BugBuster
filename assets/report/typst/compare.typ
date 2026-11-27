#import "lib.typ": *
#let s = sec("compare")

= Comparison of ARG methods
#s.facts.methods on the same samples, gene by gene, with allele numbers removed (KPC-2 and KPC-97 are both KPC). K-MARVEL reports calls; RGI bwt reports every gene with a mapped read, so an RGI gene counts when reads cover at least #s.facts.rgi_cov % of its best reference sequence, and "RGI below cutoff" marks K-MARVEL genes RGI saw with partial coverage, which is not absence. #if s.facts.args_oap [ARGs-OAP is compared through argNorm's mapping of its SARG hits to CARD genes: a gene counts when at least #s.facts.min_hits of its reads map to it with a Perfect, Strict or Manual mapping (Loose mappings, often a homolog, never count). ARGs-OAP reports any read with a filtered hit, so it lists more genes than the CARD methods, many of them below their thresholds; this is detection only, its copies per cell are not comparable with read pairs] else [ARGs-OAP uses SARG, whose names differ from CARD's, and is not compared here]#refs(s, figs: ("arg_method_upset", "arg_method_by_category"), tabs: ("overlap",)).

#fig(s, "arg_method_upset")
#fig(s, "arg_method_by_category")
#datatable(s, "overlap")
#datatable(s, "comparison")
