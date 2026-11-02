#import "lib.typ": *
#let s = sec("args")
#let f = s.facts

= #s.title
The ARG figures and matrices use #f.primary#if f.others != "" [; the calls of #f.others are in arg_long.tsv#if f.comparison [ and in the comparison of methods]]. RGI genes count when reads cover at least #f.rgi_cov % of their best reference sequence. #if f.n_tier1_samples != none [#f.n_tier1_samples samples carry a tier 1 (clinical priority) gene; #f.n_tier1_genes distinct tier 1 genes.]

#fig(s, "arg_priority_heatmap")
#fig(s, "arg_drug_class_heatmap")
#fig(s, "arg_vs_depth")
#if has(s, "tables", "tier1") [== Tier 1 calls]
#datatable(s, "tier1")
#datatable(s, "group")
== Mutation-dependent hits and tier 2
#datatable(s, "tier4")
#datatable(s, "tier2")
