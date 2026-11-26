#import "lib.typ": *
#let s = sec("sarg")
#let f = s.facts

= ARGs-OAP (SARG)
ARGs-OAP found #f.n_types SARG types and #f.n_subtypes subtypes in #f.n_samples samples. Abundance is copies per cell: SARG copies (reads weighted by the share of the reference they cover, singletons included) over the ARGs-OAP cell estimate. It is not RPKM, and SARG types are not CARD drug classes: compare samples within this section only#refs(s, figs: ("arg_sarg_type_heatmap",), tabs: ("types",)).#if f.not_normalised != "" [ No 16S or cell estimate, so nothing normalised: #f.not_normalised.]

#if f.at("card", default: false) [

argNorm maps each hit's SARG sequence to a CARD ARO, giving the hits the CARD gene and reporting tier the K-MARVEL and RGI calls use. With at least #f.min_hits hits mapped by a Perfect, Strict or Manual mapping, ARGs-OAP detects #f.n_card_genes CARD genes#if f.tier1_genes != "" [, tier 1: #f.tier1_genes]; #f.n_loose subtype x ARO rows rest on a Loose mapping (often a homolog) and get no tier. ARGs-OAP counts any read with a hit, so check a tier 1 gene against the method comparison and its read hits#refs(s, tabs: ("card",)).
]

#fig(s, "arg_sarg_type_heatmap")
#datatable(s, "types")
#datatable(s, "card")
