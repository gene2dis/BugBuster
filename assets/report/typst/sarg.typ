#import "lib.typ": *
#let s = sec("sarg")
#let f = s.facts

= ARGs-OAP (SARG)
ARGs-OAP found #f.n_types SARG types and #f.n_subtypes subtypes in #f.n_samples samples. Abundance is copies per cell: SARG copies (reads weighted by the share of the reference they cover, singletons included) over the ARGs-OAP cell estimate. It is not RPKM, and SARG types are not CARD drug classes: compare samples within this section only (#figref(s, "arg_sarg_type_heatmap"), #tabref(s, "types")).#if f.not_normalised != "" [ No 16S or cell estimate, so nothing normalised: #f.not_normalised.]

#fig(s, "arg_sarg_type_heatmap")
#datatable(s, "types")
