#import "lib.typ": *
#let s = sec("overview")
#let f = s.facts

= Overview
This summary covers #f.n_samples samples from #f.n_batches results director#if f.n_batches == 1 [y] else [ies] (#f.batches)#if f.at("median_clean", default: none) != none [, with a median of #f.median_clean clean read pairs]. #f.stages_ran of #f.stages_total pipeline stages ran (#tabref(s, "stages")).

#if s.lists.attention.len() > 0 [
  == Attention
  #for a in s.lists.attention [
    *#a.title*
    #bullets(a.items.map(i => i.text))
  ]
]

// stage, one column per results directory (sharing at most 48 %), what the summary takes; capped tables lose the last
#let sh = s.tables.stages.pdf.header
#let nb = sh.len() - 2
#datatable(s, "stages", widths: if sh.last() == "In this summary" {
  (24%,) + (calc.min(12%, 48% / nb),) * nb + (1fr,)
} else { auto })

== Samples and sequencing depth
#f.n_low of #f.n_samples samples have fewer than #f.min_pairs clean read pairs#if f.n_low > 0 [ (#f.low_samples)]#[;] a gene or taxon missing from a low-depth sample is unknown, not absent (#figref(s, "qc_read_fate"), #tabref(s, "samples")).

#fig(s, "qc_read_fate")
#if has(s, "figures", "nonpareil_curves") [
  Nonpareil estimates how much of each sample's community its reads cover, from the redundancy of its R1 reads (#figref(s, "nonpareil_curves")): #f.n_low_cov of #f.n_samples samples are below #f.min_cov#if f.n_low_cov > 0 [ (#f.low_cov_samples)]. The dashed model reaches the 95 % line at the read pairs a sample would need for 95 % coverage (samples_qc nonpareil_depth_95pct).

  #fig(s, "nonpareil_curves")
]
#datatable(s, "samples")
