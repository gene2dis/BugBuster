#import "lib.typ": *
#let s = sec("mags")
#let f = s.facts

= MAGs
#f.n_mags refined MAGs: #f.high high, #f.medium medium and #f.low low quality by MIMAG (#figref(s, "mag_quality"), #tabref(s, "mags")). The High / Mid / Low classes of the run's binning plots are another, older scheme. #figref(s, "mag_abundance_heatmap") shows the share of each sample's sequenced bases on each MAG.

#if "gunc_fail" in f [
  GUNC flags #f.gunc_fail of #f.n_mags MAGs as likely chimeras of distant lineages (`gunc_pass` false#if f.gunc_fail > 0 [, circled in #figref(s, "mag_quality")]); #f.gunc_pass pass#if f.gunc_unscored > 0 [, and #f.gunc_unscored have no GUNC result (not scored, or from a batch without GUNC)]. A fail does not change the MIMAG tier: check such a MAG before using it.
]

#if "multisample" in f [
  MetaBAT2 binned each assembly with one depth column per sample of its group (multi-sample binning: every sample aligned to the concatenation of the group's assemblies), so its bins differ from those of a run without the option; the MAG depth (#tabref(s, "abundance")) is still each sample's own alignment.
]

#fig(s, "mag_quality")
#fig(s, "mag_abundance_heatmap")
#datatable(s, "mags")
#datatable(s, "abundance")
#if "n_clusters" in f [
  == Run-level clusters
  #f.n_clustered MAGs fall into #f.n_clusters clusters at #f.derep_ani % ANI (galah; the representative of each has the highest completeness − 5 × contamination). #tabref(s, "cluster_abundance") gives each representative's share of every sample's reads (CoverM). The clusters cover this run's samples only.
] else if "derep_empty" in f [
  == Run-level clusters
  No refined MAG passed galah's quality filter (`--bin_min_completeness`, `--bin_max_contamination`), so the run has no clusters.
]
#datatable(s, "clusters")
#datatable(s, "cluster_abundance")
