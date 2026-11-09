#import "lib.typ": *
#let s = sec("mags")
#let f = s.facts

= MAGs
#f.n_mags refined MAGs: #f.high high, #f.medium medium and #f.low low quality by MIMAG (#figref(s, "mag_quality"), #tabref(s, "mags")). The High / Mid / Low classes of the run's binning plots are another, older scheme. #figref(s, "mag_abundance_heatmap") shows the share of each sample's sequenced bases on each MAG.

#fig(s, "mag_quality")
#fig(s, "mag_abundance_heatmap")
#datatable(s, "mags")
#datatable(s, "abundance")
#if "n_clusters" in f [
  == Run-level clusters
  #f.n_clustered MAGs fall into #f.n_clusters clusters at #f.derep_ani % ANI (galah; the representative of each has the highest completeness − 5 × contamination). #tabref(s, "cluster_abundance") gives each representative's share of every sample's reads (CoverM). The clusters cover this run's samples only.
]
#datatable(s, "clusters")
#datatable(s, "cluster_abundance")
