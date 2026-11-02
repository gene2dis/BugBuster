#import "lib.typ": *
#let s = sec("taxonomy")

= Taxonomy
Community composition from #s.facts.db#if s.facts.n_species != none [; #s.facts.n_species species were detected in at least one sample]. The heatmaps show the union of each sample's most abundant taxa on a log scale (#figref(s, "taxonomy_species_heatmap"), #figref(s, "taxonomy_genus_heatmap")); the bars show proportions (#figref(s, "taxonomy_genus_bars")).

#fig(s, "taxonomy_species_heatmap")
#fig(s, "taxonomy_genus_heatmap")
#fig(s, "taxonomy_genus_bars")
#if has(s, "figures", "watchlist") [
  == Watchlist
  Priority pathogens and other watchlist taxa found in the samples (#figref(s, "watchlist"), #tabref(s, "watchlist")). Detection is not culture.
]
#fig(s, "watchlist")
#datatable(s, "watchlist")
#datatable(s, "prevalence")
