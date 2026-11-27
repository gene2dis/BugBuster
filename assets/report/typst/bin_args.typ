#import "lib.typ": *
#let s = sec("bin_args")

= ARGs of MAGs
The #s.facts.n_genes ARG genes of the MAGs were clustered at 90-100 % identity; #s.facts.n_shared clusters at 90 % are found in more than one MAG#refs(s, tabs: ("sharing",)). A cluster in several MAGs is the same sequence in different genomes: consistent with transfer, not evidence of it alone.

#if "card" in s.facts [
  argNorm maps the MAGs' DeepARG hits to CARD: of the #s.facts.n_card_binned hits on a MAG, #s.facts.n_card_loose rest on a Loose mapping (often a homolog) and get no tier#if s.facts.tier1_genes != "" [; tier 1: #s.facts.tier1_genes]. Each hit's CARD gene, tier and drug classes are in arg_deeparg_card.tsv (the rows with a bin).
]

#datatable(s, "sharing")
#datatable(s, "long")
