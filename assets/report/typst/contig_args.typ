#import "lib.typ": *
#let s = sec("contig_args")

= ARGs on contigs
DeepARG found #s.facts.n_genes ARG genes on the contigs of #s.facts.n_assemblies assembl#if s.facts.n_assemblies == 1 [y] else [ies]. Each is shown with its contig's taxon (Metabuli), bin and geNomad class: physical linkage, unlike the read-level candidate hosts. These are gene counts on contigs, not read abundances#refs(s, figs: ("arg_contig_hosts",), tabs: ("hosts",)).

#let f = s.facts
#if "mobile_genes" in f [
  Of the #f.mobile_genes ARG genes#if f.mobile_partial [ in the results directories where geNomad ran], geNomad places #f.on_plasmid on plasmid contigs, #f.on_virus on virus contigs and #f.in_provirus in proviruses#refs(s, tabs: ("mobile",)). Its plasmid class includes integrative conjugative elements and transposon cargo: an ARG "on a plasmid" sits on a mobile-element contig, not necessarily a plasmid.#if "checkv_viral" in f [ Of the ARG genes on virus contigs or in proviruses (#f.checkv_viral#if f.checkv_partial [, in the results directories where CheckV ran]), CheckV rates the virus complete or high quality for #f.checkv_high, medium for #f.checkv_medium, low for #f.checkv_low and not determined for #f.checkv_not_determined (too short or too novel to judge).]
]

#if f.at("card", default: false) [
  argNorm maps each DeepARG hit's best database hit to a CARD ARO: #f.n_card_mapped hits by a Perfect, Strict or Manual mapping, which gives them the CARD gene and reporting tier of the read-level calls#if f.tier1_genes != "" [ (tier 1: #f.tier1_genes)], and #f.n_card_loose by a Loose mapping (often a homolog), which get no tier#refs(s, tabs: ("card",)).
]

#fig(s, "arg_contig_hosts")
#datatable(s, "hosts")
#datatable(s, "mobile")
#datatable(s, "card")
#datatable(s, "long")
