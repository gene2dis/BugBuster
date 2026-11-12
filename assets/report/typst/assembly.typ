#import "lib.typ": *
#let s = sec("assembly")

= Assembly
#s.facts.n_assemblies MEGAHIT assembl#if s.facts.n_assemblies == 1 [y] else [ies] (#tabref(s, "stats")). Only contigs kept by the length filter go on to binning, contig taxonomy, ARGs and functions, and reads that are not on contigs are outside every contig-based table (#tabref(s, "read_use")).

#if "n_groups" in s.facts [
  The co-assemblies are one per `--coassembly_group_by` group (#s.facts.n_groups group#if s.facts.n_groups != 1 [s]): each pools its group's samples#if has(s, "tables", "read_use") [, and #tabref(s, "read_use") gives each sample's group].
]

#if "bbnorm" in s.facts [
  BBNorm normalised the reads given to MEGAHIT (`--assembly_read_reduction normalize`): only the assembler input was reduced, and every alignment, depth and count uses the full clean reads.
]

#datatable(s, "stats")
#datatable(s, "read_use")
