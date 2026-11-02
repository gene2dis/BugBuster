#import "lib.typ": *
#show: setup

#block(below: 4pt, text(19pt, weight: "bold", d.run.title))
#text(9.5pt, fill: muted)[Generated #d.run.generated#if d.run.pipeline_version != "" [ by BugBuster #d.run.pipeline_version] from #d.run.n_samples samples (#d.run.batches.join(", ")).]
#v(6pt)
#outline(depth: 2, indent: 1em)

// one file per section present, in the order Python gives; the data dictionary is in the HTML report only
#for s in d.sections {
  if s.id != "dictionary" { include (s.id + ".typ") }
}
