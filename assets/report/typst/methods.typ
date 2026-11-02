#import "lib.typ": *
#let s = sec("methods")

= Methods
== Reading this report
#bullets(s.lists.units)
#if s.lists.arg_methods.len() > 0 [
  *ARG methods.* The ARG figures and matrices use #s.facts.primary.
  #list(..s.lists.arg_methods.map(m => [*#m.term* #m.text]))
]
#datatable(s, "tiers", widths: (8%, 14%, 30%, 1fr))

== Caveats
#bullets(s.lists.caveats)

== Inputs
#for x in s.lists.inputs [#raw(x) \ ]

== Software and databases
#if has(s, "tables", "tools") {
  let t = s.tables.tools
  let p = t.pdf
  let cites(i) = t.keys.at(i).map(k => cite(label(k))).join()
  let tool-col = p.header.position(h => h == "Tool")
  let skip = p.header.position(h => h == "Ref.")
  block(breakable: false, below: 4pt, text(8pt)[*Table #t.no.* #t.caption])
  set text(6.5pt)
  let cols = range(p.header.len()).filter(i => i != skip)
  table(columns: cols.len(), stroke: none, inset: (x: 3.5pt, y: 2.2pt),
    fill: (_, y) => if y > 0 and calc.even(y) { luma(245) },
    table.header(..cols.map(i => text(weight: "bold", p.header.at(i)))),
    table.hline(stroke: 0.5pt),
    ..p.rows.enumerate().map(((r, row)) => cols.map(i =>
      if i == tool-col [#row.at(i)#if t.keys.at(r).len() > 0 [ #cites(r)]] else [#row.at(i).replace(",", ",\u{200b}")])).flatten())
}
#datatable(s, "watchlist")
#if s.lists.references.len() > 0 {
  bibliography("references.bib", title: "References", style: "nature")
}
