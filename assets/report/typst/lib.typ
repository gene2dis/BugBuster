// Shared layout of the BugBuster summary PDF. Everything shown comes from data.json, written by
// bin/summary_report.py (numbers, captions, table rows, figure and table numbers); the section files hold the
// prose only.
#let d = json("data.json")
#let accent = rgb("#0072B2")
#let muted = luma(95)
#let sec(id) = d.sections.find(x => x.id == id)
#let has(s, kind, key) = s.at(kind).at(key, default: none) != none

#let setup(body) = {
  set document(title: d.run.title, date: none)
  set page(paper: "a4", margin: (x: 13.5mm, top: 19mm, bottom: 17mm),
    header: context { if counter(page).get().first() > 1 [
      #set text(7pt, fill: muted); #d.run.title #h(1fr) #d.run.generated ] },
    footer: context [#set text(7pt, fill: muted); #h(1fr) #counter(page).display("1 / 1", both: true)])
  set text(font: "Inter", size: 9pt, lang: "en")
  set par(justify: true, leading: 0.6em, spacing: 0.9em)
  set heading(numbering: "1.1")
  show heading.where(level: 1): it => { pagebreak(weak: true); set text(15pt); block(below: 9pt, it) }
  show heading.where(level: 2): set text(11pt)
  show heading: set block(above: 1.2em, below: 0.6em)
  show link: set text(fill: accent)
  set list(indent: 0.4em, spacing: 0.55em)
  body
}

// "Figure 3" / "Table 5" for the prose (numbers from Python, as in the HTML)
#let figref(s, key) = if has(s, "figures", key) [Figure #s.figures.at(key).no] else []
#let tabref(s, slot) = if has(s, "tables", slot) [Table #s.tables.at(slot).no] else []

#let isnum(x) = {
  let t = x.replace(",", "").replace("%", "").trim()
  t.len() > 0 and t.match(regex("^[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?$")) != none
}

// a numbered table: caption above, header repeated on every page, zebra rows, numbers right-aligned, the
// row/column cap note below; landscape page when it is wide
#let datatable(s, slot, size: 6.5pt, widths: auto) = {
  if not has(s, "tables", slot) { return }
  let t = s.tables.at(slot)
  if t.pdf == none { return }
  let p = t.pdf
  if p.rows.len() == 0 {
    // an empty table keeps its number, as in the HTML
    block(above: 8pt, below: 10pt, text(8pt)[*Table #t.no.* #t.caption #text(fill: muted)[No rows.]])
    return
  }
  let n = p.header.len()
  let body = {
    block(breakable: false, below: 4pt, text(8pt)[*Table #t.no.* #t.caption])
    set text(size)
    table(
      columns: if widths == auto { n } else { widths },
      stroke: none, inset: (x: 3.5pt, y: 2.2pt),
      fill: (_, y) => if y > 0 and calc.even(y) { luma(245) },
      align: (x, y) => if y > 0 and isnum(p.rows.at(y - 1).at(x)) { right } else { left },
      table.header(..p.header.map(h => text(weight: "bold", h.replace("_", "_\u{200b}")))),
      table.hline(stroke: 0.5pt),
      ..p.rows.flatten().map(c => c.replace("_", "_\u{200b}"))
    )
    if p.note != "" { block(above: 3pt, text(7pt, fill: muted, p.note)) }
  }
  if p.landscape { page(flipped: true, body) } else { block(above: 8pt, below: 10pt, body) }
}

// a figure from its SVG at the size it was drawn (89 or 183 mm), caption below, never split from it
#let fig(s, key) = {
  if not has(s, "figures", key) { return }
  let f = s.figures.at(key)
  let img = if f.height_mm > 222 { image("figures/" + key + ".svg", height: 222mm) } else {
    image("figures/" + key + ".svg", width: calc.min(f.width_mm, 183) * 1mm) }
  block(breakable: false, above: 10pt, below: 12pt, {
    align(center, img)
    v(3pt)
    set text(7.5pt)
    set par(justify: false)
    [*Figure #f.no. #f.title.* #f.caption]
  })
}

#let bullets(items) = list(..items.map(x => [#x]))
