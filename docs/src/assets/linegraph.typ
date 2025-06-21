#import "tensors.typ": *

#set page(width: auto, height: auto, margin: 5pt)

#grid(columns: 5, rows: 1, gutter: 1.5cm, 
align( center + horizon,
canvas(length: 1cm, {
  import draw: *
  let scircle(loc, radius, name) = {
    circle((loc.at(0)-0.1, loc.at(1)-0.1), radius: radius, fill:black, stroke: 3pt)
    circle(loc, radius: radius, stroke: (paint: black, thickness: 3pt), name: name, fill:white)
  }
  let s = 1.0
  let dy = 3.0
  let la = (0, 2*s)
  let lb = (-1.5 * s, -0.6 * s)
  let lc = (1.5 * s, -0.6 * s)
  let ld = (0, 4.3 * s)
  let le = (-3 * s, - 2.5 * s)
  let lf = (3 * s, - 2.5 * s)
  for (l, n) in ((la, "A"), (lb, "B"), (lc, "C"), (ld, "D"), (le, "E"), (lf, "F")){
    circle((l.at(0), l.at(1)-s/2), radius:0.5, name: n, stroke: (paint: {black}, thickness: 3pt))
    content((l.at(0), l.at(1)-s/2), n)
  }
  for (a, b) in (("A", "B"), ("B", "C"), ("A", "C"), ("A", "D"), ("B", "E"), ("C", "F")){
    line(a, b, stroke: 3pt)
  }
  for (i, j, indice) in ((la, lb, "i"), (la, lc, "j"), (lb, lc, "k"), (la, ld, "l"), (lb, le, "m"), (lc, lf, "n")){
    // square(, width: 10pt, height: 10pt, paint: black, thickness: 3pt)
    circle(((i.at(0) + j.at(0))/2, (i.at(1) + j.at(1))/2 -s/2 ), radius:0.42, stroke: none, fill: white, thickness: 3pt)
    content(((i.at(0) + j.at(0))/2, (i.at(1) + j.at(1))/2 -s/2 ), indice)
  }
})),
align( center + horizon,
$arrow.r.long$
),
align( center + horizon,
canvas(length: 1cm, {
  import draw: *
  let s = 1.0
  let scale = 2.5
  let li = (-0.5 * scale, 0.8 * scale)
  let lj = (0.5 * scale, 0.8 * scale)
  let lk = (0, 0)
  let lm = (-1 * scale, 0)
  let ln = (1 * scale, 0)
  let ll = (0, 1.6 * scale)
  for (l, n) in ((li, "i"), (lj, "j"), (lk, "k"), (ll, "l"), (lm, "m"), (ln, "n")){
    circle((l.at(0), l.at(1)-s/2), radius:0.5, name: n, stroke:none, paint: {black}, thickness: 3pt)
    // ((l.at(0), l.at(1)-s/2), side: 0.5, paint: black, thickness: 3pt)
    content((l.at(0), l.at(1)-s/2), n)
  }
  for (a, b) in (("i", "j"), ("i", "k"), ("j", "k"), ("i", "m"), ("j", "n"), ("i", "l"), ("k", "m"), ("k", "n"), ("j", "l")){
    line(a, b, stroke: 3pt)
  }
})),
align( center + horizon,
$arrow.r.long$
),
align( center + horizon,
canvas(length: 0.8cm, {
  import draw: *
  let scircle(loc, radius, name) = {
    circle((loc.at(0)-0.1, loc.at(1)-0.1), radius: radius, fill:black, stroke: 3pt)
    circle(loc, radius: radius, stroke: (paint: black, thickness: 3pt), name: name, fill:white)
  }
  let s = 1.0
  let scale = 4.0
  let lijk = (0, 0)
  let lijl = (0,  scale * 0.7)
  let limk = (-0.5 * scale, -0.5 * scale)
  let ljkn = (0.5 * scale, -0.5 * scale)
  for (l, n) in ((lijk, "ijk"), (lijl, "ijl"), (limk, "imk"), (ljkn, "jkn")){
    circle((l.at(0), l.at(1)-s/2), width: 0.5, height: 0.5, name: n, stroke: (paint: {black}, thickness: 3pt))
    content((l.at(0), l.at(1)-s/2), n)
  }
  for (a, b) in (("ijk", "ijl"), ("ijk", "imk"), ("ijk", "jkn")){
    line(a, b, stroke: 3pt)
  }
})),
)