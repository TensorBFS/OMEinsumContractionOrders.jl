#import "tensors.typ": *

#set page(width: auto, height: auto, margin: 5pt)

#figure(canvas(length: 0.6cm, {
  import draw: *
  let scircle(loc, radius, name) = {
    circle((loc.at(0)-0.1, loc.at(1)-0.1), radius: radius, fill:black, stroke: 1pt)
    circle(loc, radius: radius, stroke: (paint: black, thickness: 1pt), name: name, fill:white)
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
    circle((l.at(0), l.at(1)-s/2), radius:0.5, name: n, stroke: (paint: {black}, thickness: 1pt))
    content((l.at(0), l.at(1)-s/2), n)
  }
  for (a, b) in (("A", "B"), ("B", "C"), ("A", "C"), ("A", "D"), ("B", "E"), ("C", "F")){
    line(a, b, stroke: 1pt)
  }
  for (i, j, indice) in ((la, lb, "i"), (la, lc, "j"), (lb, lc, "k"), (la, ld, "l"), (lb, le, "m"), (lc, lf, "n")){
    // square(, width: 10pt, height: 10pt, paint: black, thickness: 1pt)
    circle(((i.at(0) + j.at(0))/2, (i.at(1) + j.at(1))/2 -s/2 ), radius:0.42, stroke: none, fill: white, thickness: 1pt)
    content(((i.at(0) + j.at(0))/2, (i.at(1) + j.at(1))/2 -s/2 ), indice)
  }
  set-origin((5, 0))
  content((0, 0), $arrow.r.long$)
  set-origin((5, 0))
  let scale = 2.5
  let li = (-0.5 * scale, 0.8 * scale - 1)
  let lj = (0.5 * scale, 0.8 * scale -1)
  let lk = (0, 0 -1)
  let lm = (-1 * scale, 0 - 1)
  let ln = (1 * scale, 0 - 1)
  let ll = (0, 1.6 * scale - 1)
  for (l, n) in ((li, "i"), (lj, "j"), (lk, "k"), (ll, "l"), (lm, "m"), (ln, "n")){
    circle((l.at(0), l.at(1)-s/2), radius:0.5, name: n, stroke:none, paint: {black}, thickness: 1pt)
    // ((l.at(0), l.at(1)-s/2), side: 0.5, paint: black, thickness: 1pt)
    content((l.at(0), l.at(1)-s/2), n)
  }
  for (a, b) in (("i", "j"), ("i", "k"), ("j", "k"), ("i", "m"), ("j", "n"), ("i", "l"), ("k", "m"), ("k", "n"), ("j", "l")){
    line(a, b, stroke: 1pt)
  }
  set-origin((5, 0))
  content((0, 0), $arrow.r.long$)
  set-origin((5, 0))
  let scale = 4.0
  let lijk = (0, 0)
  let lijl = (0,  scale * 0.7)
  let limk = (-0.5 * scale, -0.5 * scale)
  let ljkn = (0.5 * scale, -0.5 * scale)
  for (l, n) in ((lijk, "ijk"), (lijl, "ijl"), (limk, "imk"), (ljkn, "jkn")){
    circle((l.at(0), l.at(1)-s/2), width: 0.5, height: 0.5, name: n, stroke: (paint: {black}, thickness: 1pt))
    content((l.at(0), l.at(1)-s/2), n)
  }
  for (a, b) in (("ijk", "ijl"), ("ijk", "imk"), ("ijk", "jkn")){
    line(a, b, stroke: 1pt)
  }
}))