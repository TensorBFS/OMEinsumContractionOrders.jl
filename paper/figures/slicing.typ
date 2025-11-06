#import "@preview/cetz:0.4.0": canvas, draw, tree, coordinate
#set page(width: auto, height: auto, margin: 5pt)

#let labelnode(loc, label, name: none) = {
  import draw: *
  content(loc, text(black, label), align: center, fill:silver, frame:"rect", padding:0.07, stroke: none, name: name)
}
#canvas({
  import draw: *
  let d = 1.5
  let s(it) = text(11pt, it)
  let locs_labels = ((0.5 * d, -d), (0, -2 * d), (d, - 1.5 * d), (1.5 * d, -d), (2 * d, -2 * d), (d, -d))
  for (loc, t, name) in (((0.5 * d, -0.5 * d), s[$A$], "A"), ((1.5 * d, -0.5 * d), s[$B$], "B"), ((1.5 * d, -1.5 * d), s[$C$], "C"), ((0.5 * d, -1.5 * d), s[$D$], "D")) {
    circle(loc, radius: 0.3, name: name)
    content(loc, s[#t])
  }
  for ((loc, t), name) in locs_labels.zip((s[$i$], s[$j$], s[$k$], s[$l$], s[$m$], s[$n$])).zip(("i", "j", "k", "l", "m", "n")) {
    labelnode(loc, t, name: name)
  }
  for (src, dst) in (("A", "i"), ("D", "i"), ("D", "j"), ("D", "k"), ("C", "k"), ("B", "n"), ("B", "l"), ("C", "n"), ("C", "m"), ("C", "k"), ("C", "l"), ("D", "n")) {
    line(src, dst)
  }

  content((2.5 * d, -1.1 * d), text(16pt)[=])

  set-origin((3.5 * d, 0.0))

  content((-0.2 * d, -1.1 * d), text(12pt)[$sum_n$])

  for (loc, t, name) in (((0.5 * d, -0.5 * d), s[$A$], "A"), ((1.5 * d, -0.5 * d), s[$B_n$], "B"), ((1.5 * d, -1.5 * d), s[$C_n$], "C"), ((0.5 * d, -1.5 * d), s[$D_n$], "D")) {
    circle(loc, radius: 0.3, name: name)
    content(loc, s[#t])
  }
  for ((loc, t), name) in locs_labels.zip((s[$i$], s[$j$], s[$k$], s[$l$], s[$m$])).zip(("i", "j", "k", "l", "m")) {
    labelnode(loc, t, name: name)
  }
  for (src, dst) in (("A", "i"), ("D", "i"), ("D", "j"), ("D", "k"), ("C", "k"), ("B", "l"), ("C", "m"), ("C", "k"), ("C", "l")) {
    line(src, dst)
  }
})
