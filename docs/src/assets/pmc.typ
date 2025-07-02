#import "tensors.typ": *

#set page(width: auto, height: auto, margin: 5pt)

#figure(canvas({
  import draw: *
  let locs = ((-1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1), (1, 1), (1, 0), (1, -1))
  let labels_colors = (("A", black), ("B", red), ("C", red), ("D", black), ("E", red), ("F", black), ("G", black), ("H", black))
  for (loc, lc) in locs.zip(labels_colors) {
    let (label, color) = lc
    tensor(loc, label, label, color: black)
  }

  for (src, dst) in (("A", "B"), ("A", "C"), ("B", "F"), ("C", "D"), ("B", "G"), ("B", "C"), ("G", "F"), ("G", "H"), ("E", "H"), ("D", "E"), ("G", "E"), ("C", "E")) {
    line(src, dst)
  }

  content((2.5, 1.2), $arrow$)
  content((2.5, - 1.2), $arrow$)

  set-origin((5, 2))
  let labels_colors = (("A", black), ("B", red), ("C", red), ("D", black), ("E", red), ("F", black), ("G", black), ("H", black))
  for (loc, lc) in locs.zip(labels_colors) {
    let (label, color) = lc
    tensor(loc, label, label, color: black)
  }

  for (src, dst) in (("A", "B"), ("A", "C"), ("B", "F"), ("C", "D"), ("B", "G"), ("B", "C"), ("G", "F"), ("G", "H"), ("E", "H"), ("D", "E"), ("G", "E"), ("C", "E")) {
    line(src, dst)
  }

  line("C", "G")

  set-origin((0, - 4))
  let labels_colors = (("A", black), ("B", red), ("C", red), ("D", black), ("E", red), ("F", black), ("G", black), ("H", black))
  for (loc, lc) in locs.zip(labels_colors) {
    let (label, color) = lc
    tensor(loc, label, label, color: black)
  }

  for (src, dst) in (("A", "B"), ("A", "C"), ("B", "F"), ("C", "D"), ("B", "G"), ("B", "C"), ("G", "F"), ("G", "H"), ("E", "H"), ("D", "E"), ("G", "E"), ("C", "E")) {
    line(src, dst)
  }
  line("B", "E")
}))