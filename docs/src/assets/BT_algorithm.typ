#import "tensors.typ": *

#set page(width: auto, height: auto, margin: 5pt)

#figure(canvas({
  import draw: *
  let locs = ((-1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1), (1, 1), (1, 0), (1, -1))
  let labels_colors = (("A", black), ("B", red), ("C", red), ("D", black), ("E", red), ("F", black), ("G", black), ("H", black))
  for (loc, lc) in locs.zip(labels_colors) {
    let (label, color) = lc
    tensor(loc, label, label, color: color)
  }

  for (src, dst) in (("A", "B"), ("A", "C"), ("B", "F"), ("C", "D"), ("B", "G"), ("B", "C"), ("G", "F"), ("G", "H"), ("E", "H"), ("D", "E"), ("G", "E"), ("C", "E")) {
    line(src, dst)
  }

  content((2, 0), $arrow$)

  let sx = 4
  let sy = 1
  let locs_1 = ((-1 + sx, 1 + sy), (0 + sx, 1 + sy), (-1 + sx, 0 + sy))
  let labels_colors_1 = (("A", black), ("B", red), ("C", red))
  for (loc, lc) in locs_1.zip(labels_colors_1) {
    let (label, color) = lc
    tensor(loc, label, label, color: color)
  }

  for (src, dst) in (("A", "B"), ("A", "C"), ("B", "C")) {
    line(src, dst)
  }

  // cde
  let sx = 4
  let sy = -1
  let locs_2 = ((-1 + sx, 0 + sy), (-1 + sx, -1 + sy), (0 + sx, -1 + sy))
  let labels_colors_2 = (("C", red), ("D", black), ("E", red))
  for (loc, lc) in locs_2.zip(labels_colors_2) {
    let (label, color) = lc
    tensor(loc, label, label, color: color)
  }

  for (src, dst) in (("C", "D"), ("C", "E"), ("D", "E")) {
    line(src, dst)
  }

  // bce
  let sx = 5
  let sy = -0
  let locs_3 = ((sx, 1 + sy), (sx - 1, sy), (sx, sy - 1))
  let labels_colors_3 = (("B", red), ("C", red), ("E", red))
  for (loc, lc) in locs_3.zip(labels_colors_3) {
    let (label, color) = lc
    tensor(loc, label, label, color: color)
  }

  for (src, dst) in (("B", "C"), ("B", "E"), ("C", "E")) {
    line(src, dst)
  }

  // befgh
  let sx = 6
  let sy = 0
  let locs_4 = ((0 + sx, 1 + sy), (0 + sx, -1 + sy), (1 + sx, 1 + sy), (1 + sx, 0 + sy), (1 + sx, -1 + sy))
  let labels_colors_4 = (("B", red), ("E", red), ("F", black), ("G", black), ("H", black))
  for (loc, lc) in locs_4.zip(labels_colors_4) {
    let (label, color) = lc
    tensor(loc, label, label, color: color)
  }

  for (src, dst) in (("B", "F"), ("B", "F"), ("F", "G"), ("G", "H"), ("E", "H"), ("B", "G"), ("G", "E")) {
    line(src, dst)
  }
}))