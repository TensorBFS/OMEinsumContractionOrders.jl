#import "tensors.typ": *

#set page(width: auto, height: auto, margin: 5pt)

#canvas({
  import draw: *

  let locs = ((0, 1), (-1, -0.5), (1, -0.5))
  for (loc, name) in locs.zip(("A", "B", "C")) {
    tensor(loc, name, name)
  }

  let locs_empty = ((0, 2.4), (-2.1, -1.5), (2.1, -1.5))
  for (loc, name) in locs_empty.zip(("D", "E", "F")) {
    circle(loc, radius: 0, name: name)
  }

  for (src, dst, label) in (("A", "B", "i"), ("B", "C", "k"), ("A", "C", "j"), ("D", "A", "l"), ("B", "E", "m"), ("C", "F", "n")) {
    labeledge(src, dst, label)
  }

  content((3, 0.5), $arrow$)

  set-origin((6, 0))
  let s = 1.5
  let locs = ((0, 1 * s), (-1 * s, -0.5 * s), (1 * s, -0.5 * s), (0, 0))
  for (loc, name) in locs.zip(("A", "B", "C", "D")) {
    tensor(loc, name, name)
  }
  for (src, dst, label) in (("A", "B", "i"), ("B", "C", "k"), ("A", "C", "j"), ("D", "A", "l"), ("B", "D", "m"), ("C", "D", "n")) {
    labeledge(src, dst, label)
  }
})