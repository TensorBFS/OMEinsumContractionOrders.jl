#import "tensors.typ": *

#set page(width: auto, height: auto, margin: 5pt)

#figure(canvas({
  import draw: *
  let locs = ((-1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1), (1, 1), (1, 0), (1, -1))
  let labels = ("A", "B", "C", "D", "E", "F", "G", "H")
  for (loc, label) in locs.zip(labels) {
    tensor(loc, label, label)
  }
  for (src, dst) in (("B", "C"), ("G", "F"), ("G", "H"), ("E", "H"), ("D", "E")) {
    line(src, dst)
  }
  for (src, dst) in (("B", "C"), ("G", "F"), ("G", "H"), ("E", "H"), ("D", "E"), ("G", "E")) {
    line(src, dst)
  }
}))