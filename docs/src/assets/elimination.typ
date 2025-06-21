#import "tensors.typ": *

#set page(width: auto, height: auto, margin: 5pt)

#figure(canvas({
  import draw: *
  
  let locs = ((0, 1), (-1, -0.5), (1, -0.5), (0, 2.4), (-2.1, -1.5), (2.1, -1.5))
  for (loc, name) in locs.zip(("A", "B", "C", "D", "E", "F")) {
    tensor(loc, name, name)
  }

  for (src, dst, label) in (("A", "B", "i"), ("B", "C", "k"), ("A", "C", "j"), ("D", "A", "l"), ("B", "E", "m"), ("C", "F", "n")) {
    labeledge(src, dst, label)
  }

  content((3, 0.5), $arrow$)
  content((3, 1), $n$)

  set-origin((6, 0))
  let locs_2 = ((0, 1), (-1, -0.5), (1, -0.5), (0, 2.4), (-2.1, -1.5))
  for (loc, name) in locs_2.zip(("A", "B", "C", "D", "E")) {
    tensor(loc, name, name)
  }

  for (src, dst, label) in (("A", "B", "i"), ("B", "C", "k"), ("A", "C", "j"), ("D", "A", "l"), ("B", "E", "m")) {
    labeledge(src, dst, label)
  }

  content((2, 0.5), $arrow$)
  content((2, 1), $m$)

  set-origin((4, 0))

  let locs_3 = ((0, 1), (-1, -0.5), (1, -0.5), (0, 2.4))
  for (loc, name) in locs_3.zip(("A", "B", "C", "D")) {
    tensor(loc, name, name)
  }

  for (src, dst, label) in (("A", "B", "i"), ("B", "C", "k"), ("A", "C", "j"), ("D", "A", "l")) {
    labeledge(src, dst, label)
  }

  content((2, 0.5), $arrow$)
  content((2, 1), $k$)

  set-origin((3.5, 0))
  let locs_4 = ((0, 0.5), (0, 2), (0, -1))
  for (loc, name) in locs_4.zip(("A", "D", "B")) {
    tensor(loc, name, name)
  }
  for (src, dst, label) in (("A", "D", "l"), ("A", "B", "ij")) {
    labeledge(src, dst, label)
  }

  content((1, 0.5), $arrow$)
  content((1, 1), $i j l$)

  content((2, 0.5), $s$)
}))