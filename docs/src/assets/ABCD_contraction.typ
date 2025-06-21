#import "tensors.typ": *

#set page(width: auto, height: auto, margin: 5pt)

#figure(canvas({
  import draw: *
  tensor((1, 1), "A", "A")
  tensor((3, 1), "B", "B")
  tensor((3, 3), "C", "C")
  tensor((1, 3), "D", "D")
  labeledge("A", "B", "i")
  labeledge("B", "C", "j")
  labeledge("C", "D", "k")
  labeledge("D", "A", "l")

  content((4, 2), $arrow$)

  tensor((5, 1), "AB", "AB")
  tensor((5, 3), "CD", "CD")
  labeledge("AB", "CD", "jl")

  content((6, 2), $arrow$)

  tensor((7, 2), "s", "s")
}))