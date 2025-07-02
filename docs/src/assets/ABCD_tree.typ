#import "tensors.typ": *

#set page(width: auto, height: auto, margin: 5pt)

#align(center, canvas(length:1.0cm, {
  import draw: *
  set-origin((4, 0.35))
  let DY = 1.2
  let DX1 = 1.5
  let DX2 = 0.9
  let root = (0, DY)
  let left = (-DX1, 0)
  let right = (DX1, 0)
  let left_left = (-DX1 - DX2, -DY)
  let left_right = (-DX1 + DX2, -DY)
  let right_left = (DX1 - DX2, -DY)
  let right_right = (DX1 + DX2, -DY)

  for (l, t, lb) in ((root, [$s$], "C"), (left, [$A B$], "A"), (right, [$C D$], "B"), (left_left, [$A$], "T_1"), (left_right, [$B$], "T_2"), (right_left, [$C$], "T_3"), (right_right, [$D$], "T_4")){
    tensor(l, lb, text(11pt, t))
  }

  for (a, b) in (("C", "A"), ("C", "B"), ("A", "T_1"), ("A", "T_2"), ("B", "T_3"), ("B", "T_4")){
    line(a, b)
  }


}))