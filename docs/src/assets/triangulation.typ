#import "tensors.typ": *

#set page(width: auto, height: auto, margin: 5pt)

#canvas({
  import draw: *
  let s(it) = text(11pt, it)
  let locs_labels = ((0, -0.2), (0, -2), (1.5, -1), (2, 0), (2, -2.2), (3, -1))

  for ((loc, t), name) in locs_labels.zip((s[$a$], s[$b$], s[$c$], s[$d$], s[$e$], s[$f$])).zip(("a", "b", "c", "d", "e", "f")) {
    circle(loc, radius: 0.4, name: name)
    content(loc, s[#t])
  }

  for (src, dst) in (("a", "d"), ("b", "e"), ("c", "d"), ("c", "e"), ("d", "f"), ("e", "f"), ("b", "a")) {
    line(src, dst)
  }

  set-origin((7, 0))
  for ((loc, t), name) in locs_labels.zip((s[$a$], s[$b$], s[$c$], s[$d$], s[$e$], s[$f$])).zip(("a2", "b2", "c2", "d2", "e2", "f2")) {
    circle(loc.map(v => v), radius: 0.4, name: name)
    content(loc.map(v => v ), s[#t])
  }

  for (src, dst) in (("a2", "d2"), ("b2", "e2"), ("c2", "d2"), ("c2", "e2"), ("d2", "f2"), ("e2", "f2"), ("b2", "a2"), ("e2", "a2"), ("d2", "e2")) {
    line(src, dst)
  }
})