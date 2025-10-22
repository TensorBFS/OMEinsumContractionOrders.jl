#import "@preview/cetz:0.4.0": canvas, draw, tree, coordinate
#set page(width: auto, height: auto, margin: 5pt)

#canvas({
  import draw: *
  let dy = 1.5
  let boxed(it, width: 80pt, fill: none) = box(width: width, stroke: black, fill:fill, inset: 5pt, align(left)[#it])
  content((-5, 0), boxed(
[`EinCode`
- `ixs`
- `iy`]), name: "EinCode")
  content((0, 0), boxed([
    `NestedEinsum`
- `args`
- `tensorindex`
- `eins`
]), name: "NestedEinsum")
  content((5, 0), boxed([
    `SlicedEinsum`
- `slicing`
- `eins`
]), name: "SlicedEinsum")
  line("EinCode", "NestedEinsum", mark: (end: "straight"), name: "optimize")
  line("NestedEinsum", "SlicedEinsum", mark: (end: "straight"), name: "slice")
  content((rel: (-0.1, 0.2), to: "slice.mid"), text(8pt)[`slice_code`])
  content((rel: (-0.1, 0.2), to: "optimize.mid"), text(8pt)[`optimize_code`])
})
