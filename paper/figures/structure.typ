#import "@preview/cetz:0.4.0": canvas, draw, tree, coordinate
#set page(width: auto, height: auto, margin: 5pt)

#canvas({
  import draw: *
  let dy = 1.5
  let boxed(it, width: 100pt) = box(width: width, radius: 3pt, stroke: black, inset: 5pt, align(center)[#it])
  content((-2, -0.5 * dy), boxed([CliqueTrees]))
  content((2, -0.5 * dy), boxed([TreeWidthSolver]))
  content((0, 0), boxed([OMEinsumContractionOrders], width: 214pt), name: "OMEinsumContractionOrders")
  content((0, dy), boxed([OMEinsum]), name: "OMEinsum")
  content((4, dy), boxed([NCTSSoS]), name: "NCTSSoS")
  content((0, 3*dy), boxed([GenericTensorNetworks], width: 150pt))
  content((2, 2*dy), boxed([TensorQEC]))
  content((-2, 2*dy), boxed([TensorInference]))
  content((-2, 2.5 *dy), boxed([Yao]))
  content((2, 2.5*dy), boxed([TensorBranching]))
  rect((-4, 1.7 * dy), (4, 3.3*dy), stroke: (paint: black, dash: "dashed"), name: "Applications")
  rect((-4, -0.8 * dy), (4, 0.3*dy), stroke: (paint: black, dash: "dashed"), name: "Algorithms")
  content((-6, 2.5 * dy), [Applications])
  content((-6, -0.25 * dy), [Algorithms])
  content((-6, 1 * dy), [Modeling])
  line("Algorithms", "OMEinsum", mark: (end: "straight"))
  line("Algorithms", "NCTSSoS", mark: (end: "straight"))
  line("OMEinsum", "Applications", mark: (end: "straight"))
})
