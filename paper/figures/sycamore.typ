#import "@preview/cetz:0.4.0": canvas, draw, tree
#import "@preview/cetz-plot:0.1.2": plot
#import "report.typ": plot-compare, grouped_data
#set page(width: auto, height: auto, margin: 5pt)

//#let combined_keys = grouped_data.keys().sorted()
#figure({
  let dataset = grouped_data.at("quantumcircuit/sycamore_53_20_0")
  let a = 1
  let b = 0
  let c = 0
  canvas(length: 1cm, {
    plot-compare(dataset, a: a, b: b, c: c)
  })},
)