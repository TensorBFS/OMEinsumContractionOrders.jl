#import "@preview/cetz:0.4.0": canvas, draw, tree
#import "@preview/cetz-plot:0.1.2": plot
#import "report.typ": plot-compare, grouped_data
#set page(width: auto, height: auto, margin: 5pt)

// the cost function is a * sc + b * tc + c * rwc
#let plot-slicing(points, a: 1, b: 0, c: 0, x-label: "Log2 Space Complexity", y-label: "Log10 Time Complexity") = {
 plot.plot(
  size: (12, 8),
  x-label: x-label,
  y-label: y-label,
  x-min: 30.0,
  x-max: 52.0,
  y-min: 18.0,
  y-max: 24.0,
  legend: "inner-north-east",
    { 
      plot.add(points, style: (stroke: none, fill: red), mark: "o", mark-size: 0.15, mark-style: (fill: red, stroke: red), label: "TreeSA Slicer")
    }
  )
}

#let tcs = [23.39903483741713, 22.850373303040385, 21.286770346844808, 20.511290170386072, 19.839694209644776, 18.985439210394247, 18.6068528031432, 18.425122025308436, 18.26920496351226, 18.2026439528554, 18.13556549285669]
#let scs = [31.0, 33.0, 35.0, 37.0, 39.0, 41.0, 43.0, 45.0, 47.0, 49.0, 51.0]

#let points = ((31.0, 23.39903483741713), (33.0, 22.850373303040385), (35.0, 21.286770346844808), (37.0, 20.511290170386072), (39.0, 19.839694209644776), (41.0, 18.985439210394247), (43.0, 18.6068528031432), (45.0, 18.425122025308436), (47.0, 18.26920496351226), (49.0, 18.2026439528554), (51.0, 18.13556549285669))

// #let points = ((31.0, 33.0, 35.0, 37.0, 39.0, 41.0, 43.0, 45.0, 47.0, 49.0, 51.0), (23.39903483741713, 22.850373303040385, 21.286770346844808, 20.511290170386072, 19.839694209644776, 18.985439210394247, 18.6068528031432, 18.425122025308436, 18.26920496351226, 18.2026439528554, 18.13556549285669))


#canvas(length: 1cm, {
  import draw: bezier, content
  plot-slicing(points, a: 0, b: 1, c: 0)
})