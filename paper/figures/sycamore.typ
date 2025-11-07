#import "@preview/cetz:0.4.0": canvas, draw, tree
#import "@preview/cetz-plot:0.1.2": plot
#import "report.typ": plot-compare, grouped_data
#set page(width: auto, height: auto, margin: 5pt)

// Read the merged summary data
#let all_data = json("data/summary.json")

// Function to extract instance name from path
#let get_instance_name(instance_path) = {
  let parts = instance_path.split("/")
  let filename = parts.last()
  filename.replace(".json", "")
}

// Filter data for sycamore_53_20_0 only
#let filtered_data = all_data.filter(entry => {
  let instance_name = get_instance_name(entry.instance)
  instance_name == "sycamore_53_20_0"
})

// the cost function is a * sc + b * tc + c * rwc
#let plot-compare(dataset, x-max: auto, y-min: -3, y-max: 3.5, a: 1, b: 0, c: 0, x-label: "Contraction Cost", y-label: "Log10 Computing Time (seconds)") = {
 plot.plot(
  size: (12, 8),
  x-label: x-label,
  y-label: y-label,
  x-min: auto,
  y-min: y-min,
  y-max: y-max,
  x-max: x-max,
  legend: "inner-north-east",
    {
      // Automatically group points by optimizer
      let optimizer_groups = (:)
      
      // Collect all data points grouped by optimizer
      for entry in dataset {
        let sc = entry.contraction_complexity.sc
        let tc = entry.contraction_complexity.tc
        let rwc = entry.contraction_complexity.rwc
        let cost = a * sc + b * tc + c * rwc
        let time = calc.log(entry.time_elapsed, base: 10)
        let optimizer = entry.optimizer
        
        if optimizer not in optimizer_groups {
          optimizer_groups.insert(optimizer, ())
        }
        optimizer_groups.at(optimizer).push((cost, time))
      }
      
      // Define colors and markers for different optimizers
      let colors = (red, blue, green, purple, orange, black, aqua, gray, teal, maroon, navy, olive)
      let markers = ("o", "square", "triangle", "o", "square", "+", "x")
      
      // Plot each optimizer group (sorted by optimizer name)
      let sorted_optimizers = optimizer_groups.keys().sorted()
      let i = 0
      for optimizer in sorted_optimizers {
        let points = optimizer_groups.at(optimizer)
        let color = colors.at(calc.rem(i, colors.len()))
        let marker = markers.at(calc.rem(i, markers.len()))
        
        plot.add(
          points,
          style: (stroke: none, fill: color),
          mark: marker,
          mark-size: 0.15,
          mark-style: (fill: color, stroke: color),
          label: optimizer
        )
        i += 1
      }
    }
  )
}

#grid(columns: 2, gutter: 15pt,
  canvas(length: 1cm, {
    import draw: bezier, content
    plot-compare(filtered_data, a: 0, b: 1, c: 0, y-min: -2, y-max: 4, x-label: "Log10 FLOPs")
    bezier((0, 5), (5, 0.7), (2, 3.2), stroke: (dash: "dashed"))
    content((1.8, 2), "Pareto Front", align: center, fill:white, frame:"rect", padding:0.1, stroke: none)
  }),
  canvas(length: 1cm, {
    import draw: bezier, content
    plot-compare(filtered_data, a: 1, b: 0, c: 0, y-min: -2, y-max: 4, y-label: "Log10 Max Tensor Size")
    bezier((0, 3.3), (2, 0.7), (2, 3.2), stroke: (dash: "dashed"))
    content((3, 2.3), "Pareto Front", align: center, fill:white, frame:"rect", padding:0.1, stroke: none)
  })
)