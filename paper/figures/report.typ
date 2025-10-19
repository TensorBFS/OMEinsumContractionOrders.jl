#import "@preview/cetz:0.4.0": canvas, draw
#import "@preview/cetz-plot:0.1.2": plot
#show link: set text(blue)

// Read the merged summary data
#let all_data = json("data/summary.json")

// Function to extract instance name from path
#let get_instance_name(instance_path) = {
  let parts = instance_path.split("/")
  let filename = parts.last()
  filename.replace(".json", "")
}

// the cost function is a * sc + b * tc + c * rwc
#let plot-compare(dataset, x-max: auto, y-min: -3, y-max: 3.5, a: 1, b: 0, c: 0) = {
 plot.plot(
  size: (10, 7),
  x-label: "Log2 Contraction Cost",
  y-label: "Log10 Computing Time (seconds)",
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
        let optimizer_config = entry.optimizer_config
        
        if optimizer not in optimizer_groups {
          optimizer_groups.insert(optimizer, ())
        }
        optimizer_groups.at(optimizer).push((cost, time))
      }
      
      // Define colors and markers for different optimizers
      let colors = (red, blue, green, purple, orange, black)
      let markers = ("o", "x", "square", "triangle", "+")
      
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
          label: optimizer
        )
        i += 1
      }
    }
  )
}

// Group data by problem name + instance name from summary file
#let grouped_data = (:)

#for entry in all_data {
  let problem_name = entry.problem_name
  let instance_name = get_instance_name(entry.instance)
  let combined_key = problem_name + "/" + instance_name
  if combined_key not in grouped_data {
    grouped_data.insert(combined_key, ())
  }
  grouped_data.at(combined_key).push(entry)
}

#align(center, text(12pt)[= OMEinsumContractionOrders (v1.0.0) Benchmark Results])
#v(30pt)

Note: the `Treewidth` optimizer is greedy, only a subset of backends (`MF`, `MMD`, `AMF`) are tested (check #link("https://github.com/TensorBFS/OMEinsumContractionOrdersBenchmark/issues/2")[Issue #2]).

// Create individual plots for each problem + instance combination
#let combined_keys = grouped_data.keys().sorted()
#for combined_key in combined_keys {
  let dataset = grouped_data.at(combined_key)
  let a = 1
  let b = 0
  let c = 0
  figure(
    canvas(length: 1cm, {
      plot-compare(dataset, a: a, b: b, c: c)
    }),
    caption: [Scatter plot for *#combined_key* showing contraction cost (#a\*sc + #b\*tc + #c\*rwc) vs computing time for different optimizers.]
  )
}

// #pagebreak()

// // Summary statistics table
// Summary of benchmark results showing space complexity, computing time, and efficiency ratio for each instance and optimizer.
// #v(10pt)
// #table(
//   columns: 6,
//   stroke: 0.5pt,
//   [*Problem*], [*Instance*], [*Optimizer*], [*Space Complexity*], [*Computing Time (s)*], [*Efficiency*],
//   ..for combined_key in combined_keys {
//     let dataset = grouped_data.at(combined_key)
//     for entry in dataset {
//       let sc = entry.contraction_complexity.sc
//       let tc = entry.contraction_complexity.tc
//       let rwc = entry.contraction_complexity.rwc
//       let time = entry.time_elapsed
//       let efficiency = calc.round(sc / time, digits: 1)
//       let problem_name = entry.problem_name
//       let instance_name = get_instance_name(entry.instance)
//       let optimizer = entry.optimizer
//       (problem_name, instance_name, optimizer, str(sc), str(calc.round(time, digits: 4)), str(efficiency))
//     }
//   }
// )