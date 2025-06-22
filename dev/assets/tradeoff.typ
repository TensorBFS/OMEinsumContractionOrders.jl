#import "@preview/cetz:0.2.2": canvas, draw, tree, plot
#set page(width: auto, height: auto, margin: 5pt)

#figure(canvas(length:0.9cm, {
  import plot
  import draw: *
  let s(it) = text(10pt, it)
  hobby((1, 7.5), (1, 5), (3, 5), (3, 7), fill: blue.lighten(70%), close: true, stroke: none)
  plot.plot(size: (12,8),
    x-tick-step: none,
    y-tick-step: none,
    x-label: text(13pt)[Time to optimize contraction order],
    y-label: text(13pt)[Time to contract],
    y-max: 10,
    y-min: -2,
    x-max: 10,
    x-min: 0,
    name: "plot",
    {
      let greedy = (1, 9)
      let localsearch = (6, 1)
      let bipartition = (4.2, 3.0)
      let sabipartite = (5.5, 3.7)
      let exacttreewidth = (8.5, 0)
      let hypernd = (2, 3.0)
      let treewidth = (2.5, 8.0)
      plot.add(
        (greedy, hypernd, localsearch, exacttreewidth), style: (stroke: (paint: black, dash: "dashed")),
      )
      plot.add-anchor("greedy", greedy)
      plot.add(
        (greedy, bipartition, localsearch, hypernd, exacttreewidth, sabipartite), style: (stroke: none), mark:"o",
      )
      plot.add-anchor("greedy", greedy)
      plot.add-anchor("localsearch", localsearch)
      plot.add-anchor("bipartition", bipartition)
      plot.add-anchor("sabipartite", sabipartite)
      plot.add-anchor("exacttreewidth", exacttreewidth)
      plot.add-anchor("hypernd", hypernd)
      plot.add-anchor("treewidth", treewidth)
    }
  )
  content((rel: (0.0, -0.3), to: "plot.greedy"), s[`GreedyMethod`])
  content((rel: (0.0, -0.3), to: "plot.localsearch"), s[`TreeSA`])
  content((rel: (0.0, -0.3), to: "plot.bipartition"), box(fill: white, inset: 1pt, s[`KaHyParBipartite`]))
  content((rel: (0.0, -0.3), to: "plot.sabipartite"), s[`SABipartite`])
  content((rel: (0, -0.3), to: "plot.exacttreewidth"), box(fill: white, inset: 1pt, s[`ExactTreewidth`]))
  content((rel: (0, -0.3), to: "plot.hypernd"), box(fill: white, s[`HyperND`], inset: 1pt))
  content((rel: (0, -0.3), to: "plot.treewidth"), s[`Treewidth`])
}))


