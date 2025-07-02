#import "@preview/cetz:0.2.2": canvas, draw, tree, plot
#import "@preview/ctheorems:1.1.3": *
#import "@preview/ouset:0.2.0": ouset

#let tensor(location, name, label, color: black) = {
  import draw: *
  circle(location, radius: 10pt, name: name)
  content((), text(color, label))
}


#let labelnode(loc, label, name: none) = {
  import draw: *
  content(loc, text(black, label), align: center, fill:silver, frame:"rect", padding:0.07, stroke: none, name: name)
}
#let labeledge(from, to, label, name: none) = {
  import draw: *
  line(from, to, name: "line")
  labelnode("line.mid", label, name: name)
}

#set page(width: auto, height: auto, margin: 5pt)