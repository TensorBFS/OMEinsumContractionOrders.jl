function ein2hypergraph(args...; kwargs...)
    throw(ArgumentError("Extension `LuxorTensorPlot` not loaeded, please load it first by `using LuxorGraphPlot`."))
end

function ein2elimination(args...; kwargs...)
    throw(ArgumentError("Extension `LuxorTensorPlot` not loaeded, please load it first by `using LuxorGraphPlot`."))
end

"""
    viz_eins(code::AbstractEinsum; locs=StressLayout(), filename = nothing, kwargs...)

Visualizes an `AbstractEinsum` object by creating a tensor network graph and rendering it using GraphViz.

# Arguments
- `code::AbstractEinsum`: The `AbstractEinsum` object to visualize.

# Keyword Arguments
- `locs=StressLayout()`: The coordinates or layout algorithm to use for positioning the nodes in the graph.
- `filename = nothing`: The name of the file to save the visualization to. If `nothing`, the visualization will be displayed on the screen instead of saving to a file.
- `config = GraphDisplayConfig()`: The configuration for displaying the graph. Please refer to the documentation of [`GraphDisplayConfig`](https://giggleliu.github.io/LuxorGraphPlot.jl/dev/ref/#LuxorGraphPlot.GraphDisplayConfig) for more information.
- `kwargs...`: Additional keyword arguments to be passed to the [`GraphViz`](https://giggleliu.github.io/LuxorGraphPlot.jl/dev/ref/#LuxorGraphPlot.GraphViz) constructor.
"""
function viz_eins(args...; kwargs...)
    throw(ArgumentError("Extension `LuxorTensorPlot` not loaeded, please load it first by `using LuxorGraphPlot`."))
end

"""
    viz_contraction(code::Union{NestedEinsum, SlicedEinsum}; locs=StressLayout(), framerate=10, filename=tempname() * ".mp4", show_progress=true)

Visualize the contraction process of a tensor network.

# Arguments
- `code`: The tensor network to visualize.

# Keyword Arguments
- `locs`: The coordinates or layout algorithm to use for positioning the nodes in the graph. Default is `StressLayout()`.
- `framerate`: The frame rate of the animation. Default is `10`.
- `filename`: The name of the output file, with `.gif` or `.mp4` extension. Default is a temporary file with `.mp4` extension.
- `show_progress`: Whether to show progress information. Default is `true`.

# Returns
- the path of the generated file.
"""
function viz_contraction(args...; kwargs...)
    throw(ArgumentError("Extension `LuxorTensorPlot` not loaeded, please load it first by `using LuxorGraphPlot`."))
end