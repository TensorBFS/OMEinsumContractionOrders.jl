function LuxorGraphPlot.GraphViz(tng::TensorNetworkGraph, locs=StressLayout(); highlight::Vector=[], highlight_color = (0.0, 0.0, 255.0, 0.5), kwargs...)

    white = (255.0, 255.0, 255.0, 0.8)
    black = (0.0, 0.0, 0.0, 1.0)
    r = (255.0, 0.0, 0.0, 0.8)
    g = (0.0, 255.0, 0.0, 0.8)

    colors = Vector{typeof(r)}()
    text = Vector{String}()
    sizes = Vector{Float64}()

    for i in 1:nv(tng.graph)
        if i in keys(tng.tensors_labels)
            push!(colors, white)
            push!(text, string(tng.tensors_labels[i]))
            push!(sizes, 20.0)
        else
            push!(colors, r)
            push!(text, string(tng.indices_labels[i]))
            push!(sizes, 10.0)
        end
    end

    for oi in tng.open_indices
        id = _get_key(tng.indices_labels, oi)
        colors[id] = g
    end

    for hl in highlight
        id = _get_key(tng.indices_labels, hl)
        colors[id] = highlight_color
    end

    return GraphViz(tng.graph, locs, texts = text, vertex_colors = colors, vertex_sizes = sizes, kwargs...)
end

function _get_key(dict::Dict, value)
    for (key, val) in dict
        if val == value
            return key
        end
    end
    @error "Value not found in dictionary"
end

function OMEinsumContractionOrders.ein2hypergraph(ec::T) where{T <: AbstractEinsum}
    ixs = getixsv(ec)
    iy = getiyv(ec)

    edges = unique!([Iterators.flatten(ixs)...])
    open_edges = [iy[i] for i in 1:length(iy) if iy[i] in edges]

    rows = Int[]
    cols = Int[]
    for (i,ix) in enumerate(ixs)
        push!(rows, map(x->i, ix)...)
        push!(cols, map(x->findfirst(==(x), edges), ix)...)
    end
    adj = sparse(rows, cols, ones(Int, length(rows)))

    return LabeledHyperGraph(adj, el = edges, oe = open_edges)
end

"""
    viz_eins(ec::AbstractEinsum; locs=StressLayout(), filename = nothing, kwargs...)

Visualizes an `AbstractEinsum` object by creating a tensor network graph and rendering it using GraphViz.

## Arguments
- `ec::AbstractEinsum`: The `AbstractEinsum` object to visualize.
- `locs=StressLayout()`: The layout algorithm to use for positioning the nodes in the graph. Default is `StressLayout()`.
- `filename = nothing`: The name of the file to save the visualization to. If `nothing`, the visualization will be displayed on the screen instead of saving to a file.
- `kwargs...`: Additional keyword arguments to be passed to the `GraphViz` constructor.

"""
function OMEinsumContractionOrders.viz_eins(ec::AbstractEinsum; locs=StressLayout(), filename = nothing, kwargs...)
    tng = TensorNetworkGraph(ein2hypergraph(ec))
    gviz = GraphViz(tng, locs; kwargs...)
    return show_graph(gviz, filename = filename)
end