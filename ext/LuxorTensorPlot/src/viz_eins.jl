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

function OMEinsumContractionOrders.ein2hypergraph(code::T) where{T <: AbstractEinsum}
    ixs = getixsv(code)
    iy = getiyv(code)

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

function OMEinsumContractionOrders.viz_eins(code::AbstractEinsum; locs=StressLayout(), filename = nothing, config=LuxorTensorPlot.GraphDisplayConfig(), kwargs...)
    tng = TensorNetworkGraph(ein2hypergraph(code))
    gviz = GraphViz(tng, locs; kwargs...)
    return show_graph(gviz; filename, config)
end