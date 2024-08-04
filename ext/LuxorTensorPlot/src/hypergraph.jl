struct LabeledHyperGraph{TS, TV, TE}
    adjacency_matrix::SparseMatrixCSC{TS}
    vertex_labels::Vector{TV}
    edge_labels::Vector{TE}
    open_edges::Vector{TE}

    function LabeledHyperGraph(adjacency_matrix::SparseMatrixCSC{TS}; vl::Vector{TV} = [1:size(adjacency_matrix, 1)...], el::Vector{TE} = [1:size(adjacency_matrix, 2)...], oe::Vector = []) where{TS, TV, TE}
        if size(adjacency_matrix, 1) != length(vl)
            throw(ArgumentError("Number of vertices does not match number of vertex labels"))
        end
        if size(adjacency_matrix, 2) != length(el)
            throw(ArgumentError("Number of edges does not match number of edge labels"))
        end
        if !all(oei in el for oei in oe)
            throw(ArgumentError("Open edges must be in edge labels"))
        end
        if isempty(oe)
            oe = Vector{TE}()
        end
        new{TS, TV, TE}(adjacency_matrix, vl, el, oe)
    end
end

Base.show(io::IO, g::LabeledHyperGraph{TS, TV, TE}) where{TS,TV,TE} = print(io, "LabeledHyperGraph{$TS, $TV, $TE} \n adjacency_mat: $(g.adjacency_matrix) \n vertex: $(g.vertex_labels) \n edges: $(g.edge_labels)) \n open_edges: $(g.open_edges)")

Base.:(==)(a::LabeledHyperGraph, b::LabeledHyperGraph) = a.adjacency_matrix == b.adjacency_matrix && a.vertex_labels == b.vertex_labels && a.edge_labels == b.edge_labels && a.open_edges == b.open_edges

struct TensorNetworkGraph{TT, TI}
    graph::SimpleGraph
    tensors_labels::Dict{Int, TT}
    indices_labels::Dict{Int, TI}
    open_indices::Vector{TI}

    function TensorNetworkGraph(graph::SimpleGraph; tl::Dict{Int, TT} = Dict{Int, Int}(), il::Dict{Int, TI} = Dict{Int, Int}(), oi::Vector = []) where{TT, TI}
        if length(tl) + length(il) != nv(graph)
            throw(ArgumentError("Number of tensors + indices does not match number of vertices"))
        end
        if !all(oii in values(il) for oii in oi)
            throw(ArgumentError("Open indices must be in indices"))
        end
        if isempty(oi)
            oi = Vector{TI}()
        end
        new{TT, TI}(graph, tl, il, oi)
    end
end

Base.show(io::IO, g::TensorNetworkGraph{TT, TI}) where{TT, TI} = print(io, "TensorNetworkGraph{$TT, $TI} \n graph: {$(nv(g.graph)), $(ne(g.graph))} \n tensors: $(g.tensors_labels) \n indices: $(g.indices_labels)) \n open_indices: $(g.open_indices)")

# convert the labeled hypergraph to a tensor network graph, where vertices and edges of the hypergraph are mapped as the vertices of the tensor network graph, and the open edges are recorded.
function TensorNetworkGraph(lhg::LabeledHyperGraph{TS, TV, TE}) where{TS, TV, TE}
    graph = SimpleGraph(length(lhg.vertex_labels) + length(lhg.edge_labels))
    tensors_labels = Dict{Int, TV}()
    indices_labels = Dict{Int, TE}()

    lv = length(lhg.vertex_labels)
    for i in 1:length(lhg.vertex_labels)
        tensors_labels[i] = lhg.vertex_labels[i]
    end
    for i in 1:length(lhg.edge_labels)
        indices_labels[i + lv] = lhg.edge_labels[i]
    end

    for i in 1:size(lhg.adjacency_matrix, 1)
        for j in findall(!iszero, lhg.adjacency_matrix[i, :])
            add_edge!(graph, i, j + lv)
        end
    end
    
    TensorNetworkGraph(graph, tl=tensors_labels, il=indices_labels, oi=lhg.open_edges)
end