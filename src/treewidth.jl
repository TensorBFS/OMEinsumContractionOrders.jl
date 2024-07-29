using TreeWidthSolver
using Graphs: connected_components, induced_subgraph, SimpleGraph, add_edge!

"""
    struct ExactTreewidth{GM} <: CodeOptimizer
    ExactTreewidth(greedy_config::GM = GreedyMethod(nrepeat=1))

A optimizer using the exact tree width solver proved in TreeWidthSolver.jl, the greedy_config is the configuration for the greedy method, which is used to solve the subproblems in the tree decomposition.

# Fields
- `greedy_config::GM`: The configuration for the greedy method.

"""
Base.@kwdef struct ExactTreewidth{GM} <: CodeOptimizer 
    greedy_config::GM = GreedyMethod(nrepeat=1)
end

"""
    exact_treewidth_method(incidence_list::IncidenceList{VT,ET}, log2_edge_sizes; α::TA = 0.0, temperature::TT = 0.0, nrepeat=1) where {VT,ET,TA,TT}

This function calculates the exact treewidth of a graph using TreeWidthSolver.jl. It takes an incidence list representation of the graph (`incidence_list`) and a dictionary of logarithm base 2 edge sizes (`log2_edge_sizes`) as input. The function also accepts optional parameters `α`, `temperature`, and `nrepeat` with default values of 0.0, 0.0, and 1 respectively, which are parameter of the GreedyMethod used in the contraction process.

## Arguments
- `incidence_list`: An incidence list representation of the graph.
- `log2_edge_sizes`: A dictionary of logarithm base 2 edge sizes.

## Returns
- The function returns a `ContractionTree` representing the contraction process.

```
julia> optimizer = OMEinsumContractionOrders.ExactTreewidth()
OMEinsumContractionOrders.ExactTreewidth{GreedyMethod{Float64, Float64}}(GreedyMethod{Float64, Float64}(0.0, 0.0, 1))

julia> eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], ['a'])
ab, acd, bcef, e, df -> a

julia> size_dict = Dict([c=>(1<<i) for (i,c) in enumerate(['a', 'b', 'c', 'd', 'e', 'f'])]...)
Dict{Char, Int64} with 6 entries:
  'f' => 64
  'a' => 2
  'c' => 8
  'd' => 16
  'e' => 32
  'b' => 4

julia> optcode = optimize_code(eincode, size_dict, optimizer)
ab, ab -> a
├─ ab
└─ fac, bcf -> ab
   ├─ df, acd -> fac
   │  ├─ df
   │  └─ acd
   └─ e, bcef -> bcf
      ├─ e
      └─ bcef
```

"""
function exact_treewidth_method(incidence_list::IncidenceList{VT,ET}, log2_edge_sizes; α::TA = 0.0, temperature::TT = 0.0, nrepeat=1) where {VT,ET,TA,TT}
    indicies = collect(keys(incidence_list.e2v))
    weights = [log2_edge_sizes[e] for e in indicies]
    line_graph = il2lg(incidence_list, indicies)

    contraction_trees = Vector{Union{ContractionTree, VT}}()

    # avoid the case that the line graph is not connected
    for vertice_ids in connected_components(line_graph)
        lg = induced_subgraph(line_graph, vertice_ids)[1]
        lg_indicies = indicies[vertice_ids]
        lg_weights = weights[vertice_ids]
        labeled_graph = LabeledSimpleGraph(lg, lg_indicies, lg_weights)
        tree_decomposition = exact_treewidth(labeled_graph)
        elimination_order = EliminationOrder(tree_decomposition.tree)
        push!(contraction_trees, eo2ct(elimination_order, incidence_list, log2_edge_sizes, α, temperature, nrepeat))
    end
    
    if length(contraction_trees) == 1 # the line graph having only one connected component
        return contraction_trees[1]
    else # the line graph having multiple connected components
        contraction_tree = contraction_trees[1]
        for i in 2:length(contraction_trees)
            contraction_tree = ContractionTree(contraction_tree, contraction_trees[i])
        end
        return contraction_tree
    end
end

# transform incidence list to line graph
function il2lg(incidence_list::IncidenceList{VT, ET}, indicies::Vector{ET}) where {VT, ET}

    line_graph = SimpleGraph(length(indicies))
    
    for (i, e) in enumerate(indicies)
        for v in incidence_list.e2v[e]
            for ej in incidence_list.v2e[v]
                if e != ej add_edge!(line_graph, i, findfirst(==(ej), indicies)) end
            end
        end
    end

    return line_graph
end

# transform elimination order to contraction tree
function eo2ct(elimination_order::EliminationOrder, incidence_list::IncidenceList{VT, ET}, log2_edge_sizes, α::TA, temperature::TT, nrepeat) where {VT, ET, TA, TT}
    eo = copy(elimination_order.order)
    incidence_list = copy(incidence_list)
    contraction_tree_nodes = Vector{Union{VT, ContractionTree}}(collect(keys(incidence_list.v2e)))
    tensors_list = Dict{VT, Int}()
    for (i, v) in enumerate(contraction_tree_nodes)
        tensors_list[v] = i
    end

    flag = contraction_tree_nodes[1]

    while !isempty(eo)
        e = pop!(eo)
        if haskey(incidence_list.e2v, e)
            vs = incidence_list.e2v[e]
            if length(vs) == 2
                vi = vs[1]
                vj = vs[2]
                contract_pair!(incidence_list, vi, vj, log2_edge_sizes)
                node_i = contraction_tree_nodes[tensors_list[vi]]
                node_j = contraction_tree_nodes[tensors_list[vj]]
                contraction_tree_nodes[tensors_list[vi]] = ContractionTree(node_i, node_j)
                flag = vi
            elseif length(vs) > 2
                sub_list = IncidenceList(Dict([v => incidence_list.v2e[v] for v in vs]); openedges=incidence_list.openedges)
                sub_tree, scs, tcs = tree_greedy(sub_list, log2_edge_sizes; nrepeat=nrepeat, α=α, temperature=temperature)
                vi = contract_tree!(incidence_list, sub_tree, log2_edge_sizes, scs, tcs)
                contraction_tree_nodes[tensors_list[vi]] = sub_tree
                flag = vi
            end
        end
    end

    return contraction_tree_nodes[tensors_list[flag]]
end

"""
    optimize_greedy(optimizer, eincode, size_dict)

Optimizing the contraction order via solve the exact tree width of the line graph corresponding to the eincode and return a `NestedEinsum` object.
Check the docstring of `exact_treewidth_method` for detailed explaination of other input arguments.
"""
function optimize_exact_treewidth(optimizer::ExactTreewidth{GM}, code::EinCode{L}, size_dict::Dict) where {L,GM}
    optimize_exact_treewidth(optimizer, getixsv(code), getiyv(code), size_dict)
end
function optimize_exact_treewidth(optimizer::ExactTreewidth{GM}, ixs::AbstractVector{<:AbstractVector}, iy::AbstractVector, size_dict::Dict{L,TI}) where {L, TI, GM}
    if length(ixs) <= 2
        return NestedEinsum(NestedEinsum{L}.(1:length(ixs)), EinCode(ixs, iy))
    end
    log2_edge_sizes = Dict{L,Float64}()
    for (k, v) in size_dict
        log2_edge_sizes[k] = log2(v)
    end
    incidence_list = IncidenceList(Dict([i=>ixs[i] for i=1:length(ixs)]); openedges=iy)
    α = optimizer.greedy_config.α
    temperature = optimizer.greedy_config.temperature
    nrepeat = optimizer.greedy_config.nrepeat
    tree = exact_treewidth_method(incidence_list, log2_edge_sizes; α = α, temperature = temperature, nrepeat=nrepeat)
    parse_eincode!(incidence_list, tree, 1:length(ixs))[2]
end