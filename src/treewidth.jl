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

This function calculates the exact treewidth of a graph using TreeWidthSolver.jl. It takes an incidence list representation of the graph (`incidence_list`) and a dictionary of logarithm base 2 edge sizes (`log2_edge_sizes`) as input. The function also accepts optional parameters `α`, `temperature`, and `nrepeat` with default values of 0.0, 0.0, and 1 respectively, which are parameter of the GreedyMethod used in the contraction process as a sub optimizer.

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
    tensors = collect(keys(incidence_list.v2e))
    weights = [log2_edge_sizes[e] for e in indicies]
    line_graph = il2lg(incidence_list, indicies)

    scalars = [i for i in tensors if isempty(incidence_list.v2e[i])]
    contraction_trees = Vector{Union{ContractionTree, VT}}()

    # avoid the case that the line graph is not connected
    for vertice_ids in connected_components(line_graph)
        lg = induced_subgraph(line_graph, vertice_ids)[1]
        lg_indicies = indicies[vertice_ids]
        lg_weights = weights[vertice_ids]
        eo = elimination_order(lg, labels = lg_indicies, weights = lg_weights)
        contraction_tree = eo2ct(eo, incidence_list, log2_edge_sizes, α, temperature, nrepeat)
        push!(contraction_trees, contraction_tree)
    end

    # add the scalars back to the contraction tree
    return reduce((x,y) -> ContractionTree(x, y), contraction_trees ∪ scalars)
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
function eo2ct(elimination_order::Vector{Vector{TL}}, incidence_list::IncidenceList{VT, ET}, log2_edge_sizes, α::TA, temperature::TT, nrepeat) where {TL, VT, ET, TA, TT}
    eo = copy(elimination_order)
    incidence_list = copy(incidence_list)
    contraction_tree_nodes = Vector{Union{VT, ContractionTree}}(collect(keys(incidence_list.v2e)))
    tensors_list = Dict{VT, Int}()
    for (i, v) in enumerate(contraction_tree_nodes)
        tensors_list[v] = i
    end

    flag = contraction_tree_nodes[1]

    while !isempty(eo)
        eliminated_vertices = pop!(eo) # e is a vector of vertices, which are eliminated at the same time
        vs = unique!(vcat([incidence_list.e2v[ei] for ei in eliminated_vertices if haskey(incidence_list.e2v, ei)]...)) # the tensors to be contracted, since they are connected to the eliminated vertices
        if length(vs) >= 2
            sub_list_indices = unique!(vcat([incidence_list.v2e[v] for v in vs]...)) # the vertices connected to the tensors to be contracted
            sub_list_open_indices = setdiff(sub_list_indices, eliminated_vertices) # the vertices connected to the tensors to be contracted but not eliminated
            sub_list = IncidenceList(Dict([v => incidence_list.v2e[v] for v in vs]); openedges=sub_list_open_indices) # the subgraph of the contracted tensors
            sub_tree, scs, tcs = tree_greedy(sub_list, log2_edge_sizes; nrepeat=nrepeat, α=α, temperature=temperature) # optmize the subgraph with greedy method
            vi = contract_tree!(incidence_list, sub_tree, log2_edge_sizes, scs, tcs) # insert the contracted tensors back to the total graph
            contraction_tree_nodes[tensors_list[vi]] = st2ct(sub_tree, tensors_list, contraction_tree_nodes)
            flag = vi
        end
    end

    return contraction_tree_nodes[tensors_list[flag]]
end

function st2ct(sub_tree::Union{ContractionTree, VT}, tensors_list::Dict{VT, Int}, contraction_tree_nodes::Vector{Union{ContractionTree, VT}}) where{VT}
    if sub_tree isa ContractionTree
        return ContractionTree(st2ct(sub_tree.left, tensors_list, contraction_tree_nodes), st2ct(sub_tree.right, tensors_list, contraction_tree_nodes))
    else
        return contraction_tree_nodes[tensors_list[sub_tree]]        
    end
end

"""
    optimize_exact_treewidth(optimizer, eincode, size_dict)

Optimizing the contraction order via solve the exact tree width of the line graph corresponding to the eincode and return a `NestedEinsum` object.
Check the docstring of `exact_treewidth_method` for detailed explaination of other input arguments.
"""
function optimize_exact_treewidth(optimizer::ExactTreewidth{GM}, code::AbstractEinsum, size_dict::Dict) where {GM}
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
    # complete all open edges as a clique, connected with a dummy tensor
    incidence_list = IncidenceList(Dict([i=>ixs[i] for i=1:length(ixs)] ∪ [(length(ixs) + 1 => iy)]))

    α = optimizer.greedy_config.α
    temperature = optimizer.greedy_config.temperature
    nrepeat = optimizer.greedy_config.nrepeat
    tree = exact_treewidth_method(incidence_list, log2_edge_sizes; α = α, temperature = temperature, nrepeat=nrepeat)

    # remove the dummy tensor added for open edges
    optcode = parse_eincode!(incidence_list, tree, 1:length(ixs) + 1)[2]

    return pivot_tree(optcode, length(ixs) + 1)
end