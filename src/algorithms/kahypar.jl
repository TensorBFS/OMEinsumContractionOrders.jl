export optimize_kahypar, optimize_kahypar_auto
export KaHyParBipartite

"""
    KaHyParBipartite{RT,IT,GM}
    KaHyParBipartite(; sc_target, imbalances=collect(0.0:0.005:0.8),
        max_group_size=40, greedy_config=GreedyMethod())

Optimize the einsum code contraction order using the KaHyPar + Greedy approach.
This program first recursively cuts the tensors into several groups using KaHyPar,
with maximum group size specifed by `max_group_size` and maximum space complexity specified by `sc_target`,
Then finds the contraction order inside each group with the greedy search algorithm. Other arguments are

* `sc_target` is the target space complexity, defined as `log2(number of elements in the largest tensor)`,
* `imbalances` is a KaHyPar parameter that controls the group sizes in hierarchical bipartition,
* `max_group_size` is the maximum size that allowed to used greedy search,
* `greedy_config` is a greedy optimizer.

### References
* [Hyper-optimized tensor network contraction](https://arxiv.org/abs/2002.01935)
* [Simulating the Sycamore quantum supremacy circuits](https://arxiv.org/abs/2103.03074)
"""
Base.@kwdef struct KaHyParBipartite{RT,IT,GM} <: CodeOptimizer
    sc_target::RT
    imbalances::IT = 0.0:0.005:0.8
    max_group_size::Int = 40
    greedy_config::GM = GreedyMethod()
end

function induced_subhypergraph(s::SparseMatrixCSC, group)
    s0 = s[group,:]
    nvs = vec(sum(s0, dims=1))
    remaining_edges = findall(!iszero, nvs)
    s0[:,remaining_edges], remaining_edges
end

function convert2int(sizes::AbstractVector)
    round.(Int, sizes .* 100)
end

function bipartite_sc(bipartiter::KaHyParBipartite, adj::SparseMatrixCSC, vertices, log2_sizes)
    n_v = length(vertices)
    subgraph, remaining_edges = induced_subhypergraph(adj, vertices)
    if !isdefined(@__MODULE__, :KaHyPar)
        error("Module `KaHyPar` not found, please type `using KaHyPar` before using the `KaHyParBipartite` optimizer!")
    end
    hypergraph = KaHyPar.HyperGraph(subgraph, ones(n_v), convert2int(log2_sizes[remaining_edges]))
    local parts
    min_sc = 999999
    for imbalance in bipartiter.imbalances
        parts = @suppress KaHyPar.partition(hypergraph, 2; imbalance=imbalance, configuration=:edge_cut)
        part0 = vertices[parts .== 0]
        part1 = vertices[parts .== 1]
        sc0, sc1 = group_sc(adj, part0, log2_sizes), group_sc(adj, part1, log2_sizes)
        sc = max(sc0, sc1)
        min_sc = min(sc, min_sc)
        @debug "imbalance $imbalance: sc = $sc, group = ($(length(part0)), $(length(part1)))"
        if sc <= bipartiter.sc_target
            return part0, part1
        end
    end
    error("fail to find a valid partition for `sc_target = $(bipartiter.sc_target)`, got minimum value `$min_sc` (imbalances = $(bipartiter.imbalances))")
end

# the space complexity (external degree of freedoms) if we contract this group
function group_sc(adj, group, log2_sizes)
    degree_in = sum(adj[group,:], dims=1)
    degree_all = sum(adj, dims=1)
    sum(i->(degree_in[i]!=0 && degree_in[i]!=degree_all[i] ? Float64(log2_sizes[i]) : 0.0), 1:size(adj,2))
end

function bipartition_recursive(bipartiter, adj::SparseMatrixCSC, vertices::AbstractVector{T}, log2_sizes) where T
    if length(vertices) > bipartiter.max_group_size
        parts = bipartite_sc(bipartiter, adj, vertices, log2_sizes)
        groups = Vector{T}[]
        for part in parts
            for component in _connected_components(adj, part)
                push!(groups, component)
            end
        end
        newparts = [bipartition_recursive(bipartiter, adj, groups[i], log2_sizes) for i=1:length(groups)]
        if length(groups) > 2
            tree = coarse_grained_optimize(adj, groups, log2_sizes, bipartiter.greedy_config)
            return map_tree_to_parts(tree, newparts)
        else
            return newparts
        end
    else
        return [vertices]
    end
end

function _connected_components(adj, part::AbstractVector{T}) where T
    A = adj[part,:]
    A = A * A'   # connectivility matrix
    n = length(part)
    visit_mask = zeros(Bool, n)
    groups = Vector{T}[]
    while !all(visit_mask)
        newset = Int[]
        push_connected!(newset, visit_mask, A, findfirst(==(false), visit_mask))
        push!(groups, getindex.(Ref(part), newset))
    end
    return groups
end

function push_connected!(set, visit_mask, adj, i)
    visit_mask[i] = true
    push!(set, i)
    for v = 1:size(adj, 2)
        if !visit_mask[v] && !iszero(adj[i,v])
            push_connected!(set, visit_mask, adj, v)
        end
    end
end

function coarse_grained_optimize(adj, parts, log2_sizes, greedy_config)
    incidence_list = get_coarse_grained_graph(adj, parts)
    log2_edge_sizes = Dict([i=>log2_sizes[i] for i=1:length(log2_sizes)])
    tree, _, _ = tree_greedy(incidence_list, log2_edge_sizes; method=greedy_config.method, nrepeat=greedy_config.nrepeat)
    return tree
end

function get_coarse_grained_graph(adj, parts)
    ADJ = vcat([sum(adj[part,:], dims=1) for part in parts]...)
    degree_in = sum(ADJ, dims=1)
    degree_all = sum(adj, dims=1)
    openedges = filter(i->degree_in[i]!=0 && degree_in[i]!=degree_all[i], 1:size(adj, 2))
    v2e = Dict{Int,Vector{Int}}()
    for v=1:size(ADJ, 1)
        v2e[v] = findall(!iszero, view(ADJ,v,:))
    end
    incidence_list = IncidenceList(v2e; openedges=openedges)
    return incidence_list
end

function map_tree_to_parts(tree, parts)
    if tree isa ContractionTree
        [map_tree_to_parts(tree.left, parts), map_tree_to_parts(tree.right, parts)]
    else
        parts[tree]
    end
end

# KaHyPar
function adjacency_matrix(ixs::AbstractVector)
    rows = Int[]
    cols = Int[]
    edges = unique!([Iterators.flatten(ixs)...])
    for (i,ix) in enumerate(ixs)
        push!(rows, map(x->i, ix)...)
        push!(cols, map(x->findfirst(==(x), edges), ix)...)
    end
    return sparse(rows, cols, ones(Int, length(rows))), edges
end

# legacy interface
"""
    optimize_kahypar(code, size_dict; sc_target, max_group_size=40, imbalances=0.0:0.01:0.2, greedy_method=MinSpaceOut(), greedy_nrepeat=10)

Optimize the einsum `code` contraction order using the KaHyPar + Greedy approach. `size_dict` is a dictionary that specifies leg dimensions. 
Check the docstring of `KaHyParBipartite` for detailed explaination of other input arguments.
"""
function optimize_kahypar(code::EinCode, size_dict; sc_target, max_group_size=40, imbalances=0.0:0.01:0.2, greedy_method=MinSpaceOut(), greedy_nrepeat=10)
    bipartiter = KaHyParBipartite(; sc_target=sc_target, max_group_size=max_group_size, imbalances=imbalances, greedy_config=GreedyMethod(method=greedy_method, nrepeat=greedy_nrepeat))
    recursive_bipartite_optimize(bipartiter, code, size_dict)
end

function recursive_bipartite_optimize(bipartiter, code::EinCode, size_dict)
    ixs, iy = getixsv(code), getiyv(code)
    ixv = [ixs..., iy]
    adj, edges = adjacency_matrix(ixv)
    vertices=collect(1:length(ixs))
    parts = bipartition_recursive(bipartiter, adj, vertices, [log2(size_dict[e]) for e in edges])
    recursive_construct_nestedeinsum(ixv, iy, parts, size_dict, 0, bipartiter.greedy_config)
end

maplocs(ne::NestedEinsum{ET}, parts) where ET = isleaf(ne) ? NestedEinsum{ET}(parts[ne.tensorindex]) : NestedEinsum(maplocs.(ne.args, Ref(parts)), ne.eins)

function kahypar_recursive(ne::NestedEinsum; log2_size_dict, sc_target, min_size, imbalances=0.0:0.04:0.8)
    if length(ne.args >= min_size) && all(isleaf, ne.args)
        bipartite_eincode(adj, ne.args, ne.eins; log2_size_dict=log2_size_dict, sc_target=sc_target, min_size=min_size, imbalances=imbalances)
    end
    kahypar_recursive(ne.args; log2_size_dict, sc_target=sc_target, min_size=min_size, imbalances=imbalances)
end

recursive_flatten(obj::Tuple) = vcat(recursive_flatten.(obj)...)
recursive_flatten(obj::AbstractVector) = vcat(recursive_flatten.(obj)...)
recursive_flatten(obj) = obj

"""
    optimize_kahypar_auto(code, size_dict; max_group_size=40, greedy_method=MinSpaceOut(), greedy_nrepeat=10)

Find the optimal contraction order automatically by determining the `sc_target` with bisection.
It can fail if the tree width of your graph is larger than `100`.
"""
function optimize_kahypar_auto(code::EinCode, size_dict; max_group_size=40, effort=500, greedy_method=MinSpaceOut(), greedy_nrepeat=10)
    sc_high = 100
    sc_low = 1
    order_high = optimize_kahypar(code, size_dict; sc_target=sc_high, max_group_size=max_group_size, imbalances=0.0:0.6/effort*(sc_high-sc_low):0.6)
    _optimize_kahypar_auto(code, size_dict, sc_high, order_high, sc_low, max_group_size, effort, greedy_method, greedy_nrepeat)
end
function _optimize_kahypar_auto(code::EinCode, size_dict, sc_high, order_high, sc_low, max_group_size, effort, greedy_method, greedy_nrepeat)
    if sc_high <= sc_low + 1
        order_high
    else
        sc_mid = (sc_high + sc_low) ÷ 2
        try
            order_mid = optimize_kahypar(code, size_dict; sc_target=sc_mid, max_group_size=max_group_size, imbalances=0.0:0.6/effort*(sc_high-sc_low):0.6)
            order_high, sc_high = order_mid, sc_mid
            # `sc_target` too high
        catch
            # `sc_target` too low
            sc_low = sc_mid
        end
        _optimize_kahypar_auto(code, size_dict, sc_high, order_high, sc_low, max_group_size, effort, greedy_method, greedy_nrepeat)
    end
end

function recursive_construct_nestedeinsum(ixs::AbstractVector{<:AbstractVector}, iy::AbstractVector{L}, parts::AbstractVector, size_dict, level, greedy_config) where L
    if length(parts) == 2
        # code is a nested einsum
        code1 = recursive_construct_nestedeinsum(ixs, iy, parts[1], size_dict, level+1, greedy_config)
        code2 = recursive_construct_nestedeinsum(ixs, iy, parts[2], size_dict, level+1, greedy_config)
        AB = recursive_flatten(parts[2]) ∪ recursive_flatten(parts[1])
        inset12, outset12 = ixs[AB], ixs[setdiff(1:length(ixs), AB)]
        iy12 = Iterators.flatten(inset12) ∩  (Iterators.flatten(outset12) ∪ iy)
        iy1, iy2 = getiyv(code1.eins), getiyv(code2.eins)
        return NestedEinsum([code1, code2], EinCode([iy1, iy2], L[(level==0 ? iy : iy12)...]))
    elseif length(parts) == 1
        return recursive_construct_nestedeinsum(ixs, iy, parts[1], size_dict, level, greedy_config)
    else
        error("not a bipartition, got size $(length(parts))")
    end
end

function recursive_construct_nestedeinsum(ixs::AbstractVector{<:AbstractVector}, iy::AbstractVector{L}, parts::AbstractVector{<:Integer}, size_dict, level, greedy_config) where L
    if isempty(parts)
        error("got empty group!")
    end
    inset, outset = ixs[parts], ixs[setdiff(1:length(ixs), parts)]
    iy1 = level == 0 ? iy : Iterators.flatten(inset) ∩  (Iterators.flatten(outset) ∪ iy)
    res = optimize_greedy(inset, iy1, size_dict; method=greedy_config.method, nrepeat=greedy_config.nrepeat)
    return maplocs(res, parts)
end
