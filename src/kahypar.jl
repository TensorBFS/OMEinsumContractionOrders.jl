export uniformsize, optimize_kahypar

function uniformsize(@nospecialize(code::EinCode{ixs,iy}), size::Int) where {ixs, iy}
    Dict([c=>size for c in [Iterators.flatten(ixs)..., iy...]])
end
uniformsize(ne::NestedEinsum, size::Int) = uniformsize(Iterators.flatten(ne), size)

function induced_subhypergraph(s::SparseMatrixCSC, group)
    s0 = s[group,:]
    nvs = vec(sum(s0, dims=1))
    remaining_edges = findall(!iszero, nvs)
    s0[:,remaining_edges], remaining_edges
end

function kahypar_partitions_sc(adj::SparseMatrixCSC, vertices=collect(1:size(adj,1)); sc_target, log2_sizes, imbalances=0.02, verbose=false)
    n_v = length(vertices)
    subgraph, remaining_edges = induced_subhypergraph(adj, vertices)
    hypergraph = KaHyPar.HyperGraph(subgraph, ones(n_v), log2_sizes[remaining_edges])
    local parts
    for imbalance in imbalances
        parts = @suppress KaHyPar.partition(hypergraph, 2; imbalance=imbalance, configuration=:edge_cut)
        part0 = vertices[parts .== 0]
        part1 = vertices[parts .== 1]
        sc0, sc1 = group_sc(adj, part0), group_sc(adj, part1)
        sc = max(sc0, sc1)
        verbose && println("imbalance $imbalance: sc = $sc, group = ($(length(part0)), $(length(part1)))")
        if sc <= sc_target
            return part0, part1
        end
    end
    error("fail to find a valid partition for `sc_target = $sc_target`")
end

function group_sc(adj, group)
    degree_in = sum(adj[group,:], dims=1)
    degree_all = sum(adj, dims=1)
    count(i->degree_in[i]!=0 && degree_in[i]!=degree_all[i], 1:size(adj,2))
end

function kahypar_partitions_sc_recursive(adj::SparseMatrixCSC, vertices=collect(1:size(adj,1)); sc_target, max_size, log2_sizes, imbalances=0:0.03:0.99, verbose=false)
    if length(vertices) > max_size
        part1, part2 = kahypar_partitions_sc(adj, vertices; sc_target=sc_target, log2_sizes=log2_sizes, imbalances=imbalances, verbose=verbose)
        parts1 = kahypar_partitions_sc_recursive(adj, part1; sc_target=sc_target, max_size=max_size, log2_sizes=log2_sizes, imbalances=imbalances, verbose=verbose)
        parts2 = kahypar_partitions_sc_recursive(adj, part2; sc_target=sc_target, max_size=max_size, log2_sizes=log2_sizes, imbalances=imbalances, verbose=verbose)
        return [parts1, parts2]
    else
        return [vertices]
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

"""
    optimize_kahypar(code, size_dict; sc_target, max_group_size, imbalances=0.0:0.01:0.2, verbose=false)

Optimize the einsum code contraction order using the KaHyPar + Greedy approach.
This program first recursively cuts the tensors into several groups using KaHyPar,
with maximum group size specifed by `max_group_size` and maximum space complexity specified by `sc_target`,
Then finds the contraction order inside each group with the greedy search algorithm. Other arguments are

* `size_dict`, a dictionary that specifies leg dimensions,
* `imbalances`, KaHyPar parameter that controls the group sizes in hierarchical bipartition,
* `verbose`, showing more detail if true.

### References
* [Hyper-optimized tensor network contraction](https://arxiv.org/abs/2002.01935)
* [Simulating the Sycamore quantum supremacy circuits](https://arxiv.org/abs/2103.03074)
"""
function optimize_kahypar(@nospecialize(code::EinCode{ixs,iy}), size_dict; sc_target, max_group_size, imbalances=0.0:0.01:0.2, verbose=false) where {ixs, iy}
    ixv = collect(ixs)
    adj, edges = adjacency_matrix(ixv)
    vertices=collect(1:length(1:length(ixs)))
    parts = kahypar_partitions_sc_recursive(adj, vertices; sc_target, max_size=max_group_size, log2_sizes=[log2(size_dict[e]) for e in edges], imbalances, verbose=verbose)
    recursive_construct_nestedeinsum(ixv, collect(iy), parts, size_dict, 0)
end

function recursive_construct_nestedeinsum(ixs::AbstractVector, iy, parts::AbstractVector, size_dict, level=0)
    if length(parts) == 2
        # code is a nested einsum
        code1 = recursive_construct_nestedeinsum(ixs, iy, parts[1], size_dict, level+1)
        code2 = recursive_construct_nestedeinsum(ixs, iy, parts[2], size_dict, level+1)
        AB = recursive_flatten(parts[2]) ∪ recursive_flatten(parts[1])
        inset12, outset12 = ixs[AB], ixs[setdiff(1:length(ixs), AB)]
        iy12 = Iterators.flatten(inset12) ∩  (Iterators.flatten(outset12) ∪ iy)
        iy1, iy2 = OMEinsum.getiy(code1.eins), OMEinsum.getiy(code2.eins)
        return NestedEinsum((code1, code2), EinCode((iy1, iy2), ((level==0 ? iy : iy12)...,)))
    elseif length(parts) == 1
        return recursive_construct_nestedeinsum(ixs, iy, parts[1], size_dict, level)
    else
        error("not a bipartition, got size $(length(parts))")
    end
end

function recursive_construct_nestedeinsum(ixs::AbstractVector, iy, parts::AbstractVector{<:Integer}, size_dict, level=0)
    if isempty(parts)
        error("got empty group!")
    end
    inset, outset = ixs[parts], ixs[setdiff(1:length(ixs), parts)]
    iy1 = Iterators.flatten(inset) ∩  (Iterators.flatten(outset) ∪ iy)
    code = EinCode{(inset...,), (iy1...,)}()
    res = optimize_greedy(code, size_dict)
    return maplocs(res, parts)
end

maplocs(ne::NestedEinsum, parts) = NestedEinsum(maplocs.(ne.args, Ref(parts)), ne.eins)
maplocs(i::Int, parts) = parts[i]

function kahypar_recursive(ne::NestedEinsum; log2_size_dict, sc_target, min_size, imbalances=0.0:0.04:0.8)
    if length(ne.args >= min_size) && all(x->x isa Integer, ne.args)
        bipartite_eincode(adj, ne.args, ne.eins; log2_size_dict=log2_size_dict, sc_target=sc_target, min_size=min_size, imbalances=imbalances)
    end
    kahypar_recursive(ne.args; log2_size_dict, sc_target=sc_target, min_size=min_size, imbalances=imbalances)
end

recursive_flatten(obj::Tuple) = vcat(recursive_flatten.(obj)...)
recursive_flatten(obj::AbstractVector) = vcat(recursive_flatten.(obj)...)
recursive_flatten(obj) = obj
