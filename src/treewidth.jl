"""
    struct Treewidth{EL <: EliminationAlgorithm, GM} <: CodeOptimizer
    Treewidth(; alg::EL = SafeRules(BT(), MMW{3}(), MF()))

Tree width based solver. The solvers are implemented in [CliqueTrees.jl](https://algebraicjulia.github.io/CliqueTrees.jl/stable/) and [TreeWidthSolver.jl](https://github.com/ArrogantGao/TreeWidthSolver.jl). They include:

| Algorithm | Description | Time Complexity | Space Complexity |
|:-----------|:-------------|:----------------|:-----------------|
| `AMF` | approximate minimum fill | O(mn) | O(m + n) |
| `MF` | minimum fill | O(mn²) | - |
| `MMD` | multiple minimum degree | O(mn²) | O(m + n) |

Detailed descriptions is available in the [CliqueTrees.jl](https://algebraicjulia.github.io/CliqueTrees.jl/stable/api/#Elimination-Algorithms).

# Fields
- `alg::EL`: The algorithm to use for the treewidth calculation. Available elimination algorithms are listed above.

# Example
```jldoctest
julia> optimizer = Treewidth();

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
ab, ba -> a
├─ ab
└─ bcf, acf -> ba
   ├─ bcef, e -> bcf
   │  ├─ bcef
   │  └─ e
   └─ acd, df -> acf
      ├─ acd
      └─ df
```
"""
Base.@kwdef struct Treewidth{EL <: EliminationAlgorithm} <: CodeOptimizer 
    alg::EL = SafeRules(BT(), MMW{3}(), MF())
end

"""
    const ExactTreewidth = Treewidth{SafeRules{BT, MMW{3}(), MF}}
    ExactTreewidth() = Treewidth()

`ExactTreewidth` is a specialization of `Treewidth` for the `SafeRules` preprocessing algorithm with the `BT` elimination algorithm.
The `BT` algorithm is an exact solver for the treewidth problem that implemented in [`TreeWidthSolver.jl`](https://github.com/ArrogantGao/TreeWidthSolver.jl).
"""
const ExactTreewidth = Treewidth{SafeRules{BT, MMW{3}, MF}}
ExactTreewidth() = Treewidth()

"""
    optimize_treewidth(optimizer, eincode, size_dict)

Optimizing the contraction order via solve the exact tree width of the line graph corresponding to the eincode and return a `NestedEinsum` object.
Check the docstring of `treewidth_method` for detailed explaination of other input arguments.
"""
function optimize_treewidth(optimizer::Treewidth, code::AbstractEinsum, size_dict::Dict; binary::Bool=true)
    optimize_treewidth(optimizer, getixsv(code), getiyv(code), size_dict; binary)
end

function optimize_treewidth(optimizer::Treewidth, ixs::AbstractVector{<:AbstractVector}, iy::AbstractVector, size_dict::Dict{L, Int}; binary::Bool=true) where {L}
    log2_size_dict = _log2_size_dict(size_dict)
    marker = zeros(Int, max(length(ixs) + 1, length(size_dict)))

    # construct incidence matrix `ve`
    #           indices
    #         [         ]
    # tensors [    ve   ]
    #         [         ]
    # we only care about the sparsity pattern
    weights, ev, ve, el = einexpr_to_matrix!(marker, ixs, iy, size_dict)

    # compute a tree (forest) decomposition of `ve`
    tree = matrix_to_tree!(marker, weights, ev, ve, el, optimizer.alg)

    # transform tree decomposition in contraction tree
    code = tree_to_einexpr!(marker, tree, ve, el, ixs, iy)

    if binary
        # binarize contraction tree
        code = optimize_greedy_log2size(code, log2_size_dict; α = 0.0, temperature = 0.0)
    end

    return code
end

"""
    einexpr_to_matrix!(marker, ixs, iy, size_dict)

Construct the weighted incidence matrix correponding to an Einstein summation expression.
Returns a quadruple (weights, ev, ve, el).

Each Einstein summation expression has a set E ⊆ L of indices, a set V := {1, …, |V|} of
(inner) tensors, and an outer tensor * := |V| + 1. Each tensor v ∈ V is incident to a sequence
ixs[v] of indices, and the outer tensor is incident to the sequence iy. Note that an index
can appear multiple times in ixs[v], e.g.

    ixs[v] = ('a', 'a', 'b').

Each index l ∈ E also has a positive dimension, given by size_dict[l].

The function `einexpr_to_matrix` does two things. First of all, it enumerates the index set
E, mapping each index to a distinct natural number.

    el: {1, …, |E|} → E
    le: E → {1, …, |E|}

Next, it constructs a vector weights: {1, …, |E|} → [0, ∞) satisfying

    weights[e] := log2(size_dict[el[e]]),

and a sparse matrix ve: {1, …, |V| + 1} × {1, …, |E|} → {0, 1} satisfying

    ve[v, e] := { 1 if el[e] is incident to v
                { 0 otherwise

We can think of the pair H := (weights, ve) as an edge-weighted hypergraph
with incidence matrix ve.
"""
function einexpr_to_matrix!(marker::AbstractVector{Int}, ixs::AbstractVector{<:AbstractVector{L}}, iy::AbstractVector{L}, size_dict::AbstractDict{L}) where {L}
    m = length(size_dict)
    n = length(ixs) + 1
    o = sum(length, ixs) + length(iy)

    # construct incidence matrix `ve`
    #           indices
    #         [         ]
    # tensors [    ve   ]
    #         [         ]
    # we only care about the sparsity pattern
    le = sizehint!(Dict{L, Int}(), m); el = sizehint!(L[], m) # el ∘ le = id
    weights = sizehint!(Float64[], m)
    colptr = sizehint!(Int[1], n + 1)
    rowval = sizehint!(Int[], o)
    nzval = sizehint!(Int[], o)

    # for all tensors v...
    for (v, ix) in enumerate(ixs)
        # for each index l incident to v...
        for l in ix
            # let e := le[l]
            if haskey(le, l)
                e = le[l]
            else
                push!(weights, log2(size_dict[l]))
                push!(el, l)
                e = le[l] = length(el)
            end

            # if l has not been seen before in ixs[v]...
            if marker[e] < v
                # mark e as seen
                marker[e] = v

                # set ev[e, v] := 1
                push!(rowval, e)
                push!(nzval, 1)
            end           
        end

        push!(colptr, length(rowval) + 1)
    end

    # v is the outer tensor
    v = length(colptr)

    # for each index l incident to v...
    for l in iy
        # let e := le[l]
        if haskey(le, l)
            e = le[l]
        else
            push!(weights, log2(size_dict[l]))
            push!(el, l)
            e = le[l] = length(el)
        end

        # if l has not been seen before in iy...
        if marker[e] < v
            # mark e as seen
            marker[e] = v

            # set ev[e, v] = 1
            push!(rowval, e)
            push!(nzval, 1)
        end           
    end

    push!(colptr, length(rowval) + 1)

    m = length(el)
    n = length(colptr) - 1

    ev = SparseMatrixCSC{Int, Int}(m, n, colptr, rowval, nzval)
    ve = copy(transpose(ev))
    return weights, ev, ve, el
end

"""
    matrix_to_tree!(marker, weights, ev, ve, el, alg)

Construct a tree decomposition of an edge-weighted hypergraph using
the elimination algorithm `alg`. We ensure that the indices incident
to the outer tensor are contained in the root bag of the tree decomposition.
"""
function matrix_to_tree!(marker::AbstractVector{Int}, weights::AbstractVector{Float64}, ev::SparseMatrixCSC{Int, Int}, ve::SparseMatrixCSC{Int, Int}, el::AbstractVector{L}, alg::EliminationAlgorithm) where {L}
    n, m = size(ve); tag = n + 1

    # construct line graph `ee`
    #           indices
    #         [         ]
    # indices [    ee   ]
    #         [         ]
    # we only care about the sparsity pattern
    ee = ve' * ve

    # compute a tree (forest) decomposition of ee
    perm, tree = cliquetree(weights, ee; alg=ConnectedComponents(alg))

    # find the bag containing iy, call it root
    root = length(tree)

    # mark the indices in iy
    for e in view(rowvals(ev), nzrange(ev, n))
        marker[e] = tag
    end

    # the first bag containing an index in iy
    # must contain all of iy
    for (b, bag) in enumerate(tree)
        root < length(tree) && break

        for e in residual(bag)
            root < length(tree) && break

            if marker[perm[e]] == tag
                root = b
            end
        end
    end

    # make root a root node of the tree decomposition
    permute!(perm, cliquetree!(tree, root))

    # permute incidence matrix `ve` and label vector `el`.
    permute!(ve, axes(ve, 1), perm)
    permute!(el, perm)

    return tree
end

"""
    tree_to_einexpr!(marker, tree, ve, el, ixs, iy)

Transform a tree decomposition into a contraction tree.
"""
function tree_to_einexpr!(marker::AbstractVector{Int}, tree::CliqueTree{Int, Int}, ve::SparseMatrixCSC{Int, Int}, el::AbstractVector{L}, ixs::AbstractVector{<:AbstractVector{L}}, iy::AbstractVector{L}) where {L}
    n, m = size(ve); tag = n + 2

    # dynamic programming
    stack = NestedEinsum{L}[]

    # for each bag b...
    for (b, bag) in enumerate(tree)
        # sep is the separator at b
        sep = separator(bag)

        # res is the residual at b
        res = residual(bag)

        # code is the Einstein summation expression at b
        code = NestedEinsum(NestedEinsum{L}[], EinCode(Vector{L}[], L[]))

        for e in sep
            push!(code.eins.iy, el[e])
        end

        # for each index e in the residual...
        for e in res
            # for each tensor v indicent to e...
            for v in view(rowvals(ve), nzrange(ve, e))
                # if has not been seen before...
                if marker[v] < tag
                    # mark v as seen
                    marker[v] = tag

                    # if v is the outer tensor...
                    if v == n
                        # expose iy
                        append!(code.eins.iy, iy)
                    # if v is an inner tensor...
                    else
                        # make v a child of code
                        push!(code.args, NestedEinsum{L}(v))
                        push!(code.eins.ixs, ixs[v])
                    end
                end
            end
        end

        # for each child bag of b...
        for _ in childindices(tree, b)
            # the Einstein summation expression corresponding to
            # the child bag is at the top of the stack
            child = pop!(stack)

            # make this expression a child of code
            push!(code.args, child)
            push!(code.eins.ixs, child.eins.iy)
        end

        # push code to the stack
        push!(stack, code)
    end

    # we now have an expression for each root of the tree decomposition.
    # merge these together into a single Einstein expression code.
    if isone(length(stack))
        code = only(stack)
    else
        code = NestedEinsum(NestedEinsum{L}[], EinCode(Vector{L}[], L[]))
        append!(code.eins.iy, iy)

        while !isempty(stack)
            child = pop!(stack)
            push!(code.args, child)
            push!(code.eins.ixs, child.eins.iy)
        end
    end

    # append scalars to code
    for (v, ix) in enumerate(ixs)
        if isempty(ix)
            push!(code.args, NestedEinsum{L}(v))
            push!(code.eins.ixs, ix)
        end
    end

    return code
end

# no longer used
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
