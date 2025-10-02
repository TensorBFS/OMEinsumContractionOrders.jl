"""
    struct Treewidth{EL <: EliminationAlgorithm, GM} <: CodeOptimizer
    Treewidth(; alg::EL = SafeRules(BT(), MMW{3}(), MF()))

Tree width based solver. The solvers are implemented in [CliqueTrees.jl](https://algebraicjulia.github.io/CliqueTrees.jl/stable/) and [TreeWidthSolver.jl](https://github.com/ArrogantGao/TreeWidthSolver.jl). They include:

| Algorithm | Description | Time Complexity | Space Complexity |
|:-----------|:-------------|:----------------|:-----------------|
| `BFS` | breadth-first search | O(m + n) | O(n) |
| `MCS` | maximum cardinality search | O(m + n) | O(n) |
| `LexBFS` | lexicographic breadth-first search | O(m + n) | O(m + n) |
| `RCMMD` | reverse Cuthill-Mckee (minimum degree) | O(m + n) | O(m + n) |
| `RCMGL` | reverse Cuthill-Mckee (George-Liu) | O(m + n) | O(m + n) |
| `MCSM` | maximum cardinality search (minimal) | O(mn) | O(n) |
| `LexM` | lexicographic breadth-first search (minimal) | O(mn) | O(n) |
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
ba, ab -> a
├─ bcf, fac -> ba
│  ├─ e, bcef -> bcf
│  │  ├─ e
│  │  └─ bcef
│  └─ df, acd -> fac
│     ├─ df
│     └─ acd
└─ ab
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
function optimize_treewidth(optimizer::Treewidth{EL}, code::AbstractEinsum, size_dict::Dict; binary::Bool=true) where {EL}
    optimize_treewidth(optimizer, getixsv(code), getiyv(code), size_dict; binary)
end

function optimize_treewidth(optimizer::Treewidth{EL}, ixs::AbstractVector{<:AbstractVector}, iy::AbstractVector, size_dict::Dict{L,TI}; binary::Bool=true) where {L, TI, EL}
    le = Dict{L, Int}(); el = L[]

    marker = zeros(Int, max(length(size_dict), length(ixs) + 1))
    weights = Float64[]
    colptr = Int[1]
    rowval = Int[]
    nzval = Int[]    

    for (v, ix) in enumerate(ixs)
        for l in ix
            if haskey(le, l)
                e = le[l]
            else
                push!(weights, log2(size_dict[l]))
                push!(el, l)
                e = le[l] = length(el)
            end

            if marker[e] < v
                marker[e] = v
                push!(rowval, e)
                push!(nzval, 1)
            end           
        end

        push!(colptr, length(rowval) + 1)
    end

    v = length(colptr)

    for l in iy
        if haskey(le, l)
            e = le[l]
        else
            push!(weights, log2(size_dict[l]))
            push!(el, l)
            e = le[l] = length(el)
        end

        if marker[e] < v
            marker[e] = v
            push!(rowval, e)
            push!(nzval, 1)
        end           
    end

    push!(colptr, length(rowval) + 1)

    m = length(el)
    n = length(colptr) - 1

    ev = SparseMatrixCSC{Int, Int}(m, n, colptr, rowval, nzval)
    ve = copy(transpose(ev))

    # construct line graph `ee`
    #           indices
    #         [         ]
    # indices [    ee   ]
    #         [         ]
    # we only care about the sparsity pattern
    ee = ve' * ve

    # compute a tree (forest) decomposition of `ee`
    perm, tree = cliquetree(weights, ee; alg=ConnectedComponents(optimizer.alg))

    # find the bag containing `iy`, call it `root`
    root = length(tree)

    for e in view(rowvals(ev), nzrange(ev, n))
        marker[e] = -1
    end

    for (b, bag) in enumerate(tree)
        root < length(tree) && break

        for e in residual(bag)
            root < length(tree) && break

            if marker[perm[e]] == -1
                root = b
            end
        end
    end

    # make `root` a root node of the tree decomposition
    permute!(perm, cliquetree!(tree, root))

    # permute incidence matrix `ve`
    permute!(el, perm)
    permute!(ve, oneto(n), perm)

    # the vector `roots` maps each vertex to the root node
    # of its subtree
    roots = Vector{Int}(undef, m)

    for (b, bag) in enumerate(tree), e in residual(bag)
        roots[e] = b
    end

    # dynamic programming
    stack = NestedEinsum{L}[]

    for (b, bag) in enumerate(tree)
        sep = separator(bag)
        res = residual(bag)
        code = NestedEinsum(NestedEinsum{L}[], EinCode(Vector{L}[], L[]))

        for e in sep
            push!(code.eins.iy, el[e])
        end

        for e in res, v in view(rowvals(ve), nzrange(ve, e))
            if marker[v] != -2
                marker[v] = -2

                if v == n
                    append!(code.eins.iy, iy)
                else
                    push!(code.args, NestedEinsum{L}(v))
                    push!(code.eins.ixs, ixs[v])
                end
            end
        end

        for _ in childindices(tree, b)
            child = pop!(stack)
            push!(code.args, child)
            push!(code.eins.ixs, child.eins.iy)
        end

        push!(stack, code)
    end

    # we now have an expression for each root
    # of the tree decomposition
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

    # append scalars to the root
    for (v, ix) in enumerate(ixs)
        if isempty(ix)
            push!(code.args, NestedEinsum{L}(v))
            push!(code.eins.ixs, ix)
        end
    end

    if binary
        code = _optimize_code(code, size_dict, GreedyMethod())
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
