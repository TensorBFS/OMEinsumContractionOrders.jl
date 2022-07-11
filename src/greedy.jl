struct ContractionTree
    left
    right
end

struct MinSpaceOut end
struct MinSpaceDiff end

struct LegInfo{ET}
    l1::Vector{ET}
    l2::Vector{ET}
    l12::Vector{ET}
    l01::Vector{ET}
    l02::Vector{ET}
    l012::Vector{ET}
end

"""
    tree_greedy(incidence_list, log2_sizes; method=MinSpaceOut())

Compute greedy order, and the time and space complexities, the rows of the `incidence_list` are vertices and columns are edges.
`log2_sizes` are defined on edges.

```julia
julia> code = ein"(abc,cde),(ce,sf,j),ak->ael"
aec, ec, ak -> ael
├─ ce, sf, j -> ec
│  ├─ sf
│  ├─ j
│  └─ ce
├─ ak
└─ abc, cde -> aec
   ├─ cde
   └─ abc


julia> optimize_greedy(code, Dict([c=>2 for c in "abcdefjkls"]))
ae, ak -> ea
├─ ak
└─ aec, ec -> ae
   ├─ ce,  -> ce
   │  ├─ sf, j -> 
   │  │  ├─ j
   │  │  └─ sf
   │  └─ ce
   └─ abc, cde -> aec
      ├─ cde
      └─ abc
```
"""
function tree_greedy(incidence_list::IncidenceList{VT,ET}, log2_edge_sizes; method=MinSpaceOut(), nrepeat=10) where {VT,ET}
    @assert nrepeat >= 1
    best_tree, best_tcs, best_scs = _tree_greedy(incidence_list, log2_edge_sizes; method=method)
    best_tc, best_sc = log2sumexp2(best_tcs), maximum(best_scs)
    for _ = 1:nrepeat-1
        tree, tcs, scs = _tree_greedy(incidence_list, log2_edge_sizes; method=method)
        tc, sc = log2sumexp2(tcs), maximum(scs)
        if sc < best_sc || (sc <= best_sc && tc < best_tc)
            best_tcs, best_scs, best_tc, best_sc, best_tree = tcs, scs, tc, sc, tree
        end
    end
    return best_tree, best_tcs, best_scs
end

function _tree_greedy(incidence_list::IncidenceList{VT,ET}, log2_edge_sizes; method=MinSpaceOut()) where {VT,ET}
    incidence_list = copy(incidence_list)
    n = nv(incidence_list)
    if n == 0
        return nothing
    elseif n == 1
        return collect(vertices(incidence_list))[1]
    end
    log2_tcs = Float64[] # time complexity
    log2_scs = Float64[]

    tree = Dict{VT,Any}([v=>v for v in vertices(incidence_list)])
    cost_values = evaluate_costs(method, incidence_list, log2_edge_sizes)
    while true
        if length(cost_values) == 0
            vpool = collect(vertices(incidence_list))
            pair = minmax(vpool[1], vpool[2])  # to prevent empty intersect
        else
            pair = find_best_cost(cost_values)
        end
        log2_tc_step, sc, code = contract_pair!(incidence_list, pair..., log2_edge_sizes)
        push!(log2_tcs, log2_tc_step)
        push!(log2_scs, space_complexity(incidence_list, log2_edge_sizes))
        if nv(incidence_list) > 1
            tree[pair[1]] = ContractionTree(tree[pair[1]], tree[pair[2]])
        else
            return ContractionTree(tree[pair[1]], tree[pair[2]]), log2_tcs, log2_scs
        end
        update_costs!(cost_values, pair..., method, incidence_list, log2_edge_sizes)
    end
end

function contract_pair!(incidence_list, vi, vj, log2_edge_sizes)
    log2dim(legs) = isempty(legs) ? 0 : sum(l->log2_edge_sizes[l], legs)  # for 1.5, you need this patch because `init` kw is not allowed.
    # compute time complexity and output tensor
    legsets = analyze_contraction(incidence_list, vi, vj)
    D12,D01,D02,D012 = log2dim.(getfield.(Ref(legsets),3:6))
    tc = D12+D01+D02+D012  # dangling legs D1 and D2 do not contribute

    # einsum code
    eout = legsets.l01 ∪ legsets.l02 ∪ legsets.l012
    code = (edges(incidence_list, vi), edges(incidence_list, vj)) => eout
    sc = log2dim(eout)

    # change incidence_list
    delete_vertex!(incidence_list, vj)
    change_edges!(incidence_list, vi, eout)
    for e in eout
        replace_vertex!(incidence_list, e, vj=>vi)
    end
    remove_edges!(incidence_list, legsets.l1 ∪ legsets.l2 ∪ legsets.l12)
    return tc, sc, code
end

function evaluate_costs(method, incidence_list::IncidenceList{VT,ET}, log2_edge_sizes) where {VT,ET}
    # initialize cost values
    cost_values = Dict{Tuple{VT,VT},Float64}()
    for vi = vertices(incidence_list)
        for vj in neighbors(incidence_list, vi)
            if vj > vi
                cost_values[(vi,vj)] = greedy_loss(method, incidence_list, log2_edge_sizes, vi, vj)
            end
        end
    end
    return cost_values
end

function update_costs!(cost_values, va, vb, method, incidence_list::IncidenceList{VT,ET}, log2_edge_sizes) where {VT,ET}
    for vj in neighbors(incidence_list, va)
        vx, vy = minmax(vj, va)
        cost_values[(vx,vy)] = greedy_loss(method, incidence_list, log2_edge_sizes, vx, vy)
    end
    for k in keys(cost_values)
        if vb ∈ k
            delete!(cost_values, k)
        end
    end
end

function find_best_cost(cost_values::Dict{PT}) where PT
    length(cost_values) < 1 && error("cost value information missing")
    minval = minimum(Base.values(cost_values))
    pairs = PT[]
    for (k, v) in cost_values
        if v == minval
            push!(pairs, k)
        end
    end
    return rand(pairs)
end

function analyze_contraction(incidence_list::IncidenceList{VT,ET}, vi::VT, vj::VT) where {VT,ET}
    ei = edges(incidence_list, vi)
    ej = edges(incidence_list, vj)
    leg012,leg12,leg1,leg2,leg01,leg02 = ET[], ET[], ET[], ET[], ET[], ET[]
    # external legs
    for leg in ei ∪ ej
        isext = leg ∈ incidence_list.openedges || !all(x->x==vi || x==vj, vertices(incidence_list, leg))
        if isext
            if leg ∈ ei
                if leg ∈ ej
                    push!(leg012, leg)
                else
                    push!(leg01, leg)
                end
            else
                push!(leg02, leg)
            end
        else
            if leg ∈ ei
                if leg ∈ ej
                    push!(leg12, leg)
                else
                    push!(leg1, leg)
                end
            else
                push!(leg2, leg)
            end
        end
    end
    return LegInfo(leg1, leg2, leg12, leg01, leg02, leg012)
end

function greedy_loss(::MinSpaceOut, incidence_list, log2_edge_sizes, vi, vj)
    log2dim(legs) = isempty(legs) ? 0 : sum(l->log2_edge_sizes[l], legs)  # for 1.5, you need this patch because `init` kw is not allowed.
    legs = analyze_contraction(incidence_list, vi, vj)
    log2dim(legs.l01)+log2dim(legs.l02)+log2dim(legs.l012)
end

function greedy_loss(::MinSpaceDiff, incidence_list, log2_edge_sizes, vi, vj)
    log2dim(legs) = isempty(legs) ? 0 : sum(l->log2_edge_sizes[l], legs)  # for 1.5, you need this patch because `init` kw is not allowed.
    legs = analyze_contraction(incidence_list, vi, vj)
    D1,D2,D12,D01,D02,D012 = log2dim.(getfield.(Ref(legs), 1:6))
    exp2(D01+D02+D012) - exp2(D01+D12+D012) - exp2(D02+D12+D012)  # out - in
end

function space_complexity(incidence_list, log2_sizes)
    sc = 0.0
    for v in vertices(incidence_list)
        for e in edges(incidence_list, v)
            sc += log2_sizes[e]
        end
    end
    return sc
end

function contract_tree!(incidence_list::IncidenceList, tree::ContractionTree, log2_edge_sizes, tcs, scs)
    vi = tree.left isa ContractionTree ? contract_tree!(incidence_list, tree.left, log2_edge_sizes, tcs, scs) : tree.left
    vj = tree.right isa ContractionTree ? contract_tree!(incidence_list, tree.right, log2_edge_sizes, tcs, scs) : tree.right
    tc, sc, code = contract_pair!(incidence_list, vi, vj, log2_edge_sizes)
    push!(tcs, tc)
    push!(scs, sc)
    return vi
end

#################### parse to code ####################
function parse_eincode!(::IncidenceList{IT,LT}, tree, vertices_order, level=0) where {IT,LT}
    ti = findfirst(==(tree), vertices_order)
    ti, NestedEinsum{LT}(ti)
end

function parse_eincode!(incidence_list::IncidenceList{IT,LT}, tree::ContractionTree, vertices_order, level=0) where {IT,LT}
    ti, codei = parse_eincode!(incidence_list, tree.left, vertices_order, level+1)
    tj, codej = parse_eincode!(incidence_list, tree.right, vertices_order, level+1)
    dummy = Dict([e=>0 for e in keys(incidence_list.e2v)])
    _, _, code = contract_pair!(incidence_list, vertices_order[ti], vertices_order[tj], dummy)
    ti, NestedEinsum([codei, codej], EinCode([code.first...], level==0 ? incidence_list.openedges : code.second))
end

function parse_eincode(incidence_list::IncidenceList, tree::ContractionTree; vertices = collect(keys(incidence_list.v2e)))
    parse_eincode!(copy(incidence_list), tree, vertices)[2]
end

function parse_nested(code::EinCode{LT}, tree::ContractionTree) where LT
    ixs, iy = getixsv(code), getiyv(code)
    incidence_list = IncidenceList(Dict([i=>ixs[i] for i=1:length(ixs)]); openedges=iy)
    parse_eincode!(incidence_list, tree, 1:length(ixs))[2]
end

function parse_tree(ein, vertices)
    if isleaf(ein)
        vertices[ein.tensorindex]
    else
        if length(ein.args) != 2
            error("This eincode is not a binary tree.")
        end
        left, right = parse_tree.(ein.args, Ref(vertices))
        ContractionTree(left, right)
    end
end

"""
    optimize_greedy(eincode, size_dict; method=MinSpaceOut(), nrepeat=10)

Greedy optimizing the contraction order and return a `NestedEinsum` object. Methods are
* `MinSpaceOut`, always choose the next contraction that produces the minimum output tensor.
* `MinSpaceDiff`, always choose the next contraction that minimizes the total space.
"""
function optimize_greedy(code::EinCode{L}, size_dict::Dict; method=MinSpaceOut(), nrepeat=10) where {L}
    optimize_greedy(getixsv(code), getiyv(code), size_dict; method=method, nrepeat=nrepeat)
end
function optimize_greedy(ixs::AbstractVector{<:AbstractVector}, iy::AbstractVector, size_dict::Dict{L,TI}; method=MinSpaceOut(), nrepeat=10) where {L, TI}
    if length(ixs) <= 2
        return NestedEinsum(NestedEinsum{L}.(1:length(ixs)), EinCode(ixs, iy))
    end
    log2_edge_sizes = Dict{L,Float64}()
    for (k, v) in size_dict
        log2_edge_sizes[k] = log2(v)
    end
    incidence_list = IncidenceList(Dict([i=>ixs[i] for i=1:length(ixs)]); openedges=iy)
    tree, _, _ = tree_greedy(incidence_list, log2_edge_sizes; method=method, nrepeat=nrepeat)
    parse_eincode!(incidence_list, tree, 1:length(ixs))[2]
end
function optimize_greedy(code::NestedEinsum, size_dict; method=MinSpaceOut(), nrepeat=10)
    isleaf(code) && return code
    args = optimize_greedy.(code.args, Ref(size_dict); method=method, nrepeat=nrepeat)
    if length(code.args) > 2
        # generate coarse grained hypergraph.
        nested = optimize_greedy(code.eins, size_dict; method=method, nrepeat=nrepeat)
        replace_args(nested, args)
    else
        NestedEinsum(args, code.eins)
    end
end

function replace_args(nested::NestedEinsum{LT}, trueargs) where LT
    isleaf(nested) && return trueargs[nested.tensorindex]
    NestedEinsum(replace_args.(nested.args, Ref(trueargs)), nested.eins)
end

"""
    GreedyMethod{MT}
    GreedyMethod(; method=MinSpaceOut(), nrepeat=10)

The fast but poor greedy optimizer. Input arguments are

* `method` is `MinSpaceDiff()` or `MinSpaceOut`.
    * `MinSpaceOut` choose one of the contraction that produces a minimum output tensor size,
    * `MinSpaceDiff` choose one of the contraction that decrease the space most.
* `nrepeat` is the number of repeatition, returns the best contraction order.
"""
Base.@kwdef struct GreedyMethod{MT} <: CodeOptimizer
    method::MT = MinSpaceOut()
    nrepeat::Int = 10
end
