struct ContractionTree
    left
    right
end

struct LegInfo
    # We use number 0, 1, 2 to denote the output tensor, the first input tensor and the second input tensor,and use e.g. `l01` to denote the set of labels that appear in both the output tensor and the input tensor.
    l1::Vector{Int}
    l2::Vector{Int}
    l12::Vector{Int}
    l01::Vector{Int}
    l02::Vector{Int}
    l012::Vector{Int}
end

"""
    tree_greedy(incidence_list, log2_sizes; α = 0.0, temperature = 0.0, nrepeat=10)

Compute greedy order, and the time and space complexities, the rows of the `incidence_list` are vertices and columns are edges.
`log2_sizes` are defined on edges.
`α` is the parameter for the loss function, for pairwise interaction, L = size(out) - α * (size(in1) + size(in2))
`temperature` is the parameter for sampling, if it is zero, the minimum loss is selected; for non-zero, the loss is selected by the Boltzmann distribution, given by p ~ exp(-loss/temperature).

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
function tree_greedy(incidence_list::IncidenceList{Int,ET}, log2_edge_sizes; α::TA = 0.0, temperature::TT = 0.0, nrepeat=10) where {ET,TA,TT}
    @assert nrepeat >= 1

    results = Vector{Tuple{ContractionTree, Vector{Float64}, Vector{Float64}}}(undef, nrepeat)

    @threads for i = 1:nrepeat
        results[i] = _tree_greedy(incidence_list, log2_edge_sizes; α = α, temperature = temperature)
    end

    best_sc = minimum([maximum(r[3]) for r in results])
    possible_ids = findall(x -> maximum(x[3]) == best_sc, results)
    possible_results = results[possible_ids]

    best_tree, best_tcs, best_scs = results[argmin([log2sumexp2(r[2]) for r in possible_results])]

    return best_tree, best_tcs, best_scs
end

function _tree_greedy(incidence_list::IncidenceList{Int,Int}, log2_edge_sizes; α::TA = 0.0, temperature::TT = 0.0) where {TA,TT}
    incidence_list = copy(incidence_list)
    n = nv(incidence_list)
    if n == 0
        return nothing
    elseif n == 1
        return collect(vertices(incidence_list))[1]
    end
    log2_tcs = Float64[] # time complexity
    log2_scs = Float64[]

    tree = Dict{Int,Any}([v=>v for v in vertices(incidence_list)])
    cost_values, cost_graph = evaluate_costs(α, incidence_list, log2_edge_sizes)
    while true
        if isempty(cost_values)
            vpool = collect(vertices(incidence_list))
            pair = minmax(vpool[1], vpool[2])  # to prevent empty intersect
        else
            pair = find_best_cost!(temperature, cost_values, cost_graph)
        end
        log2_tc_step, sc, code = contract_pair!(incidence_list, pair..., log2_edge_sizes)
        push!(log2_tcs, log2_tc_step)
        push!(log2_scs, space_complexity(incidence_list, log2_edge_sizes))
        if nv(incidence_list) > 1
            tree[pair[1]] = ContractionTree(tree[pair[1]], tree[pair[2]])
        else
            return ContractionTree(tree[pair[1]], tree[pair[2]]), log2_tcs, log2_scs
        end
        update_costs!(cost_values, cost_graph, pair..., α, incidence_list, log2_edge_sizes)
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

function evaluate_costs(α::TA, incidence_list::IncidenceList{Int,Int}, log2_edge_sizes) where {TA}
    # initialize cost values
    cost_values = PriorityQueue{Tuple{Int,Int},Float64}()#Dict{Tuple{VT,VT},Float64}()
    cost_graph = SimpleGraph(nv(incidence_list))
    for vi in vertices(incidence_list)
        for vj in neighbors(incidence_list, vi)
            if vj > vi
                enqueue!(cost_values, (vi,vj), greedy_loss(α, incidence_list, log2_edge_sizes, vi, vj))
                add_edge!(cost_graph, vi, vj)
            end
        end
    end
    return cost_values, cost_graph
end

function update_costs!(cost_values, cost_graph, va, vb, α::TA, incidence_list::IncidenceList{Int,Int}, log2_edge_sizes) where {TA}
    for vj in neighbors(incidence_list, va)
        vx, vy = minmax(vj, va)
        if has_edge(cost_graph, vx, vy)
            delete!(cost_values, (vx,vy))
            enqueue!(cost_values, (vx,vy), greedy_loss(α, incidence_list, log2_edge_sizes, vx, vy))
        else
            @info "haha!"
            enqueue!(cost_values, (vx,vy), greedy_loss(α, incidence_list, log2_edge_sizes, vx, vy))
            add_edge!(cost_graph, vx, vy)
        end
    end
    for vj in copy(Graphs.neighbors(cost_graph, vb))
        vx, vy = minmax(vj, vb)
        delete!(cost_values, (vx,vy))
        rem_edge!(cost_graph, vx, vy)
    end
end

function find_best_cost!(temperature::TT, cost_values::PriorityQueue{PT}, cost_graph) where {PT,TT}
    length(cost_values) < 1 && error("cost value information missing")
    !iszero(temperature) && @warn "non-zero temperature is not supported any more, using temperature = 0.0"
    vx, vy = dequeue!(cost_values)
    rem_edge!(cost_graph, vx, vy)
    return vx, vy
    # return sample_best_cost(cost_values, temperature)
end

# function sample_best_cost(cost_values::Dict{PT}, t::T) where {PT, T}
#     length(cost_values) < 1 && error("cost value information missing")
#     vals = [v for v in values(cost_values)]
#     prob = exp.( - vals ./ t)
#     vc = [k for (k, v) in cost_values]
#     sample(vc, Weights(prob))
# end

function analyze_contraction(incidence_list::IncidenceList{Int,Int}, vi::Int, vj::Int)
    ei = edges(incidence_list, vi)
    ej = edges(incidence_list, vj)
    leg012,leg12,leg1,leg2,leg01,leg02 = Int[], Int[], Int[], Int[], Int[], Int[]
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

function greedy_loss(α, incidence_list, log2_edge_sizes, vi, vj)
    log2dim(legs) = isempty(legs) ? 0 : sum(l->log2_edge_sizes[l], legs)  # for 1.5, you need this patch because `init` kw is not allowed.
    legs = analyze_contraction(incidence_list, vi, vj)
    D1,D2,D12,D01,D02,D012 = log2dim.(getfield.(Ref(legs), 1:6))
    loss = exp2(D01+D02+D012) - α * (exp2(D01+D12+D012) + exp2(D02+D12+D012))  # out - in
    return loss
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
    optimize_greedy(eincode, size_dict; α = 0.0, temperature = 0.0, nrepeat=10)

Greedy optimizing the contraction order and return a `NestedEinsum` object.
Check the docstring of `tree_greedy` for detailed explaination of other input arguments.
"""
function optimize_greedy(code::EinCode{L}, size_dict::Dict; α::TA = 0.0, temperature::TT = 0.0, nrepeat=10) where {L,TA,TT}
    symbols = unique!(reduce(vcat, [getixsv(code)..., getiyv(code)]))
    symbol2int = Dict(symbols .=> 1:length(symbols))
    ixs = [[symbol2int[i] for i in ix] for ix in getixsv(code)]
    iy = [symbol2int[i] for i in getiyv(code)]
    size_dict = Dict([k=>size_dict[i] for (i,k) in symbol2int])
    result = optimize_greedy(ixs, iy, size_dict; α = α, temperature = temperature, nrepeat=nrepeat)
    inverse_map = Dict([v=>k for (k,v) in symbol2int])
    result = convert_label(result, inverse_map)
end
function convert_label(ne::NestedEinsum, labelmap::Dict{T1,T2}) where {T1,T2}
    isleaf(ne) && return NestedEinsum{T2}(ne.tensorindex)
    eins = EinCode([getindex.(Ref(labelmap), ix) for ix in ne.eins.ixs], getindex.(Ref(labelmap), ne.eins.iy))
    NestedEinsum([convert_label(arg, labelmap) for arg in ne.args], eins)
end

function optimize_greedy(ixs::AbstractVector{<:AbstractVector}, iy::AbstractVector, size_dict::Dict{L,TI}; α::TA = 0.0, temperature::TT = 0.0, nrepeat=10) where {L, TI, TA, TT}
    if length(ixs) <= 2
        return NestedEinsum(NestedEinsum{L}.(1:length(ixs)), EinCode(ixs, iy))
    end
    log2_edge_sizes = Dict{L,Float64}()
    for (k, v) in size_dict
        log2_edge_sizes[k] = log2(v)
    end
    incidence_list = IncidenceList(Dict([i=>ixs[i] for i=1:length(ixs)]); openedges=iy)
    tree, _, _ = tree_greedy(incidence_list, log2_edge_sizes; α = α, temperature = temperature, nrepeat=nrepeat)
    parse_eincode!(incidence_list, tree, 1:length(ixs))[2]
end
function optimize_greedy(code::NestedEinsum, size_dict; α::TA = 0.0, temperature::TT = 0.0, nrepeat=10) where {TT, TA}
    isleaf(code) && return code
    args = optimize_greedy.(code.args, Ref(size_dict); α = α, temperature = temperature, nrepeat=nrepeat)
    if length(code.args) > 2
        # generate coarse grained hypergraph.
        nested = optimize_greedy(code.eins, size_dict; α = α, temperature = temperature, nrepeat=nrepeat)
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
    GreedyMethod(; α = 0.0, temperature = 0.0, nrepeat=10)

The fast but poor greedy optimizer. Input arguments are

    * `α` is the parameter for the loss function, for pairwise interaction, L = size(out) - α * (size(in1) + size(in2))
    * `temperature` is the parameter for sampling, if it is zero, the minimum loss is selected; for non-zero, the loss is selected by the Boltzmann distribution, given by p ~ exp(-loss/temperature).
    * `nrepeat` is the number of repeatition, returns the best contraction order.
"""
Base.@kwdef struct GreedyMethod{TA, TT} <: CodeOptimizer
    α::TA = 0.0
    temperature::TT = 0.0
    nrepeat::Int = 10
end