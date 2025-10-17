struct ContractionTree
    left
    right
end

struct LegInfo{ET}
    # We use number 0, 1, 2 to denote the output tensor, the first input tensor and the second input tensor,and use e.g. `l01` to denote the set of labels that appear in both the output tensor and the input tensor.
    l1::Vector{ET}
    l2::Vector{ET}
    l12::Vector{ET}
    l01::Vector{ET}
    l02::Vector{ET}
    l012::Vector{ET}
end

"""
    tree_greedy(incidence_list, log2_sizes; α = 0.0, temperature = 0.0)

Compute greedy order, and the time and space complexities, the rows of the `incidence_list` are vertices and columns are edges.
`log2_sizes` are defined on edges.
`α` is the parameter for the loss function, for pairwise interaction, L = size(out) - α * (size(in1) + size(in2))
`temperature` is the parameter for sampling, if it is zero, the minimum loss is selected; for non-zero, the loss is selected by the Boltzmann distribution, given by p ~ exp(-loss/temperature).
"""
function tree_greedy(incidence_list::IncidenceList{Int, ET}, log2_edge_sizes; α::TA = 0.0, temperature::TT = 0.0) where {TA,TT, ET}
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
        push!(log2_scs, sc)   # NOTE: the original one is computed with `space_complexity`, why did I made such a mistake?
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

function evaluate_costs(α::TA, incidence_list::IncidenceList{Int,ET}, log2_edge_sizes) where {TA,ET}
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

function update_costs!(cost_values, cost_graph, va, vb, α::TA, incidence_list::IncidenceList{Int,ET}, log2_edge_sizes) where {TA,ET}
    for vj in neighbors(incidence_list, va)
        vx, vy = minmax(vj, va)
        if has_edge(cost_graph, vx, vy)
            delete!(cost_values, (vx,vy))
            enqueue!(cost_values, (vx,vy), greedy_loss(α, incidence_list, log2_edge_sizes, vx, vy))
        else
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
    pair, cost = dequeue_pair!(cost_values)
    if iszero(temperature) || isempty(cost_values)
        rem_edge!(cost_graph, pair...)
        return pair
    else
        pair2, cost2 = dequeue_pair!(cost_values)
        if rand() < exp(-(cost2 - cost) / temperature)   # pick 2
            enqueue!(cost_values, pair, cost)
            rem_edge!(cost_graph, pair2...)
            return pair2
        else
            enqueue!(cost_values, pair2, cost2)
            rem_edge!(cost_graph, pair...)
            return pair
        end
    end
end

function analyze_contraction(incidence_list::IncidenceList{Int,ET}, vi::Int, vj::Int) where {ET}
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
function parse_eincode!(::IncidenceList{IT,LT}, tree, vertices_order, size_dict, level=0) where {IT,LT}
    ti = findfirst(==(tree), vertices_order)::Int
    return ti, NestedEinsum{LT}(ti)
end

function parse_eincode!(incidence_list::IncidenceList{IT,LT}, tree::ContractionTree, vertices_order, size_dict, level=0) where {IT,LT}
    ti, codei = parse_eincode!(incidence_list, tree.left, vertices_order, size_dict, level+1)
    tj, codej = parse_eincode!(incidence_list, tree.right, vertices_order, size_dict, level+1)
    _, _, code = contract_pair!(incidence_list, vertices_order[ti], vertices_order[tj], size_dict)
    return ti, NestedEinsum([codei, codej], EinCode([code.first...], level==0 ? incidence_list.openedges : code.second))
end

function parse_eincode(incidence_list::IncidenceList, tree::ContractionTree; vertices = collect(keys(incidence_list.v2e)))
    size_dict = Dict([v=>0 for (i,v) in enumerate(vertices)])  # dummy size_dict
    parse_eincode!(copy(incidence_list), tree, vertices, size_dict)[2]
end

function parse_nested(code::EinCode{LT}, tree::ContractionTree) where LT
    ixs, iy = getixsv(code), getiyv(code)
    incidence_list = IncidenceList(Dict([i=>ixs[i] for i=1:length(ixs)]); openedges=iy)
    parse_eincode!(incidence_list, tree, 1:length(ixs), size_dict)[2]
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
    optimize_greedy(eincode, size_dict; α, temperature)

Greedy optimizing the contraction order and return a `NestedEinsum` object.
Check the docstring of `tree_greedy` for detailed explaination of other input arguments.
"""
function optimize_greedy(code::EinCode{L}, size_dict::Dict{L, T2}; α, temperature) where {L, T2}
    optimize_greedy(getixsv(code), getiyv(code), size_dict; α, temperature)
end
function convert_label(ne::NestedEinsum, labelmap::Dict{T1,T2}) where {T1,T2}
    isleaf(ne) && return NestedEinsum{T2}(ne.tensorindex)
    eins = EinCode([getindex.(Ref(labelmap), ix) for ix in ne.eins.ixs], getindex.(Ref(labelmap), ne.eins.iy))
    NestedEinsum([convert_label(arg, labelmap) for arg in ne.args], eins)
end

function optimize_greedy(ixs::AbstractVector{<:AbstractVector}, iy::AbstractVector, size_dict::Dict{L}; α, temperature) where {L}
    if length(ixs) <= 2
        return NestedEinsum(NestedEinsum{L}.(1:length(ixs)), EinCode(ixs, iy))
    end
    log2_edge_sizes = Dict{L,Float64}()
    for (k, v) in size_dict
        log2_edge_sizes[k] = log2(v)
    end
    incidence_list = IncidenceList(Dict([i=>ixs[i] for i=1:length(ixs)]); openedges=iy)
    tree, _, _ = tree_greedy(incidence_list, log2_edge_sizes; α, temperature)
    parse_eincode!(incidence_list, tree, 1:length(ixs), size_dict)[2]
end

function optimize_greedy(code::E, size_dict; α, temperature) where E <: NestedEinsum
    # construct first-child next-sibling representation of `code`
    queue = [code]
    child = [0]
    brother = [0]

    for (i, code) in enumerate(queue)
        if !isleaf(code)
            for arg in code.args
                push!(queue, arg)
                push!(brother, child[i])
                push!(child, 0)
                child[i] = length(queue)
            end
        end
    end

    # construct postordering of `code`
    order = similar(queue)
    stack = [1]

    for i in eachindex(queue)
        j = pop!(stack); k = child[j]

        while !iszero(k)
            push!(stack, j); child[j] = brother[k]; j = k; k = child[j]
        end

        order[i] = queue[j]
    end

    # optimize `code` using dynamic programming
    empty!(queue)

    for code in order
        if isleaf(code)
            push!(queue, code)
        else
            args = E[]

            for _ in code.args
                push!(args, pop!(queue))
            end

            if length(args) > 2
                code = replace_args(optimize_greedy(code.eins, size_dict; α, temperature), args)
            else
                code = NestedEinsum(args, code.eins)
            end

            push!(queue, code)
        end
    end

    return only(queue)
end

function replace_args(nested::NestedEinsum{LT}, trueargs) where LT
    isleaf(nested) && return trueargs[nested.tensorindex]
    NestedEinsum(replace_args.(nested.args, Ref(trueargs)), nested.eins)
end

"""
    GreedyMethod{MT}
    GreedyMethod(; α = 0.0, temperature = 0.0)

It may not be optimal, but it is fast.

# Fields
- `α` is the parameter for the loss function, for pairwise interaction, L = size(out) - α * (size(in1) + size(in2))
- `temperature` is the parameter for sampling, if it is zero, the minimum loss is selected; for non-zero, the loss is selected by the Boltzmann distribution, given by p ~ exp(-loss/temperature).
"""
Base.@kwdef struct GreedyMethod{TA, TT} <: CodeOptimizer
    α::TA = 0.0
    temperature::TT = 0.0
end
