using OMEinsum.ContractionOrder: ContractionTree, log2sumexp2

export optimize_tree

struct ExprInfo
    out_dims::Vector{Int}
end
struct LeafNode
    tensorid::Int
    labels::Vector{Int}
end

mutable struct ExprTree
    left::Union{ExprTree,LeafNode}
    right::Union{ExprTree,LeafNode}
    info::ExprInfo
end
function print_expr(io::IO, expr::ExprTree, level=0)
    print(io, " "^(2*level), "(\n")
    print_expr(io, expr.left, level+1)
    print("\n")
    print_expr(io, expr.right, level+1)
    print("\n")
    print(io, " "^(2*level), ") := ", expr.info.out_dims)
end
print_expr(io::IO, expr::LeafNode, level=0) = print(io, " "^(2*level), expr.labels, " ($(expr.tensorid))")
Base.show(io::IO, expr::ExprTree) = print_expr(io, expr, 0)
Base.show(io::IO, ::MIME"text/plain", expr::ExprTree) = show(io, expr)
siblings(t::ExprTree) = Any[t.left, t.right]
siblings(::LeafNode) = Any[]
Base.copy(t::ExprTree) = ExprTree(copy(t.left), copy(t.right), ExprInfo(copy(t.info.out_dims)))
Base.copy(t::LeafNode) = LeafNode(t.tensorid, copy(t.labels))
labels(t::ExprTree) = t.info.out_dims
labels(t::LeafNode) = t.labels
maxlabel(t::ExprTree) = max(isempty(labels(t)) ? 0 : maximum(labels(t)), maxlabel(t.left), maxlabel(t.right))
maxlabel(t::LeafNode) = maximum(isempty(labels(t)) ? 0 : labels(t))
Base.:(==)(t1::ExprTree, t2::ExprTree) = _equal(t1, t2)
Base.:(==)(t1::ExprInfo, t2::ExprInfo) = _equal(t1.out_dims, t2.out_dims)
_equal(t1::ExprTree, t2::ExprTree) = _equal(t1.left, t2.left) && _equal(t1.right, t2.right) && t1.info == t2.info
_equal(t1::LeafNode, t2::LeafNode) = t1.tensorid == t2.tensorid
_equal(t1::Vector, t2::Vector) = Set(t1) == Set(t2)
_equal(a, b) = false

"""
    optimize_tree(code, size_dict; sc_target=20, βs=0.1:0.1:10, ntrials=2, niters=100, sc_weight=1.0, rw_weight=1.0, initializer=:greedy, greedy_method=OMEinsum.MinSpaceOut(), greedy_nrepeat=1)

Optimize the einsum contraction pattern specified by `code`, and edge sizes specified by `size_dict`. Key word arguments are

* `sc_target` is the target space complexity,
* `ntrails`, `βs` and `niters` are annealing parameters, doing `ntrails` indepedent annealings, each has inverse tempteratures specified by `βs`, in each temperature, do `niters` updates of the tree.
* `sc_weight` is the relative importance factor of space complexity in the loss compared with the time complexity.
* `rw_weight` is the relative importance factor of memory read and write in the loss compared with the time complexity.
* `initializer` specifies how to determine the initial configuration, it can be `:greedy` or `:random`. If it is using `:greedy` method to generate the initial configuration, it also uses two extra arguments `greedy_method` and `greedy_nrepeat`.
"""
function optimize_tree(code, size_dict; sc_target=20, βs=0.1:0.1:10, ntrials=20, niters=100, sc_weight=1.0, rw_weight=1.0, initializer=:greedy, greedy_method=OMEinsum.MinSpaceOut(), greedy_nrepeat=1)
    labels = _label_dict(OMEinsum.flatten(code))  # label to int
    inverse_map = Dict([v=>k for (k,v) in labels])
    log2_sizes = [log2.(size_dict[inverse_map[i]]) for i=1:length(labels)]
    best_tree = _initializetree(code, size_dict, initializer; greedy_method=greedy_method, greedy_nrepeat=greedy_nrepeat)
    best_tc, best_sc, best_rw = tree_timespace_complexity(best_tree, log2_sizes)
    for t = 1:ntrials
        tree = _initializetree(code, size_dict, initializer; greedy_method=greedy_method, greedy_nrepeat=greedy_nrepeat)
        optimize_tree_sa!(tree, log2_sizes; sc_target=sc_target, βs=βs, niters=niters, sc_weight=sc_weight, rw_weight=rw_weight)
        tc, sc, rw = tree_timespace_complexity(tree, log2_sizes)
        @debug "trial $t, time complexity = $tc, space complexity = $sc."
        if sc < best_sc || (sc == best_sc && exp2(tc) + rw_weight * exp2(rw) < exp2(best_tc) + rw_weight * exp2(rw))
            best_tree, best_tc, best_sc, best_rw = tree, tc, sc, rw
        end
    end
    @debug "best space complexities = $best_tc, time complexity = $best_sc, read-right complexity $best_rw."
    if best_sc > sc_target
        @warn "target space complexity not found, got: $best_sc, with time complexity $best_tc, read-right complexity $best_rw."
    end
    return NestedEinsum(best_tree, inverse_map)
end

function _initializetree(@nospecialize(code), size_dict, method; greedy_method, greedy_nrepeat)
    flatcode = OMEinsum.flatten(code)
    if method == :greedy
        labels = _label_dict(flatcode)  # label to int
        return _exprtree(optimize_greedy(flatcode, size_dict; method=greedy_method, nrepeat=greedy_nrepeat), labels)
    elseif method == :random
        return random_exprtree(flatcode)
    elseif method == :specified
        labels = _label_dict(flatcode)  # label to int
        return _exprtree(code, labels)
    else
        throw(ArgumentError("intializier `$method` is not defined!"))
    end
end

function optimize_tree_sa!(tree::ExprTree, log2_sizes; βs, niters, sc_target, sc_weight, rw_weight)
    for β in βs
        global_tc, sc, rw = tree_timespace_complexity(tree, log2_sizes)  # recompute the time complexity
        @debug "β = $β, tc = $global_tc, sc = $sc, rw = $rw"
        for _ = 1:niters
            optimize_subtree!(tree, global_tc, β, log2_sizes, sc_target, sc_weight, rw_weight)  # single sweep
        end
    end
    return tree
end

function tree_timespace_complexity(tree::LeafNode, log2_sizes)
    -Inf, sum(i->log2_sizes[i], tree.labels), -Inf
end
function tree_timespace_complexity(tree::ExprTree, log2_sizes)
    tcl, scl, rwl = tree_timespace_complexity(tree.left, log2_sizes)
    tcr, scr, rwr = tree_timespace_complexity(tree.right, log2_sizes)
    tc, sc, rw = tcscrw(labels(tree.left), labels(tree.right), labels(tree), log2_sizes)
    return fast_log2sumexp2(tc, tcl, tcr), max(sc, scl, scr), fast_log2sumexp2(rw, rwl, rwr)
end
@inline function tcscrw(ix1, ix2, iy, log2_sizes::Vector{T}) where T
    l1, l2, l3 = ix1, ix2, iy
    sc1 = isempty(l1) ? zero(T) : sum(i->(@inbounds log2_sizes[i]), l1)
    sc2 = isempty(l2) ? zero(T) : sum(i->(@inbounds log2_sizes[i]), l2)
    sc = isempty(l3) ? zero(T) : sum(i->(@inbounds log2_sizes[i]), l3)
    tc = sc
    # Note: assuming labels in `l1` being unique
    @inbounds for l in l1
        if l ∈ l2 && l ∉ l3
            tc += log2_sizes[l]
        end
    end
    rw = fast_log2sumexp2(sc, sc1, sc2)
    return tc, sc, rw
end

function random_exprtree(@nospecialize(code::EinCode{ixs, iy})) where {ixs, iy}
    labels = _label_dict(code)
    return random_exprtree([Int[labels[l] for l in ix] for ix in ixs], Int[labels[l] for l in iy], length(labels))
end

function random_exprtree(ixs::Vector{Vector{Int}}, iy::Vector{Int}, nedge::Int)
    outercount = zeros(Int, nedge)
    allcount = zeros(Int, nedge)
    for l in iy
        outercount[l] += 1
        allcount[l] += 1
    end
    for ix in ixs
        for l in ix
            allcount[l] += 1
        end
    end
    _random_exprtree(ixs, collect(1:length(ixs)), outercount, allcount)
end
function _random_exprtree(ixs::Vector{Vector{Int}}, xindices, outercount::Vector{Int}, allcount::Vector{Int})
    n = length(ixs)
    if n == 1
        return LeafNode(xindices[1], ixs[1])
    end
    mask = rand(Bool, n)
    if all(mask) || !any(mask)  # prevent invalid partition
        i = rand(1:n)
        mask[i] = ~(mask[i])
    end
    info = ExprInfo([i for i=1:length(outercount) if outercount[i]!=allcount[i] && outercount[i]!=0])
    outercount1, outercount2 = copy(outercount), copy(outercount)
    for i=1:n
        counter = mask[i] ? outercount2 : outercount1
        for l in ixs[i]
            counter[l] += 1
        end
    end
    return ExprTree(_random_exprtree(ixs[mask], xindices[mask], outercount1, allcount), _random_exprtree(ixs[(!).(mask)], xindices[(!).(mask)], outercount2, allcount), info)
end

function optimize_subtree!(tree, global_tc, β, log2_sizes, sc_target, sc_weight, rw_weight)
    rst = ruleset(tree)
    if !isempty(rst)
        rule = rand(rst)
        tc0, tc1, dsc, rw0, rw1, subout = tcsc_diff(tree, rule, log2_sizes)
        #dtc = (exp2(tc1) - exp2(tc0)) / exp2(global_tc)  # note: contribution to total tc, seems not good.
        #dtc = tc1 - tc0
        dtc = log2(exp2(tc1) + rw_weight * exp2(rw1)) - log2(exp2(tc0) + rw_weight * exp2(rw0))
        #log2(α*RW + tc) is the original `tc` term, which also optimizes read-write overheads.
        sc = _sc(tree, rule, log2_sizes)
        dE = (max(sc, sc+dsc) > sc_target ? sc_weight : 0) * dsc + dtc
        if rand() < exp(-β*dE)
            update_tree!(tree, rule, subout)
        end
        for subtree in siblings(tree)
            optimize_subtree!(subtree, global_tc, β, log2_sizes, sc_target, sc_weight, rw_weight)
        end
    end
end
_sc(tree, rule, log2_sizes) = max(__sc(tree, log2_sizes), __sc((rule == 1 || rule == 2) ? tree.left : tree.right, log2_sizes))
__sc(tree, log2_sizes) = length(labels(tree))==0 ? 0.0 : sum(l->log2_sizes[l], labels(tree))

ruleset(::LeafNode) = 1:-1
@inline function ruleset(tree::ExprTree)
    if tree.left isa ExprTree && tree.right isa ExprTree
        return 1:4
    elseif tree.left isa ExprTree
        return 1:2
    elseif tree.right isa ExprTree
        return 3:4
    else
        return 1:0
    end
end

function tcsc_diff(tree::ExprTree, rule, log2_sizes)
    if rule == 1 # (a,b), c -> (a,c),b
        return abcacb(labels(tree.left.left), labels(tree.left.right), labels(tree.left), labels(tree.right), labels(tree), log2_sizes)
    elseif rule == 2 # (a,b), c -> (c,b),a
        return abcacb(labels(tree.left.right), labels(tree.left.left), labels(tree.left), labels(tree.right), labels(tree), log2_sizes)
    elseif rule == 3 # a,(b,c) -> b,(a,c)
        return abcacb(labels(tree.right.right), labels(tree.right.left), labels(tree.right), labels(tree.left), labels(tree), log2_sizes)
    else  # a,(b,c) -> c,(b,a)
        return abcacb(labels(tree.right.left), labels(tree.right.right), labels(tree.right), labels(tree.left), labels(tree), log2_sizes)
    end
end

function abcacb(a, b, ab, c, d, log2_sizes)
    tc0, sc0, rw0, ab0 = _tcsc_merge(a, b, ab, c, d, log2_sizes)
    ac = Int[]
    for l in a
        if l ∈ b || l ∈ d  # suppose no repeated indices
            push!(ac, l)
        end
    end
    for l in c
        if l ∉ a && (l ∈ b || l ∈ d)  # suppose no repeated indices
            push!(ac, l)
        end
    end
    tc1, sc1, rw1, ab1 = _tcsc_merge(a, c, ac, b, d, log2_sizes)
    return tc0, tc1, sc1-sc0, rw0, rw1, ab1  # Note: this tc diff does not make much sense
end

function _tcsc_merge(a, b, ab, c, d, log2_sizes)
    tcl, scl, rwl = tcscrw(a, b, ab, log2_sizes)  # this is correct
    tcr, scr, rwr = tcscrw(ab, c, d, log2_sizes)
    fast_log2sumexp2(tcl, tcr), max(scl, scr), fast_log2sumexp2(rwl, rwr), ab
end

function update_tree!(tree::ExprTree, rule::Int, subout)
    if rule == 1 # (a,b), c -> (a,c),b
        b, c = tree.left.right, tree.right
        tree.left.right = c
        tree.right = b
        tree.left.info = ExprInfo(subout)
    elseif rule == 2 # (a,b), c -> (c,b),a
        a, c = tree.left.left, tree.right
        tree.left.left = c
        tree.right = a
        tree.left.info = ExprInfo(subout)
    elseif rule == 3 # a,(b,c) -> b,(a,c)
        a, b = tree.left, tree.right.left
        tree.left = b
        tree.right.left = a
        tree.right.info = ExprInfo(subout)
    else  # a,(b,c) -> c,(b,a)
        a, c = tree.left, tree.right.right
        tree.left = c
        tree.right.right = a
        tree.right.info = ExprInfo(subout)
    end
    return tree
end

# from label to integer.
function _label_dict(@nospecialize(code::EinCode{ixs, iy})) where {ixs, iy}
    ixsv, iyv = collect.(ixs), collect(iy)
    v = unique(vcat(ixsv..., iyv))
    labels = Dict(zip(v, 1:length(v)))
    return labels
end

ExprTree(code::NestedEinsum) = _exprtree(code, _label_dict(OMEinsum.flatten(code)))
function _exprtree(code::NestedEinsum, labels)
    @assert length(code.args) == 2
    ExprTree(map(enumerate(code.args)) do (i,arg)
        if arg isa Int
            LeafNode(arg, [labels[i] for i=OMEinsum.getixs(code.eins)[i]])
        else
            res = _exprtree(arg, labels)
        end
    end..., ExprInfo([labels[i] for i=OMEinsum.getiy(code.eins)]))
end

OMEinsum.NestedEinsum(expr::ExprTree) = _nestedeinsum(expr, 1:maxlabel(expr))
OMEinsum.NestedEinsum(expr::ExprTree, labelmap) = _nestedeinsum(expr, labelmap)
function _nestedeinsum(tree::ExprTree, lbs)
    eins = EinCode(((getindex.(Ref(lbs), labels(tree.left))...,), (getindex.(Ref(lbs), labels(tree.right))...,)), (getindex.(Ref(lbs), labels(tree))...,))
    NestedEinsum((_nestedeinsum(tree.left, lbs), _nestedeinsum(tree.right, lbs)), eins)
end
_nestedeinsum(tree::LeafNode, lbs) = tree.tensorid

@inline function fast_log2sumexp2(a, b)
    mm, ms = minmax(a, b)
    return log2(exp2(mm - ms) + 1) + ms
end

@inline function fast_log2sumexp2(a, b, c)
    if a > b
        if a > c
            m1, m2, ms = b, c, a
        else
            m1, m2, ms = a, b, c
        end
    else
        if b > c
            m1, m2, ms = c, a, b
        else
            m1, m2, ms = b, a, c
        end
    end
    return Base.FastMath.log2(Base.FastMath.exp2(m1 - ms) + Base.FastMath.exp2(m2 - ms) + 1) + ms
end

