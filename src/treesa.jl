################### expression tree ###################
# `ExprInfo` stores the node information.
# * `out_dims` is the output dimensions of this tree/subtree.
# * `tensorid` specifies the tensor index for leaf nodes. It is `-1` is for non-leaf node.
struct ExprInfo
    out_dims::Vector{Int}
    tensorid::Int
end
ExprInfo(out_dims::Vector{Int}) = ExprInfo(out_dims, -1)

# `ExprTree` is the expression tree for tensor contraction (or contraction tree), it is a binary tree (including leaf nodes without siblings).
# `left` and `right` are left and right branches, they are either both specified (non-leaf) or both unspecified (leaf), see [`isleaf`](@ref) function.
# `ExprTree()` for constructing a leaf node,
# `ExprTree(left, right, info)` for constructing a non-leaf node.
mutable struct ExprTree
    left::ExprTree
    right::ExprTree
    info::ExprInfo
    ExprTree(info) = (res = new(); res.info=info; res)
    ExprTree(left, right, info) = new(left, right, info)
end
function print_expr(io::IO, expr::ExprTree, level=0)
    isleaf(expr) && return print(io, " "^(2*level), labels(expr), " ($(expr.info.tensorid))")
    print(io, " "^(2*level), "(\n")
    print_expr(io, expr.left, level+1)
    print("\n")
    print_expr(io, expr.right, level+1)
    print("\n")
    print(io, " "^(2*level), ") := ", labels(expr))
end
# if `expr` is a leaf, it should have `left` and `right` fields both unspecified.
isleaf(expr::ExprTree) = !isdefined(expr, :left)
Base.show(io::IO, expr::ExprTree) = print_expr(io, expr, 0)
Base.show(io::IO, ::MIME"text/plain", expr::ExprTree) = show(io, expr)
siblings(t::ExprTree) = isleaf(t) ? ExprTree[] : ExprTree[t.left, t.right]
Base.copy(t::ExprTree) = isleaf(t) ? ExprTree(t.info) : ExprTree(copy(t.left), copy(t.right), copy(t.info))
Base.copy(info::ExprInfo) = ExprInfo(copy(info.out_dims), info.tensorid)
# output tensor labels
labels(t::ExprTree) = t.info.out_dims
# find the maximum label recursively, this is a helper function for converting an expression tree back to einsum.
maxlabel(t::ExprTree) = isleaf(t) ? maximum(isempty(labels(t)) ? 0 : labels(t)) : max(isempty(labels(t)) ? 0 : maximum(labels(t)), maxlabel(t.left), maxlabel(t.right))
# comparison between `ExprTree`s, mainly for testing
Base.:(==)(t1::ExprTree, t2::ExprTree) = _equal(t1, t2)
Base.:(==)(t1::ExprInfo, t2::ExprInfo) = _equal(t1.out_dims, t2.out_dims) && t1.tensorid == t2.tensorid
function _equal(t1::ExprTree, t2::ExprTree)
    isleaf(t1) != isleaf(t2) && return false
    isleaf(t1) ? t1.info == t2.info : _equal(t1.left, t2.left) && _equal(t1.right, t2.right) && t1.info == t2.info
end
_equal(t1::Vector, t2::Vector) = Set(t1) == Set(t2)


############# random expression tree ###############
function random_exprtree(code::EinCode)
    ixs, iy = getixsv(code), getiyv(code)
    labels = _label_dict(ixs, iy)
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
        return ExprTree(ExprInfo(ixs[1], xindices[1]))
    end
    mask = rand(Bool, n)
    if all(mask) || !any(mask)  # prevent invalid partition
        i = rand(1:n)
        mask[i] = ~(mask[i])
    end
    info = ExprInfo(Int[i for i=1:length(outercount) if outercount[i]!=allcount[i] && outercount[i]!=0])
    outercount1, outercount2 = copy(outercount), copy(outercount)
    for i=1:n
        counter = mask[i] ? outercount2 : outercount1
        for l in ixs[i]
            counter[l] += 1
        end
    end
    return ExprTree(_random_exprtree(ixs[mask], xindices[mask], outercount1, allcount), _random_exprtree(ixs[(!).(mask)], xindices[(!).(mask)], outercount2, allcount), info)
end


##################### convert a contraction tree back to a nested einsum ####################
NestedEinsum(expr::ExprTree) = _nestedeinsum(expr, 1:maxlabel(expr))
NestedEinsum(expr::ExprTree, labelmap) = _nestedeinsum(expr, labelmap)
function _nestedeinsum(tree::ExprTree, lbs::Union{AbstractVector{LT}, Dict{Int,LT}}) where LT
    isleaf(tree) && return NestedEinsum{LT}(tree.info.tensorid)
    eins = EinCode([getindex.(Ref(lbs), labels(tree.left)), getindex.(Ref(lbs), labels(tree.right))], getindex.(Ref(lbs), labels(tree)))
    NestedEinsum([_nestedeinsum(tree.left, lbs), _nestedeinsum(tree.right, lbs)], eins)
end


##################### The main program ##############################
"""
    TreeSA{IT} <: CodeOptimizer
    TreeSA(; βs=collect(0.01:0.05:15), ntrials=10, niters=50, initializer=:greedy, score=ScoreFunction())

Optimize the einsum contraction pattern using the simulated annealing on tensor expression tree.

# Fields
- `ntrials`, `βs` and `niters` are annealing parameters, doing `ntrials` indepedent annealings, each has inverse tempteratures specified by `βs`, in each temperature, do `niters` updates of the tree.
- `initializer` specifies how to determine the initial configuration, it can be `:greedy`, `:random` or `:specified`. If the initializer is `:specified`, the input `code` should be a `NestedEinsum` object.
- `score` specifies the score function to evaluate the quality of the contraction tree, it is a function of time complexity, space complexity and read-write complexity.

# References
- [Recursive Multi-Tensor Contraction for XEB Verification of Quantum Circuits](https://arxiv.org/abs/2108.05665)

# Breaking changes:
- `nslices` is removed, since the slicing part is now separated from the optimization part, see `slice_code` function and `TreeSASlicer`.
- `greedy_method` is removed. If you want to have detailed control of the initializer, please pre-optimize the code with another method and then use `:specified` to initialize the tree.
"""
Base.@kwdef struct TreeSA{IT} <: CodeOptimizer
    βs::IT = 0.01:0.05:15
    ntrials::Int = 10
    niters::Int = 50
    initializer::Symbol = :greedy
    score::ScoreFunction = ScoreFunction()
end

# this is the main function
"""
    optimize_tree(code, size_dict; βs, ntrials, niters, initializer, score)

Optimize the einsum contraction pattern specified by `code`, and edge sizes specified by `size_dict`.
Check the docstring of [`TreeSA`](@ref) for detailed explaination of other input arguments.
"""
function optimize_tree(code::AbstractEinsum, size_dict::Dict{LT,Int}; βs, ntrials, niters, initializer, score) where LT
    # get input labels (`getixsv`) and output labels (`getiyv`) in the einsum code.
    ixs, iy = getixsv(code), getiyv(code)
    ninputs = length(ixs)  # number of input tensors
    if ninputs <= 2  # number of input tensors ≤ 2, can not be optimized
        return NestedEinsum(NestedEinsum{LT}.(1:ninputs), EinCode(ixs, iy))
    end

    ###### Stage 1: preprocessing ######
    labels = _label_dict(ixs, iy)  # map labels to integers
    inverse_map = Dict([v=>k for (k,v) in labels])  # the inverse transformation, map integers to labels
    log2_sizes = [log2.(size_dict[inverse_map[i]]) for i=1:length(labels)]   # use `log2` sizes in computing time

    if ntrials <= 0  # no optimization at all, then 1). initialize an expression tree and 2). convert back to nested einsum.
        best_tree = _initializetree(code, size_dict, initializer)
        return NestedEinsum(best_tree, inverse_map)
    end

    ###### Stage 2: computing ######
    # create vectors to store optimized 1). expression tree, 2). time complexities, 3). space complexities, 4). read-write complexities.
    local best_tree, best_tc, best_sc, best_rw
    @threads for t = 1:ntrials  # multi-threading on different trials, use `JULIA_NUM_THREADS=5 julia xxx.jl` for setting number of threads.

        # 1). random/greedy initialize a contraction tree.
        tree = _initializetree(code, size_dict, initializer)

        # 2). optimize the `tree` in a inplace manner.
        optimize_tree_sa!(tree, log2_sizes; βs, niters, score)

        # 3). evaluate time-space-readwrite complexities.
        tc, sc, rw = tree_timespace_complexity(tree, log2_sizes)
        @debug "trial $t, time complexity = $tc, space complexity = $sc, read-write complexity = $rw."
        if t == 1 || score(tc, sc, rw) < score(best_tc, best_sc, best_rw)
            best_tree, best_tc, best_sc, best_rw = tree, tc, sc, rw
        end
    end

    ###### Stage 3: postprocessing ######
    @debug "best space complexities = $best_tc, time complexity = $best_sc, read-write complexity $best_rw."
    return NestedEinsum(best_tree, inverse_map)
end

# initialize a contraction tree
function _initializetree(code, size_dict, method)
    if method == :greedy
        labels = _label_dict(code)  # label to int
        return _exprtree(optimize_greedy(code, size_dict; α = 0.0, temperature = 0.0), labels)
    elseif method == :random
        return random_exprtree(code)
    elseif method == :specified
        labels = _label_dict(code)  # label to int
        return _exprtree(code, labels)
    else
        throw(ArgumentError("intializier `$method` is not defined!"))
    end
end

# use simulated annealing to optimize a contraction tree
function optimize_tree_sa!(tree::ExprTree, log2_sizes; βs, niters, score)
    log2rw_weight = log2(score.rw_weight)

    tc, sc, rw = tree_timespace_complexity(tree, log2_sizes)
    @debug "Initial tc = $tc, sc = $sc, rw = $rw"
    
    for β in βs
        for _ = 1:niters
            optimize_subtree!(tree, β, log2_sizes, score.sc_target, score.sc_weight, log2rw_weight)  # single sweep
        end
    end

    tc, sc, rw = tree_timespace_complexity(tree, log2_sizes)
    @debug "After optimization: tc = $tc, sc = $sc, rw = $rw"

    return tree
end

# the time-space-readwrite complexity of a contraction tree
function tree_timespace_complexity(tree::ExprTree, log2_sizes)
    isleaf(tree) && return (-Inf, isempty(labels(tree)) ? 0.0 : sum(i->log2_sizes[i], labels(tree)), -Inf)
    tcl, scl, rwl = tree_timespace_complexity(tree.left, log2_sizes)
    tcr, scr, rwr = tree_timespace_complexity(tree.right, log2_sizes)
    tc, sc, rw = tcscrw(labels(tree.left), labels(tree.right), labels(tree), log2_sizes, true)
    return (fast_log2sumexp2(tc, tcl, tcr), max(sc, scl, scr), fast_log2sumexp2(rw, rwl, rwr))
end

# returns time complexity, space complexity and read-write complexity (0 if `compute_rw` is false)
# `ix1` and `ix2` are vectors of labels for the first and second input tensors.
# `iy` is a vector of labels for the output tensors.
# `log2_sizes` is the log2 size of labels (note labels are integers, we do not need dict to index label sizes).\
@inline function tcscrw(ix1, ix2, iy, log2_sizes::Vector{T}, compute_rw) where T
    l1, l2, l3 = ix1, ix2, iy
    sc1 = (!compute_rw || isempty(l1)) ? zero(T) : sum(i->(@inbounds log2_sizes[i]), l1)
    sc2 = (!compute_rw || isempty(l2)) ? zero(T) : sum(i->(@inbounds log2_sizes[i]), l2)
    sc = isempty(l3) ? zero(T) : sum(i->(@inbounds log2_sizes[i]), l3)
    tc = sc
    # Note: assuming labels in `l1` being unique
    @inbounds for l in l1
        if l ∈ l2 && l ∉ l3
            tc += log2_sizes[l]
        end
    end
    rw = compute_rw ? fast_log2sumexp2(sc, sc1, sc2) : 0.0
    return tc, sc, rw
end

# optimize a contraction tree recursively
function optimize_subtree!(tree, β, log2_sizes, sc_target, sc_weight, log2rw_weight)
    # find appliable local rules, at most 4 rules can be applied.
    # Sometimes, not all rules are applicable because either left or right sibling do not have siblings.
    rst = ruleset(tree)
    if !isempty(rst)
        # propose a random update rule, TODO: can we have a better selector?
        rule = rand(rst)
        optimize_rw = log2rw_weight != -Inf
        # difference in time, space and read-write complexity if the selected rule is applied
        tc0, tc1, dsc, rw0, rw1, subout = tcsc_diff(tree, rule, log2_sizes, optimize_rw)
        dtc = optimize_rw ? fast_log2sumexp2(tc1, log2rw_weight + rw1) - fast_log2sumexp2(tc0, log2rw_weight + rw0) : tc1 - tc0
        sc = _sc(tree, rule, log2_sizes)  # current space complexity

        # update the loss function
        dE = (max(sc, sc+dsc) > sc_target ? sc_weight : 0) * dsc + dtc
        if rand() < exp(-β*dE)  # ACCEPT
            update_tree!(tree, rule, subout)
        end
        for subtree in siblings(tree)  # RECURSE
            optimize_subtree!(subtree, β, log2_sizes, sc_target, sc_weight, log2rw_weight)
        end
    end
end
# if rule ∈ [1, 2], left sibling will be updated, otherwise, right sibling will be updated.
# we need to compute space complexity for current node and the updated sibling, and return the larger one.
_sc(tree, rule, log2_sizes) = max(__sc(tree, log2_sizes), __sc((rule == 1 || rule == 2) ? tree.left : tree.right, log2_sizes))
__sc(tree, log2_sizes) = length(labels(tree))==0 ? 0.0 : sum(l->log2_sizes[l], labels(tree)) # space complexity of current node

@inline function ruleset(tree::ExprTree)
    if isleaf(tree) || (isleaf(tree.left) && isleaf(tree.right))
        return 1:0
    elseif isleaf(tree.right)
        return 1:2
    elseif isleaf(tree.left)
        return 3:4
    else
        return 1:4
    end
end

function tcsc_diff(tree::ExprTree, rule, log2_sizes, optimize_rw)
    if rule == 1 # (a,b), c -> (a,c),b
        return abcacb(labels(tree.left.left), labels(tree.left.right), labels(tree.left), labels(tree.right), labels(tree), log2_sizes, optimize_rw)
    elseif rule == 2 # (a,b), c -> (c,b),a
        return abcacb(labels(tree.left.right), labels(tree.left.left), labels(tree.left), labels(tree.right), labels(tree), log2_sizes, optimize_rw)
    elseif rule == 3 # a,(b,c) -> b,(a,c)
        return abcacb(labels(tree.right.right), labels(tree.right.left), labels(tree.right), labels(tree.left), labels(tree), log2_sizes, optimize_rw)
    else  # a,(b,c) -> c,(b,a)
        return abcacb(labels(tree.right.left), labels(tree.right.right), labels(tree.right), labels(tree.left), labels(tree), log2_sizes, optimize_rw)
    end
end

# compute the time complexity, space complexity and read-write complexity information for the contraction update rule "((a,b),c) -> ((a,c),b)"
function abcacb(a, b, ab, c, d, log2_sizes, optimize_rw)
    tc0, sc0, rw0 = _tcsc_merge(a, b, ab, c, d, log2_sizes, optimize_rw)
    ac = Int[]  # labels for contraction result of (a, c)
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
    tc1, sc1, rw1 = _tcsc_merge(a, c, ac, b, d, log2_sizes, optimize_rw)
    return tc0, tc1, sc1-sc0, rw0, rw1, ac
end

# compute complexity for a two-step contraction: (a, b) -> ab, (ab, c) -> d
function _tcsc_merge(a, b, ab, c, d, log2_sizes, optimize_rw)
    tcl, scl, rwl = tcscrw(a, b, ab, log2_sizes, optimize_rw)  # this is correct
    tcr, scr, rwr = tcscrw(ab, c, d, log2_sizes, optimize_rw)
    fast_log2sumexp2(tcl, tcr), max(scl, scr), (optimize_rw ? fast_log2sumexp2(rwl, rwr) : 0.0)
end

# apply the update rule
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

# map labels to integers.
_label_dict(code) = _label_dict(getixsv(code), getiyv(code))
function _label_dict(ixsv::AbstractVector{<:AbstractVector{LT}}, iyv::AbstractVector{LT}) where LT
    v = unique(vcat(ixsv..., iyv))
    labels = Dict{LT,Int}(zip(v, 1:length(v)))
    return labels
end

# construct the contraction tree recursively from a nested einsum.
function ExprTree(code::NestedEinsum)
    _exprtree(code, _label_dict(code))
end
function _exprtree(code::NestedEinsum, labels)
    @assert length(code.args) == 2 "einsum contraction not in the binary form, got number of arguments: $(length(code.args))"
    ExprTree(map(enumerate(code.args)) do (i,arg)
        if isleaf(arg)  # leaf nodes
            ExprTree(ExprInfo(getindex.(Ref(labels), getixsv(code.eins)[i]), arg.tensorindex))
        else
            res = _exprtree(arg, labels)
        end
    end..., ExprInfo(Int[labels[i] for i=getiyv(code.eins)]))
end

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

