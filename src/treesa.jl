using OMEinsum.ContractionOrder: ContractionTree, log2sumexp2
using Base.Threads

export optimize_tree

"""
    TreeSA{RT,IT,GM}
    TreeSA(; sc_target=20, βs=collect(0.01:0.05:15), ntrials=10, niters=50,
        sc_weight=1.0, rw_weight=0.2, initializer=:greedy, greedy_config=GreedyMethod(; nrepeat=1))

Optimize the einsum contraction pattern using the simulated annealing on tensor expression tree.

* `sc_target` is the target space complexity,
* `ntrials`, `βs` and `niters` are annealing parameters, doing `ntrials` indepedent annealings, each has inverse tempteratures specified by `βs`, in each temperature, do `niters` updates of the tree.
* `sc_weight` is the relative importance factor of space complexity in the loss compared with the time complexity.
* `rw_weight` is the relative importance factor of memory read and write in the loss compared with the time complexity.
* `initializer` specifies how to determine the initial configuration, it can be `:greedy` or `:random`. If it is using `:greedy` method to generate the initial configuration, it also uses two extra arguments `greedy_method` and `greedy_nrepeat`.
* `nslices` is the number of sliced legs, default is 0.

### References
* [Recursive Multi-Tensor Contraction for XEB Verification of Quantum Circuits](https://arxiv.org/abs/2108.05665)
"""
Base.@kwdef struct TreeSA{RT,IT,GM} <: CodeOptimizer
    sc_target::RT = 20
    βs::IT = 0.01:0.05:15
    ntrials::Int = 10
    niters::Int = 50
    sc_weight::Float64 = 1.0
    rw_weight::Float64 = 0.2
    initializer::Symbol = :greedy
    nslices::Int = 0
    # configure greedy method
    greedy_config::GM = GreedyMethod(nrepeat=1)
end

struct ExprInfo
    out_dims::Vector{Int}
    tensorid::Int
end
ExprInfo(out_dims::Vector{Int}) = ExprInfo(out_dims, -1)

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
OMEinsum.isleaf(expr::ExprTree) = !isdefined(expr, :left)
Base.show(io::IO, expr::ExprTree) = print_expr(io, expr, 0)
Base.show(io::IO, ::MIME"text/plain", expr::ExprTree) = show(io, expr)
siblings(t::ExprTree) = isleaf(t) ? ExprTree[] : ExprTree[t.left, t.right]
Base.copy(t::ExprTree) = isleaf(t) ? ExprTree(t.info) : ExprTree(copy(t.left), copy(t.right), copy(t.info))
Base.copy(info::ExprInfo) = ExprInfo(copy(info.out_dims), info.tensorid)
labels(t::ExprTree) = t.info.out_dims
maxlabel(t::ExprTree) = isleaf(t) ? maximum(isempty(labels(t)) ? 0 : labels(t)) : max(isempty(labels(t)) ? 0 : maximum(labels(t)), maxlabel(t.left), maxlabel(t.right))
Base.:(==)(t1::ExprTree, t2::ExprTree) = _equal(t1, t2)
Base.:(==)(t1::ExprInfo, t2::ExprInfo) = _equal(t1.out_dims, t2.out_dims) && t1.tensorid == t2.tensorid
function _equal(t1::ExprTree, t2::ExprTree)
    isleaf(t1) != isleaf(t2) && return false
    isleaf(t1) ? t1.info == t2.info : _equal(t1.left, t2.left) && _equal(t1.right, t2.right) && t1.info == t2.info
end
_equal(t1::Vector, t2::Vector) = Set(t1) == Set(t2)

"""
    optimize_tree(code, size_dict; sc_target=20, βs=0.1:0.1:10, ntrials=2, niters=100, sc_weight=1.0, rw_weight=0.2, initializer=:greedy, greedy_method=OMEinsum.MinSpaceOut(), greedy_nrepeat=1)

Optimize the einsum contraction pattern specified by `code`, and edge sizes specified by `size_dict`. Key word arguments are
Check the docstring of `TreeSA` for detailed explaination of other input arguments.
"""
function optimize_tree(code, size_dict; nslices::Int=0, sc_target=20, βs=0.1:0.1:10, ntrials=20, niters=100, sc_weight=1.0, rw_weight=0.2, initializer=:greedy, greedy_method=OMEinsum.MinSpaceOut(), greedy_nrepeat=1)
    flatten_code = OMEinsum.flatten(code)
    ninputs = length(OMEinsum.getixs(flatten_code))
    if ninputs <= 2
        return NestedEinsum(ntuple(i->i, ninputs), DynamicEinCode(flatten_code))
    end
    labels = _label_dict(flatten_code)  # label to int
    inverse_map = Dict([v=>k for (k,v) in labels])
    log2_sizes = [log2.(size_dict[inverse_map[i]]) for i=1:length(labels)]
    if ntrials <= 0
        best_tree = _initializetree(code, size_dict, initializer; greedy_method=greedy_method, greedy_nrepeat=greedy_nrepeat)
        return NestedEinsum(best_tree, inverse_map)
    end
    trees, tcs, scs, rws, slicers = Vector{ExprTree}(undef, ntrials), zeros(ntrials), zeros(ntrials), zeros(ntrials), Vector{Slicer}(undef, ntrials)
    @threads for t = 1:ntrials
        tree = _initializetree(code, size_dict, initializer; greedy_method=greedy_method, greedy_nrepeat=greedy_nrepeat)
        slicer = Slicer(log2_sizes, nslices)
        optimize_tree_sa!(tree, log2_sizes, slicer; sc_target=sc_target, βs=βs, niters=niters, sc_weight=sc_weight, rw_weight=rw_weight)
        tc, sc, rw = tree_timespace_complexity(tree, log2_sizes)
        @debug "trial $t, time complexity = $tc, space complexity = $sc, read-write complexity = $rw."
        trees[t], tcs[t], scs[t], rws[t], slicers[t] = tree, tc, sc, rw, slicer
    end
    best_tree, best_tc, best_sc, best_rw, best_slicer = first(trees), first(tcs), first(scs), first(rws), first(slicers)
    for i=2:ntrials
        if scs[i] < best_sc || (scs[i] == best_sc && exp2(tcs[i]) + rw_weight * exp2(rws[i]) < exp2(best_tc) + rw_weight * exp2(rws[i]))
            best_tree, best_tc, best_sc, best_rw, best_slicer = trees[i], tcs[i], scs[i], rws[i], slicers[i]
        end
    end
    @debug "best space complexities = $best_tc, time complexity = $best_sc, read-write complexity $best_rw."
    if best_sc > sc_target
        @warn "target space complexity not found, got: $best_sc, with time complexity $best_tc, read-write complexity $best_rw."
    end
    return SlicedEinsum(Slicing(best_slicer, inverse_map), NestedEinsum(best_tree, inverse_map))
end

function _initializetree(code, size_dict, method; greedy_method, greedy_nrepeat)
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

function optimize_tree_sa!(tree::ExprTree, log2_sizes, slicer::Slicer; βs, niters, sc_target, sc_weight, rw_weight)
    @assert rw_weight >= 0
    @assert sc_weight >= 0
    log2rw_weight = log2(rw_weight)
    for β in βs
        @debug begin
            tc, sc, rw = tree_timespace_complexity(tree, log2_sizes)
            "β = $β, tc = $tc, sc = $sc, rw = $rw"
        end
        if slicer.max_size > 0
            # find legs that reduce the dimension the most
            scs, lbs = Float64[], Vector{Int}[]
            tensor_sizes!(tree, slicer.log2_sizes, scs, lbs)
            best_labels = _best_labels(scs, lbs)

            best_not_sliced_labels = filter(x->!haskey(slicer.legs, x), best_labels)
            if !isempty(best_not_sliced_labels)
                best_not_sliced_label = rand(best_not_sliced_labels)
                if length(slicer) < slicer.max_size
                    push!(slicer, best_not_sliced_label)
                else
                    #worst_sliced_labels = filter(x->haskey(slicer.legs, x), setdiff(log2_sizes, best_labels))
                    legs = collect(keys(slicer.legs))
                    score = [count(==(l), best_labels) for l in legs]
                    replace!(slicer, legs[argmin(score)]=>best_not_sliced_label)
                end
            end
            @debug begin
                tc, sc, rw = tree_timespace_complexity(tree, log2_sizes)
                "after slicing: β = $β, tc = $tc, sc = $sc, rw = $rw"
            end
        end
        for _ = 1:niters
            optimize_subtree!(tree, β, slicer.log2_sizes, sc_target, sc_weight, log2rw_weight)  # single sweep
        end
    end
    return tree, slicer
end

function _best_labels(scs, lbs)
    max_sc = maximum(scs)
    return vcat(lbs[scs .> max_sc-0.99]...)
end

function tensor_sizes!(tree::ExprTree, log2_sizes, scs, lbs)
    sc = isempty(labels(tree)) ? 0.0 : sum(i->log2_sizes[i], labels(tree))
    push!(scs, sc)
    push!(lbs, labels(tree))
    isleaf(tree) && return
    tensor_sizes!(tree.left, log2_sizes, scs, lbs)
    tensor_sizes!(tree.right, log2_sizes, scs, lbs)
end

function tree_timespace_complexity(tree::ExprTree, log2_sizes)
    isleaf(tree) && return (-Inf, isempty(labels(tree)) ? 0.0 : sum(i->log2_sizes[i], labels(tree)), -Inf)
    tcl, scl, rwl = tree_timespace_complexity(tree.left, log2_sizes)
    tcr, scr, rwr = tree_timespace_complexity(tree.right, log2_sizes)
    tc, sc, rw = tcscrw(labels(tree.left), labels(tree.right), labels(tree), log2_sizes, true)
    return (fast_log2sumexp2(tc, tcl, tcr), max(sc, scl, scr), fast_log2sumexp2(rw, rwl, rwr))
end
@inline function tcscrw(ix1, ix2, iy, log2_sizes::Vector{T}, optimize_rw) where T
    l1, l2, l3 = ix1, ix2, iy
    sc1 = (!optimize_rw || isempty(l1)) ? zero(T) : sum(i->(@inbounds log2_sizes[i]), l1)
    sc2 = (!optimize_rw || isempty(l2)) ? zero(T) : sum(i->(@inbounds log2_sizes[i]), l2)
    sc = isempty(l3) ? zero(T) : sum(i->(@inbounds log2_sizes[i]), l3)
    tc = sc
    # Note: assuming labels in `l1` being unique
    @inbounds for l in l1
        if l ∈ l2 && l ∉ l3
            tc += log2_sizes[l]
        end
    end
    rw = optimize_rw ? fast_log2sumexp2(sc, sc1, sc2) : 0.0
    return tc, sc, rw
end

function random_exprtree(code::EinCode)
    labels = _label_dict(code)
    return random_exprtree([Int[labels[l] for l in ix] for ix in getixsv(code)], Int[labels[l] for l in getiyv(code)], length(labels))
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

function optimize_subtree!(tree, β, log2_sizes, sc_target, sc_weight, log2rw_weight)
    rst = ruleset(tree)
    if !isempty(rst)
        rule = rand(rst)
        optimize_rw = log2rw_weight != -Inf
        tc0, tc1, dsc, rw0, rw1, subout = tcsc_diff(tree, rule, log2_sizes, optimize_rw)
        dtc = optimize_rw ? fast_log2sumexp2(tc1, log2rw_weight + rw1) - fast_log2sumexp2(tc0, log2rw_weight + rw0) : tc1 - tc0
        sc = _sc(tree, rule, log2_sizes)
        dE = (max(sc, sc+dsc) > sc_target ? sc_weight : 0) * dsc + dtc
        if rand() < exp(-β*dE)
            update_tree!(tree, rule, subout)
        end
        for subtree in siblings(tree)
            optimize_subtree!(subtree, β, log2_sizes, sc_target, sc_weight, log2rw_weight)
        end
    end
end
_sc(tree, rule, log2_sizes) = max(__sc(tree, log2_sizes), __sc((rule == 1 || rule == 2) ? tree.left : tree.right, log2_sizes))
__sc(tree, log2_sizes) = length(labels(tree))==0 ? 0.0 : sum(l->log2_sizes[l], labels(tree))

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

function abcacb(a, b, ab, c, d, log2_sizes, optimize_rw)
    tc0, sc0, rw0, ab0 = _tcsc_merge(a, b, ab, c, d, log2_sizes, optimize_rw)
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
    tc1, sc1, rw1, ab1 = _tcsc_merge(a, c, ac, b, d, log2_sizes, optimize_rw)
    return tc0, tc1, sc1-sc0, rw0, rw1, ab1  # Note: this tc diff does not make much sense
end

function _tcsc_merge(a, b, ab, c, d, log2_sizes, optimize_rw)
    tcl, scl, rwl = tcscrw(a, b, ab, log2_sizes, optimize_rw)  # this is correct
    tcr, scr, rwr = tcscrw(ab, c, d, log2_sizes, optimize_rw)
    fast_log2sumexp2(tcl, tcr), max(scl, scr), (optimize_rw ? fast_log2sumexp2(rwl, rwr) : 0.0), ab
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
function _label_dict(code::EinCode)
    LT = OMEinsum.labeltype(code)
    ixsv, iyv = getixsv(code), getiyv(code)
    v = unique(vcat(ixsv..., iyv))
    labels = Dict{LT,Int}(zip(v, 1:length(v)))
    return labels
end

function ExprTree(code::NestedEinsum)
    _exprtree(code, _label_dict(OMEinsum.flatten(code)))
end
function _exprtree(code::NestedEinsum, labels)
    @assert length(code.args) == 2
    ExprTree(map(enumerate(code.args)) do (i,arg)
        if isleaf(arg)
            ExprTree(ExprInfo(Int[labels[i] for i=OMEinsum.getixs(code.eins)[i]], arg.tensorindex))
        else
            res = _exprtree(arg, labels)
        end
    end..., ExprInfo(Int[labels[i] for i=OMEinsum.getiy(code.eins)]))
end

OMEinsum.NestedEinsum(expr::ExprTree) = _nestedeinsum(expr, 1:maxlabel(expr))
OMEinsum.NestedEinsum(expr::ExprTree, labelmap) = _nestedeinsum(expr, labelmap)
function _nestedeinsum(tree::ExprTree, lbs)
    isleaf(tree) && return tree.info.tensorid
    eins = EinCode([getindex.(Ref(lbs), labels(tree.left)), getindex.(Ref(lbs), labels(tree.right))], getindex.(Ref(lbs), labels(tree)))
    NestedEinsum((_nestedeinsum(tree.left, lbs), _nestedeinsum(tree.right, lbs)), eins)
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

