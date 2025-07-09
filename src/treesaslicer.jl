
############# Slicer ######################
struct Slicer
    log2_sizes::Vector{Float64}   # the size dict after slicing
    legs::Dict{Int,Float64}   # sliced leg and its original size
    fixed_slices::Vector{Int}     # number of fixed legs
end
function Slicer(log2_sizes::AbstractVector{Float64}, fixed_slices::AbstractVector)
    slicer = Slicer(collect(log2_sizes), Dict{Int,Float64}(), collect(Int,fixed_slices))
    for l in fixed_slices
        push!(slicer, l)
    end
    return slicer
end
Base.length(s::Slicer) = length(s.legs)
function Base.replace!(slicer::Slicer, pair::Pair)
    worst, best = pair
    @assert worst ∉ slicer.fixed_slices
    @assert haskey(slicer.legs, worst)
    @assert !haskey(slicer.legs, best)
    slicer.log2_sizes[worst] = slicer.legs[worst]       # restore worst size
    slicer.legs[best] = slicer.log2_sizes[best]  # add best to legs
    slicer.log2_sizes[best] = 0.0
    delete!(slicer.legs, worst)                  # remove worst from legs
    return slicer
end

function Base.push!(slicer::Slicer, best)
    @assert !haskey(slicer.legs, best)
    slicer.legs[best] = slicer.log2_sizes[best]  # add best to legs
    slicer.log2_sizes[best] = 0.0
    return slicer
end

function Base.delete!(slicer::Slicer, worst)
    @assert haskey(slicer.legs, worst)
    slicer.log2_sizes[worst] = slicer.legs[worst]
    delete!(slicer.legs, worst)
    return slicer
end

# convert the slicer to a vector of sliced labels
function get_slices(s::Slicer, inverse_map::Dict{Int,LT}) where LT
    # we want to keep the order of input fixed slices!
    LT[[inverse_map[l] for l in s.fixed_slices]..., [inverse_map[l] for (l, sz) in s.legs if l ∉ s.fixed_slices]...]
end

"""
    TreeSASlicer{RT, IT, LT} <: CodeSlicer

A structure for configuring the Tree Simulated Annealing (TreeSA) slicing algorithm.

# Fields
- `sc_target::RT`: The target space complexity. Default is 30.
- `ntrials`, `βs` and `niters` are annealing parameters, doing `ntrials` indepedent annealings, each has inverse tempteratures specified by `βs`, in each temperature, do `niters` updates of the tree.
- `sc_weight::Float64`: The weight for space complexity in the optimization. Default is 1.0.
- `rw_weight::Float64`: The weight for read-write complexity in the optimization. Default is 0.2.
- `fixed_slices::Vector{LT}`: A vector of fixed slices that should not be altered. Default is an empty vector.
- `optimization_ratio::Float64`: The ratio used for optimization. Default is 2.0.

# References
- [Recursive Multi-Tensor Contraction for XEB Verification of Quantum Circuits](https://arxiv.org/abs/2108.05665)
"""
Base.@kwdef struct TreeSASlicer{RT,IT,LT} <: CodeSlicer
    sc_target::RT = 30
    βs::IT = 14:0.05:15
    ntrials::Int = 10
    niters::Int = 10
    sc_weight::Float64 = 1.0
    rw_weight::Float64 = 0.2
    fixed_slices::Vector{LT} = []
    optimization_ratio::Float64 = 2.0
end

function slice_tree(code::NestedEinsum, size_dict::Dict{LT,Int}; sc_target=30, βs=14:0.05:15, ntrials=10, niters=10, sc_weight=1.0, rw_weight=0.2, fixed_slices=LT[], optimization_ratio=2.0) where LT
    ixs, iy = getixsv(code), getiyv(code)
    ninputs = length(ixs)
    if ninputs <= 2
        return SlicedEinsum(LT[], code)
    end

    ###### Stage 1: preprocessing ######
    labels = _label_dict(ixs, iy)  # map labels to integers
    inverse_map = Dict([v=>k for (k,v) in labels])  # the inverse transformation, map integers to labels
    log2_sizes = [log2.(size_dict[inverse_map[i]]) for i=1:length(labels)]   # use `log2` sizes in computing time
    labels = _label_dict(code)

    if ntrials <= 0
        return SlicedEinsum(fixed_slices, code)
    end

    ###### Stage 2: computing ######
    trees, tcs, scs, rws, slicers = Vector{ExprTree}(undef, ntrials), zeros(ntrials), zeros(ntrials), zeros(ntrials), Vector{Slicer}(undef, ntrials)
    @threads for t = 1:ntrials  # multi-threading on different trials, use `julia -t 5 xxx.jl` for setting number of threads.
        tree = _exprtree(code, labels)
        slicer = Slicer(log2_sizes, Int[labels[l] for l in fixed_slices])
        
        treesa_slice!(tree, log2_sizes, slicer; βs=βs, niters=niters, sc_target=sc_target, sc_weight=sc_weight, rw_weight=rw_weight, optimization_ratio=optimization_ratio)
        tc, sc, rw = tree_timespace_complexity(tree, log2_sizes)
        @debug "trial $t, time complexity = $tc, space complexity = $sc, read-write complexity = $rw."
        trees[t], tcs[t], scs[t], rws[t], slicers[t] = tree, tc, sc, rw, slicer
    end

    ###### Stage 3: postprocessing ######
    best_tree, best_tc, best_sc, best_rw, best_slicer = first(trees), first(tcs), first(scs), first(rws), first(slicers)
    for i=2:ntrials
        if exp2(tcs[i]) + rw_weight * exp2(rws[i]) < exp2(best_tc) + rw_weight * exp2(rws[i])
            best_tree, best_tc, best_sc, best_rw, best_slicer = trees[i], tcs[i], scs[i], rws[i], slicers[i]
        end
    end
    @debug "best space complexities = $best_tc, time complexity = $best_sc, read-write complexity $best_rw."

    # returns a sliced einsum we need to map the sliced dimensions back from integers to labels.
    return SlicedEinsum(get_slices(best_slicer, inverse_map), NestedEinsum(best_tree, inverse_map))
end

function treesa_slice!(tree::ExprTree, log2_sizes, slicer::Slicer; βs, niters, sc_target, sc_weight, rw_weight, optimization_ratio)
    @assert rw_weight >= 0
    @assert sc_weight >= 0
    log2rw_weight = log2(rw_weight)
    tc, sc, rw = tree_timespace_complexity(tree, slicer.log2_sizes)
    @debug "Initial tc = $tc, sc = $sc, rw = $rw"

    # determine the number of iterations for slicing
    optimization_length = Int(ceil(optimization_ratio * (sc - sc_target)))
    slicing_loop = 0

    while slicing_loop < optimization_length || sc > sc_target

        ###### Stage 1: add one slice at each loop  ######
        # 1). find legs that reduce the dimension the most
        scs, lbs = Float64[], Vector{Int}[]
        # space complexities and labels of all intermediate tensors
        tensor_sizes!(tree, slicer.log2_sizes, scs, lbs)
        # the set of (intermediate) tensor labels that producing maximum space complexity
        best_labels = _best_labels(scs, lbs)

        # 2). slice the best not sliced label (it must appear in largest tensors)
        best_not_sliced_labels = filter(x->!haskey(slicer.legs, x), best_labels)
        if !isempty(best_not_sliced_labels)
            #best_not_sliced_label = rand(best_not_sliced_labels)  # random or best
            best_not_sliced_label = best_not_sliced_labels[findmax(l->count(==(l), best_not_sliced_labels), best_not_sliced_labels)[2]]
            if sc > sc_target  # if has not reached maximum number of slices, add one slice
                push!(slicer, best_not_sliced_label)
            elseif length(slicer.legs) > 0                                 # otherwise replace one slice
                legs = [l for l in keys(slicer.legs) if l ∉ slicer.fixed_slices]  # only slice over not fixed legs
                if !isempty(legs)
                    score = [count(==(l), best_labels) for l in legs]
                    # replace!(slicer, legs[argmin(score)]=>best_not_sliced_label)
                    delete!(slicer, legs[argmin(score)])
                end
            end
        end

        ###### Stage 2: refine the tree with the revised slices  ######
        for β in βs
            for _ = 1:niters
                optimize_subtree!(tree, β, slicer.log2_sizes, sc_target, sc_weight, log2rw_weight)  # single sweep
            end
        end
        tc, sc, rw = tree_timespace_complexity(tree, slicer.log2_sizes)
        @debug "After optimization: tc = $tc, sc = $sc, rw = $rw, slicing_loop = $slicing_loop, number of slices = $(length(slicer.legs))"
        slicing_loop += 1
    end
    return tree, slicer
end


# here "best" means giving maximum space complexity
function _best_labels(scs, lbs)
    max_sc = maximum(scs)
    return vcat(lbs[scs .> max_sc-0.99]...)
end

# find tensor sizes and their corresponding labels of all intermediate tensors
function tensor_sizes!(tree::ExprTree, log2_sizes, scs, lbs)
    sc = isempty(labels(tree)) ? 0.0 : sum(i->log2_sizes[i], labels(tree))
    push!(scs, sc)
    push!(lbs, labels(tree))
    isleaf(tree) && return
    tensor_sizes!(tree.left, log2_sizes, scs, lbs)
    tensor_sizes!(tree.right, log2_sizes, scs, lbs)
end