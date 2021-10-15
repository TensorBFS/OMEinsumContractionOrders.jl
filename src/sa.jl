using SparseArrays
using Base: RefValue
using BetterExp

export SABipartite, optimize_sa

Base.@kwdef struct SABipartite{RT,BT}
    sc_target::RT = 25
    ntrials::Int = 50  # number of trials
    βs::BT = 0.1:0.2:15.0  # temperatures
    niters::Int = 1000  # number of iterations in each temperature
    max_group_size::Int = 40
    # configure greedy algorithm
    greedy_method = OMEinsum.MinSpaceOut()
    greedy_nrepeat::Int = 10
    initializer::Symbol = :random
end

struct PartitionState
    config::Vector{Int}
    loss::RefValue{Float64}

    group_sizes::Vector{Int}     # group sizes
    group_scs::Vector{Float64}   # space complexities
    group_degrees::Matrix{Int}   # degree of edges
end

function partition_state(adj, group, config, log2_sizes)
    group1 = [group[i] for (i,c) in enumerate(config) if c==1]
    group2 = [group[i] for (i,c) in enumerate(config) if c==2]
    group_sizes = [length(group1), length(group2)]   # group size in terms of number of vertices.
    group_scs = [group_sc(adj, group1, log2_sizes), group_sc(adj, group2, log2_sizes)]  # space complexity of each group.
    group_degrees = hcat(sum(adj[group1,:]; dims=1)', sum(adj[group2,:]; dims=1)')
    loss = compute_loss(group_scs..., group_sizes...)
    return PartitionState(config, Ref(loss), group_sizes, group_scs, group_degrees)
end

function bipartite_sc(bipartiter::SABipartite, adj::SparseMatrixCSC, vertices, log2_sizes)
    degrees_all = sum(adj, dims=1)
    adjt = SparseMatrixCSC(adj')
    config = _initialize(bipartiter.initializer,adj, vertices, log2_sizes)
    best = partition_state(adj, vertices, config, log2_sizes)  # this is the `state` of current partition.

    for _ = 1:bipartiter.ntrials
        config = _initialize(bipartiter.initializer,adj, vertices, log2_sizes)
        state = partition_state(adj, vertices, config, log2_sizes)  # this is the `state` of current partition.
        if state.group_sizes[1]==0 || state.group_sizes[2] == 0
            continue
        end

        @inbounds for β in bipartiter.βs, iter = 1:bipartiter.niters
            idxi = rand(1:length(vertices))
            ti = state.config[idxi]
            state.group_sizes[ti] <= 1 && continue

            sc_ti, sc_tinew = space_complexity_singlestep_update(state, adjt, degrees_all, log2_sizes, vertices, idxi)
            newloss = compute_loss(sc_ti, sc_tinew, state.group_sizes[ti]-1, state.group_sizes[3-ti]+1)
            sc_ti0, sc_tinew0 = state.group_scs[ti], state.group_scs[3-ti]
            accept = if max(sc_ti0, sc_tinew0) <= bipartiter.sc_target
                max(sc_ti, sc_tinew) <= bipartiter.sc_target && (rand() < BetterExp.exp2(β*(state.loss[] - newloss)))
            else
                rand() < BetterExp.exp2(-β*(max(sc_ti, sc_tinew) - max(sc_ti0, sc_tinew0)))
            #elseif sc_ti > sc_ti0 && sc_tinew > sc_tinew0
            #    false
            #elseif sc_ti <= sc_ti0 && sc_tinew <= sc_tinew0
            #    true
            #elseif sc_tinew >= sc_tinew0 && sc_tinew <  bipartiter.sc_target
            #    true
            #elseif sc_ti >= sc_ti0 && sc_ti < bipartiter.sc_target
            #    true
            #elseif max(sc_ti, sc_tinew) < max(sc_ti0, sc_tinew0)
            #    true
            #else
            #    false
            end
            accept && update_state!(state, adjt, vertices, idxi, sc_ti, sc_tinew, newloss)
        end
        tc, sc1, sc2 = timespace_complexity_singlestep(state.config, adj, vertices, log2_sizes)
        @assert state.group_scs ≈ [sc1, sc2]  # sanity check
        if maximum(state.group_scs) <= max(bipartiter.sc_target, maximum(best.group_scs)) && (maximum(best.group_scs) >= bipartiter.sc_target || state.loss[] < best.loss[])
            best = state
        end
    end
    best_tc, = timespace_complexity_singlestep(best.config, adj, vertices, log2_sizes)
    @debug "best loss = $(round(best.loss[]; digits=3)) space complexities = $(best.group_scs) time complexity = $(best_tc) groups_sizes = $(best.group_sizes)"
    if maximum(best.group_scs) > bipartiter.sc_target
        @warn "target space complexity not found, got: $(maximum(best.group_scs)), with time complexity $best_tc."
    end
    return vertices[findall(==(1), best.config)], vertices[findall(==(2), best.config)]
end

function timespace_complexity_singlestep(config, adj, group, log2_sizes)
    g1 = group[findall(==(1), config)]
    g2 = group[findall(==(2), config)]
    d1 = sum(adj[g1,:], dims=1)
    d2 = sum(adj[g2,:], dims=1)
    dall = sum(adj, dims=1)
    sc1 = sum(i->(d1[i]!=0 && d1[i]!=dall[i] ? Float64(log2_sizes[i]) : 0.0), 1:size(adj,2))
    sc2 = sum(i->(d2[i]!=0 && d2[i]!=dall[i] ? Float64(log2_sizes[i]) : 0.0), 1:size(adj,2))
    tc = sum(i->((d2[i]!=0 || d1[i]!=0) && (d2[i]!=dall[i] && d1[i]!=dall[i]) ? Float64(log2_sizes[i]) : 0.0), 1:size(adj,2))
    return tc, sc1, sc2
end

function space_complexity_singlestep_update(state, adjt, degrees_all, log2_sizes, group, idxi)
    @inbounds begin
        vertex = group[idxi]
        ti = state.config[idxi]
        tinew = 3-ti
        δsc_ti = δsc(-, view(state.group_degrees, :, ti), adjt, vertex, degrees_all, log2_sizes)
        δsc_tinew = δsc(+, view(state.group_degrees, :, tinew), adjt, vertex, degrees_all, log2_sizes)
        sc_ti = state.group_scs[ti] + δsc_ti
        sc_tinew = state.group_scs[tinew] + δsc_tinew
    end
    return sc_ti, sc_tinew
end

@inline function δsc(f, group_degrees, adjt, vertex, degrees_all, log2_sizes)
    res = 0.0
    @inbounds for k in nzrange(adjt, vertex)
        i = adjt.rowval[k]
        d0 = group_degrees[i]
        D = degrees_all[i]
        d = f(d0, adjt.nzval[k])
        if d0 == D || d0 == 0     # absent
            if d != D && d != 0   # absent
                res += Float64(log2_sizes[i])
            end
        else                      # not absent
            if d == D || d == 0   # absent
                res -= Float64(log2_sizes[i])
            end
        end
    end
    return res
end


@inline function compute_loss(sc1, sc2, gs1, gs2)
    small = min(gs1, gs2)
    max(sc1, sc2) / small * (gs1 + gs2)
end

function update_state!(state, adjt, group, idxi, sc_ti, sc_tinew, newloss)
    @inbounds begin
        ti = state.config[idxi]
        tinew = 3-ti
        state.group_scs[tinew] = sc_tinew
        state.group_scs[ti] = sc_ti
        state.config[idxi] = tinew
        state.group_sizes[ti] -= 1
        state.group_sizes[tinew] += 1
        for i = nzrange(adjt,group[idxi])
            state.group_degrees[adjt.rowval[i], ti] -= adjt.nzval[i]
            state.group_degrees[adjt.rowval[i], tinew] += adjt.nzval[i]
        end
        state.loss[] = newloss
    end
    return state
end

_initialize(method, adj, vertices, log2_sizes) = if method == :random
    initialize_random(length(vertices))
elseif method == :greedy
    initialize_greedy(adj, vertices, log2_sizes)
else
    error("initializer not implemented: `$method`")
end

function initialize_greedy(adj, vertices, log2_sizes)
    adjt = SparseMatrixCSC(adj')
    indegrees = sum(adj[vertices,:], dims=1)
    all = sum(adj, dims=1)
    openedges = findall(i->all[i] > indegrees[i] > 0, 1:size(adj, 2))
    v2e = Dict{Int,Vector{Int}}()
    for v=1:size(adj, 1)
        v2e[v] = adjt.rowval[nzrange(adjt, v)]
    end
    incidence_list = OMEinsum.IncidenceList(v2e; openedges=openedges)
    log2_edge_sizes = Dict([i=>log2_sizes[i] for i=1:length(log2_sizes)])
    # nrepeat=3 because there are overheads
    tree, _, _ = OMEinsum.tree_greedy(incidence_list, log2_edge_sizes; method=OMEinsum.MinSpaceOut(), nrepeat=3)

    # build configuration from the tree
    res = ones(Int, size(adj, 1))
    res[get_vertices!(Int[],tree.right)] .= 2
    return res[vertices]
end
initialize_random(n::Int) = [rand() < 0.5 ? 1 : 2 for _ = 1:n]

function get_vertices!(out, tree)
    if tree isa Integer
        push!(out, tree)
    else
        get_vertices!(out, tree.left)
        get_vertices!(out, tree.right)
    end
    return out
end



## interfaces
"""
    optimize_sa(code, size_dict; sc_target, max_group_size=40, βs=0.1:0.2:15.0, niters=1000, ntrials=50,
            greedy_method=OMEinsum.MinSpaceOut(), greedy_nrepeat=10, initializer=:random)

Optimize the einsum code contraction order using the Simulated Annealing + Greedy approach.
This program first recursively cuts the tensors into several groups using simulated annealing,
with maximum group size specifed by `max_group_size` and maximum space complexity specified by `sc_target`,
Then finds the contraction order inside each group with the greedy search algorithm. Other arguments are

* `size_dict`, a dictionary that specifies leg dimensions,
* `sc_target` is the target space complexity, defined as `log2(number of elements in the largest tensor)`,
* `max_group_size` is the maximum size that allowed to used greedy search,
* `βs` is a list of inverse temperature `1/T`,
* `niters` is the number of iteration in each temperature,
* `ntrials` is the number of repetition (with different random seeds),
* `greedy_method` and `greedy_nrepeat` are for configuring the greedy method,
* `initializer`, the partition configuration initializer, one can choose `:random` or `:greedy` (slow but better).

### References
* [Hyper-optimized tensor network contraction](https://arxiv.org/abs/2002.01935)
"""
function optimize_sa(@nospecialize(code::EinCode), size_dict; sc_target, max_group_size=40,
             βs=0.01:0.02:15.0, niters=1000, ntrials=50, greedy_method=OMEinsum.MinSpaceOut(), greedy_nrepeat=10,
             initializer=:random)
    bipartiter = SABipartite(; sc_target=sc_target, βs=βs, niters=niters, ntrials=ntrials, greedy_method=greedy_method, greedy_nrepeat=greedy_nrepeat, max_group_size=max_group_size, initializer=initializer)
    
    recursive_bipartite_optimize(bipartiter, code, size_dict)
end
