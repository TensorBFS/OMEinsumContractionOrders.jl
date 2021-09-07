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
end

struct PartitionState
    config::Vector{Int}
    loss::RefValue{Float64}

    group_sizes::Vector{Int}     # group sizes
    group_scs::Vector{Float64}   # space complexities
    group_degrees::Matrix{Int}   # degree of edges
end

null_state() = PartitionState(Int[],  Ref(NaN), Int[], Float64[], zeros(Int, 0, 0))
isnull(state::PartitionState) = isnan(state.loss[])

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
    best = null_state() # it stores the configurations that satisfies the sc constraint, while minimizes cutsize, which is tc in bi-section
    degrees_all = sum(adj, dims=1)

    for _ = 1:bipartiter.ntrials
        config = [rand() < 0.5 ? 1 : 2 for i = 1:length(vertices)]
        state = partition_state(adj, vertices, config, log2_sizes)  # this is the `state` of current partition.
        if state.group_sizes[1]==0 || state.group_sizes[2] == 0
            continue
        end

        @inbounds for β in bipartiter.βs, iter = 1:bipartiter.niters
            idxi = rand(1:length(vertices))
            ti = state.config[idxi]
            state.group_sizes[ti] <= 1 && continue

            sc_ti, sc_tinew = space_complexity_singlestep_update(state, adj, degrees_all, log2_sizes, vertices, idxi)
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
            accept && update_state!(state, adj, vertices, idxi, sc_ti, sc_tinew, newloss)
        end
        tc, sc1, sc2 = timespace_complexity_singlestep(state, adj, vertices, log2_sizes)
        @assert state.group_scs ≈ [sc1, sc2]  # sanity check
        if isnull(best) || (maximum(state.group_scs) <= max(bipartiter.sc_target, maximum(best.group_scs)) && (maximum(best.group_scs) >= bipartiter.sc_target || state.loss[] < best.loss[]))
            best = state
        end
    end
    best_tc, = timespace_complexity_singlestep(best, adj, vertices, log2_sizes)
    @debug "best loss = $(round(best.loss[]; digits=3)) space complexities = $(best.group_scs) tc = $(best_tc) groups_sizes = $(best.group_sizes)"
    if maximum(best.group_scs) > bipartiter.sc_target
        @warn "target group size not found, got: $(maximum(best.group_scs)), with time complexity $best_tc."
    end
    return vertices[findall(==(1), best.config)], vertices[findall(==(2), best.config)]
end

function timespace_complexity_singlestep(state, adj, group, log2_sizes)
    g1 = group[findall(==(1), state.config)]
    g2 = group[findall(==(2), state.config)]
    d1 = sum(adj[g1,:], dims=1)
    d2 = sum(adj[g2,:], dims=1)
    dall = sum(adj, dims=1)
    sc1 = sum(i->(d1[i]!=0 && d1[i]!=dall[i] ? Float64(log2_sizes[i]) : 0.0), 1:size(adj,2))
    sc2 = sum(i->(d2[i]!=0 && d2[i]!=dall[i] ? Float64(log2_sizes[i]) : 0.0), 1:size(adj,2))
    tc = sum(i->((d2[i]!=0 || d1[i]!=0) && (d2[i]!=dall[i] && d1[i]!=dall[i]) ? Float64(log2_sizes[i]) : 0.0), 1:size(adj,2))
    return tc, sc1, sc2
end

function space_complexity_singlestep_update(state, adj, degrees_all, log2_sizes, group, idxi)
    begin
        vertex = group[idxi]
        ti = state.config[idxi]
        tinew = 3-ti
        degree_changes = adj[vertex, :]
        δsc_ti = δsc(-, view(state.group_degrees, :, ti), degree_changes, degrees_all, log2_sizes)
        δsc_tinew = δsc(+, view(state.group_degrees, :, tinew), degree_changes, degrees_all, log2_sizes)
        sc_ti = state.group_scs[ti] + δsc_ti
        sc_tinew = state.group_scs[tinew] + δsc_tinew
    end
    return sc_ti, sc_tinew
end

function δsc(f, group_degrees, degree_changes, degrees_all, log2_sizes)
    nnz(degree_changes) == 0 && return 0.0
    return sum(zip(findnz(degree_changes)...)) do (i,v)
        d0 = group_degrees[i]
        D = degrees_all[i]
        d = f(d0, v)
        if d0 == D || d0 == 0     # absent
            if d == D || d == 0   # absent
                0.0
            else
                Float64(log2_sizes[i])
            end
        else                      # not absent
            if d == D || d == 0   # absent
                -Float64(log2_sizes[i])
            else
                0.0
            end
        end
    end
end


@inline function compute_loss(sc1, sc2, gs1, gs2)
    small = min(gs1, gs2)
    max(sc1, sc2) / small * (gs1 + gs2)
end

function update_state!(state, adj, group, idxi, sc_ti, sc_tinew, newloss)
    begin
        ti = state.config[idxi]
        tinew = 3-ti
        dd = adj[group[idxi],:]
        state.group_scs[tinew] = sc_tinew
        state.group_scs[ti] = sc_ti
        state.config[idxi] = tinew
        state.group_sizes[ti] -= 1
        state.group_sizes[tinew] += 1
        state.group_degrees[:,ti] .-= dd
        state.group_degrees[:,tinew] .+= dd
        state.loss[] = newloss
    end
    return state
end


## interfaces
"""
    optimize_sa(code, size_dict; sc_target, βs=0.1:0.2:15.0, niters=1000, ntrails=50, greedy_method=OMEinsum.MinSpaceOut(), greedy_nrepeat=10)

Optimize the einsum code contraction order using the Simulated Annealing + Greedy approach.
This program first recursively cuts the tensors into several groups using simulated annealing,
with maximum group size specifed by `max_group_size` and maximum space complexity specified by `sc_target`,
Then finds the contraction order inside each group with the greedy search algorithm. Other arguments are

* `size_dict`, a dictionary that specifies leg dimensions,
* `sc_target` is the target space complexity, defined as `log2(number of elements in the largest tensor)`,
* `max_group_size` is the maximum size that allowed to used greedy search,
* `βs` is a list of inverse temperature `1/T`,
* `niters` is the number of iteration in each temperature,
* `ntrails` is the number of repetition (with different random seeds),
* `greedy_method` and `greedy_nrepeat` are for configuring the greedy method.

### References
* [Hyper-optimized tensor network contraction](https://arxiv.org/abs/2002.01935)
"""
function optimize_sa(@nospecialize(code::EinCode{ixs,iy}), size_dict; sc_target, βs=0.1:0.2:15.0, niters=1000, ntrails=50, greedy_method=OMEinsum.MinSpaceOut(), greedy_nrepeat=10) where {ixs, iy}
    bipartiter = SABipartite(; sc_target=sc_target, βs=βs, niters=niters, ntrials=ntrails, greedy_method=greedy_method, greedy_nrepeat=greedy_nrepeat)
    recursive_bipartite_optimize(bipartiter, code, size_dict)
end