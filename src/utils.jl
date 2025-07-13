function log2sumexp2(s)
    ms = maximum(s)
    return log2(sum(x->exp2(x - ms), s)) + ms
end

# three different cases:
# sc_target = 0: get the tree with the smallest space complexity. (reducing space complexity)
# sc_target is a normal number: get the tree with `sc <= sc_target` first, then consider the time-space-readwrite complexity. (optimization with a target)
# sc_target = Inf: get the tree with the smallest time-space-readwrite complexity. (slicing)
function best_tree(trees, tcs, scs, rws, ntrials, rw_weight, sc_target)
    best_tree, best_tc, best_sc, best_rw = first(trees), first(tcs), first(scs), first(rws)
    
    for i in 2:ntrials
        if (best_sc > sc_target) && (scs[i] < best_sc)
            best_tree, best_tc, best_sc, best_rw = trees[i], tcs[i], scs[i], rws[i]
        elseif exp2(tcs[i]) + rw_weight * exp2(rws[i]) < exp2(best_tc) + rw_weight * exp2(rws[i])
            best_tree, best_tc, best_sc, best_rw = trees[i], tcs[i], scs[i], rws[i]
        end
    end

    return best_tree, best_tc, best_sc, best_rw
end