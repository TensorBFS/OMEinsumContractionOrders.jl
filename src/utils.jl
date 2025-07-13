function log2sumexp2(s)
    ms = maximum(s)
    return log2(sum(x->exp2(x - ms), s)) + ms
end

function best_tree(trees, tcs, scs, rws, ntrials, rw_weight, sc_priority::Bool)
    best_tree, best_tc, best_sc, best_rw = first(trees), first(tcs), first(scs), first(rws)
    for i=2:ntrials
        if (sc_priority && (scs[i] < best_sc)) || (!sc_priority && exp2(tcs[i]) + rw_weight * exp2(rws[i]) < exp2(best_tc) + rw_weight * exp2(rws[i]))
            best_tree, best_tc, best_sc, best_rw = trees[i], tcs[i], scs[i], rws[i]
        end
    end
    return best_tree, best_tc, best_sc, best_rw
end