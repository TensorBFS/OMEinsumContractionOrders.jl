function log2sumexp2(s)
    ms = maximum(s)
    return log2(sum(x->exp2(x - ms), s)) + ms
end

# three different cases:
# sc_target = 0: get the tree with the smallest space complexity. (reducing space complexity)
# sc_target is a normal number: get the tree with `sc <= sc_target` first, then consider the time-space-readwrite complexity. (optimization with a target)
# sc_target = Inf: get the tree with the smallest time-space-readwrite complexity. (slicing)
function find_best_tree(tcs, scs, rws, ntrials, rw_weight, sc_target)
    best_id = 1
    for i in 2:ntrials
        if (scs[best_id] > sc_target) && (scs[i] < scs[best_id])
            best_id = i
        elseif exp2(tcs[i]) + rw_weight * exp2(rws[i]) < exp2(tcs[best_id]) + rw_weight * exp2(rws[best_id])
            best_id = i
        end
    end

    return best_id
end