using OMEinsum.ContractionOrder: ContractionTree

struct ExprInfo
    out_dims::Vector{Int}
end

mutable struct ExprTree
    left::Union{ExprTree,Vector{Int}}
    right::Union{ExprTree,Vector{Int}}
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
print_expr(io::IO, expr::Vector, level=0) = print(io, " "^(2*level), expr)
Base.show(io::IO, expr::ExprTree) = print_expr(io, expr, 0)
Base.show(io::IO, ::MIME"text/plain", expr::ExprTree) = show(io, expr)

space_complexity(t::ExprTree, log2_sizes) = OMEinsum.log2sumexp2(getindex.(Ref(log2_sizes), t.out_dims))
tcsc(t::ExprTree, log2_sizes) = OMEinsum._timespace_complexity([labels(t.left), labels(t.right)], labels(t), log2_sizes)
siblings(t::ExprTree) = Any[t.left, t.right]
siblings(::Vector{Int}) = Any[]
Base.copy(t::ExprTree) = ExprTree(copy(t.left), copy(t.right), ExprInfo(copy(t.info.out_dims)))
labels(t::ExprTree) = t.info.out_dims
labels(t::Vector{Int}) = t
Base.:(==)(t1::ExprTree, t2::ExprTree) = _equal(t1, t2)
_equal(t1::ExprTree, t2::ExprTree) = _equal(t1.left, t2.left) && _equal(t1.right, t2.right) && t1.info == t2.info
_equal(t1::Vector, t2::Vector) = Set(t1) == Set(t2)
_equal(a, b) = false
Base.:(==)(t1::ExprInfo, t2::ExprInfo) = _equal(t1.out_dims, t2.out_dims)

function optimize_tree_sa(tree::ExprTree, log2_sizes; βs, niters, ntrials, sc_target)
    best_tree = tree
    best_tc, best_sc = timespace_complexity(tree, log2_sizes)
    for _ = 1:ntrials
        ctree = copy(tree)
        if state.group_sizes[1]==0 || state.group_sizes[2] == 0
            continue
        end
        @inbounds for β in βs, _ = 1:niters
            _optimize_subtree!(ctree, β, log2_sizes)  # single sweep
        end
        tc, sc = timespace_complexity(ctree, log2_sizes)
        if sc < best_sc || (sc <= best_sc && tc < best_tc)
            best_tree, best_tc, best_sc = ctree, tc, sc
        end
    end
    @debug "best space complexities = $best_tc time complexity = $best_sc"
    if sc > sc_target
        @warn "target space complexity not found, got: $best_sc, with time complexity $best_tc."
    end
    return tree
end

function random_exprtree(@nospecialize(code::EinCode{ixs, iy})) where {ixs, iy}
    labels = _label_dict(code)
    return random_exprtree([[labels[l] for l in ix] for ix in ixs], [labels[l] for l in iy], length(labels))
end
function _label_dict(@nospecialize(code::EinCode{ixs, iy})) where {ixs, iy}
    ixsv, iyv = collect.(ixs), collect(iy)
    v = unique(vcat(ixsv..., iyv))
    labels = Dict(zip(v, 1:length(v)))
    return labels
end

exprtree(code::NestedEinsum) = _exprtree(code, _label_dict(OMEinsum.flatten(code)))
function _exprtree(code::NestedEinsum, labels)
    @assert length(code.args) == 2
    ExprTree(map(enumerate(code.args)) do (i,arg)
        if arg isa Int
            collect([labels[i] for i=OMEinsum.getixs(code.eins)[i]])
        else
            res = _exprtree(arg, labels)
        end
    end..., ExprInfo([labels[i] for i=OMEinsum.getiy(code.eins)]))
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
    _random_exprtree(ixs, outercount, allcount)
end
function _random_exprtree(ixs::Vector{Vector{Int}}, outercount::Vector{Int}, allcount::Vector{Int})
    n = length(ixs)
    if n == 1
        return ixs[1]
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
    return ExprTree(_random_exprtree(ixs[mask], outercount1, allcount), _random_exprtree(ixs[(!).(mask)], outercount2, allcount), info)
end

function _optimize_subtree!(tree, β, log2_sizes)
    rst = ruleset(tree)
    if !isempty(rst)
        rule = rand(rst)
        sc = space_complexity(tree, log2_sizes)
        dtc, dsc = tcsc_diff(tree, rule, log2_sizes)
        #log2(α*RW + tc) is the original `tc` term, which also optimizes read-write overheads.
        dE = (max(sc, sc+dsc) > sc_target ? sc_weight : 0) * dsc + dtc
        if rand() < exp(-β*dE)
            update_tree!(tree, rule)
        end
        for subtree in siblings(tree)
            _optimize_subtree!(subtree, β, log2_sizes)
        end
    end
    return tree
end

ruleset(::Vector{Int}) = 1:-1
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
        tc0, sc0 = tcsc(tree.left, log2_sizes) .+ tcsc(tree, log2_sizes)
        labels_ac = (labels(tree.left.left) ∪ labels(tree.right)) ∩ (labels(tree.left.right) ∪ labels(tree))
        newtree = build_tree(build_tree(tree.left.left, tree.right, labels_ac), tree.left.right, labels(tree))
        tc1, sc1 = tcsc(newtree.left, log2_sizes) .+ tcsc(newtree, log2_sizes)
        return tc1-tc0, sc1-sc0
    elseif rule == 2 # (a,b), c -> (c,b),a
        tc0, sc0 = tcsc(tree.left, log2_sizes) .+ tcsc(tree, log2_sizes)
        labels_cb = (labels(tree.left.right) ∪ labels(tree.right)) ∩ (labels(tree.left.left) ∪ labels(tree))
        newtree = build_tree(build_tree(tree.right, tree.left.right, labels_cb), tree.left.left, labels(tree))
        tc1, sc1 = tcsc(newtree.left, log2_sizes) .+ tcsc(newtree, log2_sizes)
        return tc1-tc0, sc1-sc0
    elseif rule == 3 # a,(b,c) -> b,(a,c)
        tc0, sc0 = tcsc(tree.right, log2_sizes) .+ tcsc(tree, log2_sizes)
        labels_ac = (labels(tree.left) ∪ labels(tree.right.right)) ∩ (labels(tree.right.left) ∪ labels(tree))
        newtree = build_tree(tree.right.left, build_tree(tree.left, tree.right.right, labels_ac), labels(tree))
        tc1, sc1 = tcsc(newtree.right, log2_sizes) .+ tcsc(newtree, log2_sizes)
        return tc1-tc0, sc1-sc0
    else  # a,(b,c) -> c,(b,a)
        tc0, sc0 = tcsc(tree.right, log2_sizes) .+ tcsc(tree, log2_sizes)
        labels_ba = (labels(tree.left) ∪ labels(tree.right.left)) ∩ (labels(tree.right.right) ∪ labels(tree))
        newtree = build_tree(tree.right.right, build_tree(tree.right.left, tree.left, labels_ba), labels(tree))
        tc1, sc1 = tcsc(newtree.right, log2_sizes) .+ tcsc(newtree, log2_sizes)
        return tc1-tc0, sc1-sc0
    end
end

function build_tree(left, right, out_dims)
    ExprTree(left, right, ExprInfo(out_dims))
end

function update_tree!(tree::ExprTree, rule::Int)
    if rule == 1 # (a,b), c -> (a,c),b
        b, c = tree.left.right, tree.right
        tree.left.right = c
        tree.right = b
        tree.left.info = ExprInfo((labels(tree.left.left) ∪ labels(tree.left.right)) ∩ (labels(tree) ∪ labels(tree.right)))
    elseif rule == 2 # (a,b), c -> (c,b),a
        a, c = tree.left.left, tree.right
        tree.left.left = c
        tree.right = a
        tree.left.info = ExprInfo((labels(tree.left.left) ∪ labels(tree.left.right)) ∩ (labels(tree) ∪ labels(tree.right)))
    elseif rule == 3 # a,(b,c) -> b,(a,c)
        a, b = tree.left, tree.right.left
        tree.left = b
        tree.right.left = a
        tree.right.info = ExprInfo((labels(tree.right.left) ∪ labels(tree.right.right)) ∩ (labels(tree) ∪ labels(tree.left)))
    else  # a,(b,c) -> c,(b,a)
        a, c = tree.left, tree.right.right
        tree.left = c
        tree.right.right = a
        tree.right.info = ExprInfo((labels(tree.right.left) ∪ labels(tree.right.right)) ∩ (labels(tree) ∪ labels(tree.left)))
    end
    return tree
end
