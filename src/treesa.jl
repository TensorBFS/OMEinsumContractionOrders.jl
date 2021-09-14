using OMEinsum.ContractionOrder: ContractionTree, log2sumexp2

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

space_complexity(t::ExprTree, log2_sizes) = length(labels(t))==0 ? 0.0 : log2sumexp2(getindex.(Ref(log2_sizes), labels(t)))
tcsc(t::ExprTree, log2_sizes) = OMEinsum._timespace_complexity([labels(t.left), labels(t.right)], labels(t), log2_sizes)
siblings(t::ExprTree) = Any[t.left, t.right]
siblings(::LeafNode) = Any[]
Base.copy(t::ExprTree) = ExprTree(copy(t.left), copy(t.right), ExprInfo(copy(t.info.out_dims)))
Base.copy(t::LeafNode) = LeafNode(t.tensorid, copy(t.labels))
labels(t::ExprTree) = t.info.out_dims
labels(t::LeafNode) = t.labels
maxlabel(t::ExprTree) = max(isempty(labels(t)) ? 0 : maximum(labels(t)), maxlabel(t.left), maxlabel(t.right))
maxlabel(t::LeafNode) = maximum(isempty(labels(t)) ? 0 : labels(t))
Base.:(==)(t1::ExprTree, t2::ExprTree) = _equal(t1, t2)
_equal(t1::ExprTree, t2::ExprTree) = _equal(t1.left, t2.left) && _equal(t1.right, t2.right) && t1.info == t2.info
_equal(t1::LeafNode, t2::LeafNode) = t1.tensorid == t2.tensorid
_equal(t1::Vector, t2::Vector) = Set(t1) == Set(t2)
_equal(a, b) = false
Base.:(==)(t1::ExprInfo, t2::ExprInfo) = _equal(t1.out_dims, t2.out_dims)

function optimize_tree_sa(tree::ExprTree, log2_sizes; βs, niters, ntrials, sc_target, sc_weight)
    best_tree = tree
    best_tc, best_sc = tree_timespace_complexity(tree, log2_sizes)
    for _ = 1:ntrials
        ctree = copy(tree)
        @inbounds for β in βs, _ = 1:niters
            optimize_subtree!(ctree, β, log2_sizes, sc_target, sc_weight)  # single sweep
        end
        tc, sc = tree_timespace_complexity(ctree, log2_sizes)
        if sc < best_sc || (sc <= best_sc && tc < best_tc)
            best_tree, best_tc, best_sc = ctree, tc, sc
        end
    end
    @debug "best space complexities = $best_tc time complexity = $best_sc"
    @show best_sc
    if best_sc > sc_target
        @warn "target space complexity not found, got: $best_sc, with time complexity $best_tc."
    end
    return best_tree
end

function OMEinsum.timespace_complexity(tree::ExprTree, size_vec)
    tree_timespace_complexity(tree, log2.(size_vec))
end
function tree_timespace_complexity(tree::LeafNode, log2_sizes)
    -Inf, log2sumexp2(getindex.(Ref(log2_sizes), tree.labels))
end
function tree_timespace_complexity(tree::ExprTree, log2_sizes)
    tcl, scl = tree_timespace_complexity(tree.left, log2_sizes)
    tcr, scr = tree_timespace_complexity(tree.right, log2_sizes)
    tc, sc = tcsc(tree, log2_sizes)
    return log2sumexp2([tc, tcl, tcr]), max(sc, scl, scr)
end

function random_exprtree(@nospecialize(code::EinCode{ixs, iy})) where {ixs, iy}
    labels = _label_dict(code)
    return random_exprtree([[labels[l] for l in ix] for ix in ixs], [labels[l] for l in iy], length(labels))
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
    _random_exprtree(ixs, outercount, allcount, Ref(0))
end
function _random_exprtree(ixs::Vector{Vector{Int}}, outercount::Vector{Int}, allcount::Vector{Int}, k)
    n = length(ixs)
    if n == 1
        k[] += 1
        return LeafNode(k[], ixs[1])
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
    return ExprTree(_random_exprtree(ixs[mask], outercount1, allcount, k), _random_exprtree(ixs[(!).(mask)], outercount2, allcount, k), info)
end

function optimize_subtree!(tree, β, log2_sizes, sc_target, sc_weight)
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
            optimize_subtree!(subtree, β, log2_sizes, sc_target, sc_weight)
        end
    end
    return tree
end

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
