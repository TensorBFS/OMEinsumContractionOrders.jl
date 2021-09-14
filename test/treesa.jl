using OMEinsumContractionOrders, Test, Random
using OMEinsumContractionOrders: random_exprtree, ExprTree, ExprInfo, ruleset, update_tree!, tcsc, optimize_subtree!, LeafNode, optimize_tree_sa, labels
using OMEinsum, LightGraphs

@testset "random expr tree" begin
    function random_regular_eincode(n, k)
        g = LightGraphs.random_regular_graph(n, k)
        ixs = [minmax(e.src,e.dst) for e in LightGraphs.edges(g)]
        return EinCode((ixs..., [(i,) for i in LightGraphs.vertices(g)]...), ())
    end
    Random.seed!(2)
    tree = random_exprtree([[1,2,5], [2,3], [2,4]], [5], 5)
    @test tree isa OMEinsumContractionOrders.ExprTree
    tree2 = random_exprtree(EinCode(((1,2,5), (2,3), (2,4)), (5,)))
    @test tree isa OMEinsumContractionOrders.ExprTree
    code = random_regular_eincode(20, 3)
    optcode = optimize_greedy(code, uniformsize(code, 2))
    tree3 = ExprTree(optcode)
    @test tree isa OMEinsumContractionOrders.ExprTree
    labelmap = Dict([v=>k for (k,v) in OMEinsumContractionOrders._label_dict(code)])
    optcode_reconstruct = NestedEinsum(tree3, labelmap)
    @test optcode == optcode_reconstruct
end

@testset "rules" begin
    t1 = ExprTree(LeafNode(3, [1,2]), ExprTree(LeafNode(1,[2,3]), LeafNode(2,[1,4]), ExprInfo([1,2])), ExprInfo([2]))
    t2 = ExprTree(ExprTree(LeafNode(1, [2,3]), LeafNode(2, [1,4]), ExprInfo([1,2,3])), LeafNode(3,[1,2]), ExprInfo([2]))
    t3 = ExprTree(LeafNode(1,[2,3]), LeafNode(2, [1,2]), ExprInfo([2]))
    t4 = ExprTree(ExprTree(LeafNode(1, [2,3]), LeafNode(2, [1,4]), ExprInfo([1,2])), ExprTree(LeafNode(4,[5,1]), LeafNode(3,[1]), ExprInfo([1])), ExprInfo([2]))
    @test ruleset(t1) == 3:4
    @test ruleset(t2) == 1:2
    @test ruleset(t3) == 1:0
    @test ruleset(t4) == 1:4
    log2_sizes = ones(5)
    _tcsc(t, l) = tcsc(labels(t.left), labels(t.right), labels(t), l)
    @test _tcsc(t1, log2_sizes) == (2.0, 1.0)
    @test _tcsc(t2, log2_sizes) == (2.0, 1.0)
    @test _tcsc(t3, log2_sizes) == (1.0, 1.0)
    @test _tcsc(t4, log2_sizes) == (2.0, 1.0)
    t11 = update_tree!(copy(t1), 3)
    @test t11 == ExprTree(LeafNode(1,[2,3]), ExprTree(LeafNode(3,[1,2]), LeafNode(2,[1,4]), ExprInfo([2])), ExprInfo([2]))
    t11_ = update_tree!(copy(t1), 4)
    @test t11_ == ExprTree(LeafNode(2,[1,4]), ExprTree(LeafNode(1,[2,3]), LeafNode(3,[1,2]), ExprInfo([1,2])), ExprInfo([2]))
    t22 = update_tree!(copy(t2), 1)
    @test t22 == ExprTree(ExprTree(LeafNode(1,[2,3]), LeafNode(3,[1,2]), ExprInfo([1,2])), LeafNode(2,[1,4]), ExprInfo([2]))
    t22_ = update_tree!(copy(t2), 2)
    @test t22_ == ExprTree(ExprTree(LeafNode(3, [1,2]), LeafNode(2,[1,4]), ExprInfo([2])), LeafNode(1,[2,3]), ExprInfo([2]))
    t44 = update_tree!(copy(t4), 1)
    @test t44 == ExprTree(ExprTree(LeafNode(1,[2,3]), ExprTree(LeafNode(4,[5,1]), LeafNode(3,[1]), ExprInfo([1])), ExprInfo([1,2])), LeafNode(2,[1,4]), ExprInfo([2]))
end

@testset "optimization" begin
    function random_regular_eincode(n, k)
        g = LightGraphs.random_regular_graph(n, k)
        ixs = [minmax(e.src,e.dst) for e in LightGraphs.edges(g)]
        return EinCode((ixs..., [(i,) for i in LightGraphs.vertices(g)]...), ())
    end
    Random.seed!(2)
    n = 20
    log2_sizes = rand(n+n÷2) * 2
    code = random_regular_eincode(20, 3)
    optcode = optimize_greedy(code, uniformsize(code, 2))
    tree = ExprTree(optcode)
    tc0, sc0 = OMEinsum.timespace_complexity(tree, exp2.(log2_sizes))
    size_dict = Dict([j=>exp2(log2_sizes[j]) for j=1:length(log2_sizes)])
    tc0_, sc0_ = OMEinsum.timespace_complexity(NestedEinsum(tree), size_dict)
    @test tc0 ≈ tc0_ && sc0 ≈ sc0_
    opt_tree = optimize_subtree!(copy(tree), 100.0, log2_sizes, 5, 2.0)
    tc1, sc1 = OMEinsum.timespace_complexity(opt_tree, exp2.(log2_sizes))
    @test sc1 < sc0
end

@testset "optimize tree sa" begin
    function random_regular_eincode(n, k)
        g = LightGraphs.random_regular_graph(n, k)
        ixs = [minmax(e.src,e.dst) for e in LightGraphs.edges(g)]
        return EinCode((ixs..., [(i,) for i in LightGraphs.vertices(g)]...), ())
    end
    Random.seed!(3)
    n = 60
    ne = n + n÷2
    log2_sizes = ones(ne)
    code = random_regular_eincode(n, 3)
    optcode = optimize_greedy(code, uniformsize(code, 2))
    tree = ExprTree(optcode)
    tc0, sc0 = OMEinsum.timespace_complexity(tree, exp2.(log2_sizes))
    opttree = optimize_tree_sa(tree, log2_sizes; sc_target=sc0-2.0, βs=0.1:0.1:10, ntrials=2, niters=100, sc_weight=3.0)
    tc1, sc1 = OMEinsum.timespace_complexity(opttree, exp2.(log2_sizes))
    @test sc1 < sc0 || (sc1 == sc0 && tc1 < tc0)
end

@testset "sa tree" begin
    function random_regular_eincode(n, k)
        g = LightGraphs.random_regular_graph(n, k)
        ixs = [minmax(e.src,e.dst) for e in LightGraphs.edges(g)]
        return EinCode((ixs..., [(i,) for i in LightGraphs.vertices(g)]...), ())
    end
    Random.seed!(2)
    g = random_regular_graph(220, 3)
    code = random_regular_eincode(220, 3)
    res = optimize_greedy(code,uniformsize(code, 2))
    tc, sc = OMEinsum.timespace_complexity(res, uniformsize(code, 2))

    optcode = optimize_tree(res,uniformsize(code, 2); sc_target=32, βs=0.1:0.1:10, ntrials=2, niters=100, sc_weight=3.0)
    tc, sc = OMEinsum.timespace_complexity(optcode, uniformsize(code, 2))
    @test sc <= 32

    # contraction test
    code = random_regular_eincode(50, 3)
    codek = optimize_greedy(code, uniformsize(code, 2))
    codeg = optimize_tree(codek, uniformsize(code, 2))
    tc, sc = OMEinsum.timespace_complexity(codek, uniformsize(code, 2))
    @test sc <= 12
    xs = [[2*randn(2, 2) for i=1:75]..., [randn(2) for i=1:50]...]
    resg = codeg(xs...)
    resk = codek(xs...)
    @test resg ≈ resk
end

