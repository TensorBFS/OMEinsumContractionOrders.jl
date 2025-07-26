using OMEinsumContractionOrders, Test, Random
using OMEinsumContractionOrders: random_exprtree, ExprTree, ExprInfo,
    ruleset, update_tree!, tcscrw, optimize_subtree!, optimize_tree_sa!, labels, tree_timespace_complexity, fast_log2sumexp2,
    ExprTree, optimize_greedy, _label_dict, Slicer, optimize_tree
using Graphs
using OMEinsum: decorate

@testset "slicer" begin
    log2_sizes = [1, 2,3 ,4.0]
    s = Slicer(log2_sizes, [])
    push!(s, 1)
    @test_throws AssertionError push!(s, 1)
    push!(s, 2)
    push!(s, 3)
    
    replace!(s, 1=>4)
    @test s.log2_sizes == [1, 0.0, 0.0, 0.0]
    @test s.legs == Dict(2=>2.0, 3=>3.0, 4=>4.0)
    @test_throws AssertionError replace!(s, 1=>4)
end

@testset "random expr tree" begin
    function random_regular_eincode(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [[minmax(e.src,e.dst)...] for e in Graphs.edges(g)]
        return OMEinsumContractionOrders.EinCode([ixs..., [[i] for i in Graphs.vertices(g)]...], Int[])
    end
    Random.seed!(2)
    tree = random_exprtree([[1,2,5], [2,3], [2,4]], [5], 5)
    @test tree isa ExprTree
    tree2 = random_exprtree(OMEinsumContractionOrders.EinCode([[1,2,5], [2,3], [2,4]], [5]))
    @test tree isa ExprTree
    code = random_regular_eincode(20, 3)
    optcode = optimize_greedy(code, uniformsize(code, 2); α=0.0, temperature=0.0)
    tree3 = ExprTree(optcode)
    @test tree isa ExprTree
    println(tree)
    labelmap = Dict([v=>k for (k,v) in _label_dict(code)])
    optcode_reconstruct = OMEinsumContractionOrders.NestedEinsum(tree3, labelmap)
    @test optcode == optcode_reconstruct
end

@testset "rules" begin
    LeafNode(id, labels) = ExprTree(ExprInfo(labels, id))
    t1 = ExprTree(LeafNode(3, [1,2]), ExprTree(LeafNode(1,[2,3]), LeafNode(2,[1,4]), ExprInfo([1,2])), ExprInfo([2]))
    t2 = ExprTree(ExprTree(LeafNode(1, [2,3]), LeafNode(2, [1,4]), ExprInfo([1,2,3])), LeafNode(3,[1,2]), ExprInfo([2]))
    t3 = ExprTree(LeafNode(1,[2,3]), LeafNode(2, [1,2]), ExprInfo([2]))
    t4 = ExprTree(ExprTree(LeafNode(1, [2,3]), LeafNode(2, [1,4]), ExprInfo([1,2])), ExprTree(LeafNode(4,[5,1]), LeafNode(3,[1]), ExprInfo([1])), ExprInfo([2]))
    @test ruleset(t1) == 3:4
    @test ruleset(t2) == 1:2
    @test ruleset(t3) == 1:0
    @test ruleset(t4) == 1:4
    log2_sizes = ones(5)
    _tcsc(t, l) = tcscrw(labels(t.left), labels(t.right), labels(t), l, true)
    @test all(_tcsc(t1, log2_sizes) .≈ (2.0, 1.0, log2(10)))
    @test all(_tcsc(t2, log2_sizes) .≈ (2.0, 1.0, log2(14)))
    @test all(_tcsc(t3, log2_sizes) .≈ (1.0, 1.0, log2(10)))
    @test all(_tcsc(t4, log2_sizes) .≈ (2.0, 1.0, log2(8)))
    t11 = update_tree!(copy(t1), 3, [2])
    @test t11 == ExprTree(LeafNode(1,[2,3]), ExprTree(LeafNode(3,[1,2]), LeafNode(2,[1,4]), ExprInfo([2])), ExprInfo([2]))
    t11_ = update_tree!(copy(t1), 4, [1,2])
    @test t11_ == ExprTree(LeafNode(2,[1,4]), ExprTree(LeafNode(1,[2,3]), LeafNode(3,[1,2]), ExprInfo([1,2])), ExprInfo([2]))
    t22 = update_tree!(copy(t2), 1, [1,2])
    @test t22 == ExprTree(ExprTree(LeafNode(1,[2,3]), LeafNode(3,[1,2]), ExprInfo([1,2])), LeafNode(2,[1,4]), ExprInfo([2]))
    t22_ = update_tree!(copy(t2), 2, [2])
    @test t22_ == ExprTree(ExprTree(LeafNode(3, [1,2]), LeafNode(2,[1,4]), ExprInfo([2])), LeafNode(1,[2,3]), ExprInfo([2]))
    t44 = update_tree!(copy(t4), 1, [1,2])
    @test t44 == ExprTree(ExprTree(LeafNode(1,[2,3]), ExprTree(LeafNode(4,[5,1]), LeafNode(3,[1]), ExprInfo([1])), ExprInfo([1,2])), LeafNode(2,[1,4]), ExprInfo([2]))
end

@testset "optimization" begin
    function random_regular_eincode(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [[minmax(e.src,e.dst)...] for e in Graphs.edges(g)]
        return OMEinsumContractionOrders.EinCode([ixs..., [[i] for i in Graphs.vertices(g)]...], Int[])
    end
    Random.seed!(2)
    n = 40
    log2_sizes = rand(n+n÷2) * 2
    code = random_regular_eincode(n, 3)
    optcode = optimize_greedy(code, uniformsize(code, 2); α=0.0, temperature=0.0)
    println(code)
    println(optcode)
    tree = ExprTree(optcode)
    tc0, sc0, rw0 = tree_timespace_complexity(tree, log2_sizes)
    size_dict = Dict([j=>exp2(log2_sizes[j]) for j=1:length(log2_sizes)])
    cc0 = contraction_complexity(OMEinsumContractionOrders.NestedEinsum(tree), size_dict)
    tc0_, sc0_ = cc0.tc, cc0.sc
    @test tc0 ≈ tc0_ && sc0 ≈ sc0_
    opt_tree = copy(tree)
    optimize_subtree!(opt_tree, 100.0, log2_sizes, 5, 2.0, 1.0)
    tc1, sc1, rw0 = tree_timespace_complexity(opt_tree, log2_sizes)
    @test sc1 < sc0 || (sc1 == sc0 && tc1 < tc0)
end

@testset "optimize tree sa" begin
    function random_regular_eincode(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [[minmax(e.src,e.dst)...] for e in Graphs.edges(g)]
        return OMEinsumContractionOrders.EinCode([ixs..., [[i] for i in Graphs.vertices(g)]...], Int[])
    end
    Random.seed!(3)
    n = 60
    ne = n + n÷2
    log2_sizes = ones(ne)
    code = random_regular_eincode(n, 3)
    optcode = optimize_greedy(code, uniformsize(code, 2); α=0.0, temperature=0.0)
    tree = ExprTree(optcode)
    tc0, sc0, rw0 = tree_timespace_complexity(tree, log2_sizes)
    opttree = copy(tree)
    optimize_tree_sa!(opttree, log2_sizes; βs=0.1:0.1:10.0, niters=100, score=ScoreFunction(sc_target=sc0-2.0))
    tc1, sc1, rw1 = tree_timespace_complexity(opttree, log2_sizes)
    @test sc1 < sc0 || (sc1 == sc0 && tc1 < tc0)
end

@testset "sa tree" begin
    function random_regular_eincode(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [[minmax(e.src,e.dst)...] for e in Graphs.edges(g)]
        return OMEinsumContractionOrders.EinCode([ixs..., [[i] for i in Graphs.vertices(g)]...], Int[])
    end
    Random.seed!(2)
    g = random_regular_graph(220, 3)
    code = random_regular_eincode(220, 3)
    res = optimize_greedy(code,uniformsize(code, 2); α=0.0, temperature=0.0)
    cc = contraction_complexity(res, uniformsize(code, 2))
    tc, sc = cc.tc, cc.sc

    @test optimize_tree(res, uniformsize(code, 2); βs=0.1:0.05:20.0, ntrials=0, niters=10, initializer=:greedy, score=ScoreFunction(sc_target=32)) isa OMEinsumContractionOrders.NestedEinsum
    optcode = optimize_tree(res, uniformsize(code, 2); βs=0.1:0.05:20.0, ntrials=2, niters=10, initializer=:greedy, score=ScoreFunction(sc_target=32))
    cc = contraction_complexity(optcode, uniformsize(code, 2))
    @test cc.sc <= 32

    # contraction test
    code = random_regular_eincode(50, 3)
    codek = optimize_greedy(code, uniformsize(code, 2); α=0.0, temperature=0.0)
    codeg = optimize_tree(code, uniformsize(code, 2); initializer=:random, βs=0.1:0.05:20.0, ntrials=2, niters=10, score=ScoreFunction(sc_target=12))
    cc = contraction_complexity(codek, uniformsize(code, 2))
    @test cc.sc <= 12
    xs = [[2*randn(2, 2) for i=1:75]..., [randn(2) for i=1:50]...]
    resg = decorate(codeg)(xs...)
    resk = decorate(codek)(xs...)
    @test resg ≈ resk

    # contraction test
    code = random_regular_eincode(50, 3)
    codek = optimize_greedy(code, uniformsize(code, 2); α=0.0, temperature=0.0)
    codeg = optimize_tree(codek, uniformsize(code, 2); initializer=:specified, βs=0.1:0.05:20.0, ntrials=2, niters=10, score=ScoreFunction(sc_target=12))
    cc = contraction_complexity(codek, uniformsize(code, 2))
    @test cc.sc <= 12
    xs = [[2*randn(2, 2) for i=1:75]..., [randn(2) for i=1:50]...]
    resg = decorate(codeg)(xs...)
    resk = decorate(codek)(xs...)
    @test resg ≈ resk
end

@testset "fast log2sumexp2" begin
    a, b, c = randn(3)
    @test fast_log2sumexp2(a, b) ≈ log2(sum(exp2.([a,b])))
    @test fast_log2sumexp2(a, b, c) ≈ log2(sum(exp2.([a,b,c])))
end

@testset "path decomposition" begin
    function random_regular_eincode(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [[minmax(e.src,e.dst)...] for e in Graphs.edges(g)]
        return OMEinsumContractionOrders.EinCode([ixs..., [[i] for i in Graphs.vertices(g)]...], Int[])
    end
    Random.seed!(2)
    code = random_regular_eincode(50, 3)
    codek = optimize_greedy(code, uniformsize(code, 2); α=0.0, temperature=0.0)
    codeg = optimize_tree(code, uniformsize(code, 2); initializer=:random, βs=0.1:0.05:20.0, ntrials=2, niters=10, score=ScoreFunction(sc_target=12), decomposition_type=PathDecomp())
    @test !OMEinsumContractionOrders.is_path_decomposition(codek)
    @test OMEinsumContractionOrders.is_path_decomposition(codeg)
    @test OMEinsumContractionOrders.is_path_decomposition(OMEinsumContractionOrders.SlicedEinsum([1,2], codeg))
    cc = contraction_complexity(codek, uniformsize(code, 2))
    @test cc.sc <= 12
    cc = contraction_complexity(codeg, uniformsize(code, 2))
    @test cc.sc <= 12
    xs = [[2*randn(2, 2) for i=1:75]..., [randn(2) for i=1:50]...]
    resg = decorate(codeg)(xs...)
    resk = decorate(codek)(xs...)
    @test resg ≈ resk

    # no optimization
    codeg = optimize_tree(code, uniformsize(code, 2); initializer=:random, βs=0.1:0.05:-20.0, ntrials=2, niters=10, score=ScoreFunction(sc_target=12), decomposition_type=PathDecomp())
    cc = contraction_complexity(codeg, uniformsize(code, 2))
    @test OMEinsumContractionOrders.is_path_decomposition(codeg)
    @test cc.sc > 12
end