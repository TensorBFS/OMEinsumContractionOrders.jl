using OMEinsumContractionOrders, Test, Random
using OMEinsumContractionOrders: random_exprtree, ExprTree, ExprInfo, ruleset, update_tree!, tcscrw, optimize_subtree!, optimize_tree_sa!, labels, tree_timespace_complexity, fast_log2sumexp2
using OMEinsum, Graphs

@testset "random expr tree" begin
    function random_regular_eincode(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [minmax(e.src,e.dst) for e in Graphs.edges(g)]
        return EinCode((ixs..., [(i,) for i in Graphs.vertices(g)]...), ())
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
        ixs = [minmax(e.src,e.dst) for e in Graphs.edges(g)]
        return EinCode((ixs..., [(i,) for i in Graphs.vertices(g)]...), ())
    end
    Random.seed!(2)
    n = 40
    log2_sizes = rand(n+n÷2) * 2
    code = random_regular_eincode(n, 3)
    optcode = optimize_greedy(code, uniformsize(code, 2))
    tree = ExprTree(optcode)
    tc0, sc0, rw0 = tree_timespace_complexity(tree, log2_sizes)
    size_dict = Dict([j=>exp2(log2_sizes[j]) for j=1:length(log2_sizes)])
    tc0_, sc0_ = OMEinsum.timespace_complexity(NestedEinsum(tree), size_dict)
    @test tc0 ≈ tc0_ && sc0 ≈ sc0_
    opt_tree = copy(tree)
    optimize_subtree!(opt_tree, 100.0, log2_sizes, 5, 2.0, 1.0)
    tc1, sc1, rw0 = tree_timespace_complexity(opt_tree, log2_sizes)
    @test sc1 < sc0 || (sc1 == sc0 && tc1 < tc0)
end

@testset "optimize tree sa" begin
    function random_regular_eincode(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [minmax(e.src,e.dst) for e in Graphs.edges(g)]
        return EinCode((ixs..., [(i,) for i in Graphs.vertices(g)]...), ())
    end
    Random.seed!(3)
    n = 60
    ne = n + n÷2
    log2_sizes = ones(ne)
    code = random_regular_eincode(n, 3)
    optcode = optimize_greedy(code, uniformsize(code, 2))
    tree = ExprTree(optcode)
    tc0, sc0, rw0 = tree_timespace_complexity(tree, log2_sizes)
    opttree = copy(tree)
    optimize_tree_sa!(opttree, log2_sizes, Slicer(log2_sizes, 0, []); sc_target=sc0-2.0, βs=0.1:0.1:10.0, niters=100, sc_weight=1.0, rw_weight=1.0)
    tc1, sc1, rw1 = tree_timespace_complexity(opttree, log2_sizes)
    @test sc1 < sc0 || (sc1 == sc0 && tc1 < tc0)

    slicer = Slicer(log2_sizes, 5, [])
    optimize_tree_sa!(opttree, log2_sizes, slicer; sc_target=sc0-2.0, βs=0.1:0.1:10.0, niters=100, sc_weight=1.0, rw_weight=1.0)
    tc2, sc2, rw2 = tree_timespace_complexity(opttree, slicer.log2_sizes)
    @test tc2 <= tc1 + 3
    @test sc2 <= sc1 + 3
    @test length(slicer) == 5
    @test all(l->(slicer.log2_sizes[l]==1 && !haskey(slicer.legs, l)) || (slicer.log2_sizes[l]==0 && haskey(slicer.legs, l)), 1:length(log2_sizes))
    @test sc1 < sc0 || (sc1 == sc0 && tc1 < tc0)
end

@testset "sa tree" begin
    function random_regular_eincode(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [minmax(e.src,e.dst) for e in Graphs.edges(g)]
        return EinCode((ixs..., [(i,) for i in Graphs.vertices(g)]...), ())
    end
    Random.seed!(2)
    g = random_regular_graph(220, 3)
    code = random_regular_eincode(220, 3)
    res = optimize_greedy(code,uniformsize(code, 2))
    tc, sc = OMEinsum.timespace_complexity(res, uniformsize(code, 2))

    @test optimize_tree(res,uniformsize(code, 2); sc_target=32, βs=0.1:0.05:20.0, ntrials=0, niters=10, sc_weight=1.0, rw_weight=1.0) isa SlicedEinsum
    optcode = optimize_tree(res,uniformsize(code, 2); sc_target=32, βs=0.1:0.05:20.0, ntrials=2, niters=10, sc_weight=1.0, rw_weight=1.0)
    tc, sc = OMEinsum.timespace_complexity(optcode, uniformsize(code, 2))
    @test sc <= 32
    @test length(optcode.slicing) == 0

    # contraction test
    code = random_regular_eincode(50, 3)
    codek = optimize_greedy(code, uniformsize(code, 2))
    codeg = optimize_tree(code, uniformsize(code, 2); initializer=:random)
    tc, sc = OMEinsum.timespace_complexity(codek, uniformsize(code, 2))
    @test sc <= 12
    xs = [[2*randn(2, 2) for i=1:75]..., [randn(2) for i=1:50]...]
    resg = codeg(xs...)
    resk = codek(xs...)
    @test resg ≈ resk

    # contraction test
    code = random_regular_eincode(50, 3)
    codek = optimize_greedy(code, uniformsize(code, 2))
    codeg = optimize_tree(codek, uniformsize(code, 2); initializer=:specified)
    tc, sc = OMEinsum.timespace_complexity(codek, uniformsize(code, 2))
    @test sc <= 12
    xs = [[2*randn(2, 2) for i=1:75]..., [randn(2) for i=1:50]...]
    resg = codeg(xs...)
    resk = codek(xs...)
    @test resg ≈ resk
end

@testset "fast log2sumexp2" begin
    a, b, c = randn(3)
    @test fast_log2sumexp2(a, b) ≈ log2(sum(exp2.([a,b])))
    @test fast_log2sumexp2(a, b, c) ≈ log2(sum(exp2.([a,b,c])))
end

@testset "slicing" begin
    function random_regular_eincode(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [minmax(e.src,e.dst) for e in Graphs.edges(g)]
        return EinCode((ixs..., [(i,) for i in Graphs.vertices(g)]...), ())
    end
    # contraction test
    Random.seed!(2)
    code = random_regular_eincode(100, 3)
    code0 = optimize_greedy(code, uniformsize(code, 2))
    codek = optimize_tree(code0, uniformsize(code, 2); initializer=:specified, nslices=0, niters=5)
    codeg = optimize_tree(code0, uniformsize(code, 2); initializer=:specified, nslices=5, niters=5)
    tc0, sc0 = OMEinsum.timespace_complexity(codek, uniformsize(code, 2))
    tc, sc = OMEinsum.timespace_complexity(codeg, uniformsize(code, 2))
    @show tc, sc, tc0, sc0
    @test sc <= sc0 - 4
    xs = [[2*randn(2, 2) for i=1:150]..., [randn(2) for i=1:100]...]
    resg = codeg(xs...)
    resk = codek(xs...)
    @test resg ≈ resk

    # with open edges
    Random.seed!(2)
    code = EinCode(random_regular_eincode(100, 3).ixs, [3,81,2])
    codek = optimize_tree(code0, uniformsize(codek, 2); initializer=:specified, nslices=0, niters=5)
    codeg = optimize_tree(code0, uniformsize(codeg, 2); initializer=:specified, nslices=5, niters=5)
    tc0, sc0 = OMEinsum.timespace_complexity(codek, uniformsize(code, 2))
    tc, sc = OMEinsum.timespace_complexity(codeg, uniformsize(code, 2))
    fl = OMEinsum.flop(codeg, uniformsize(code, 2))
    @test tc ≈ log2(fl)
    @show tc, sc, tc0, sc0
    @test sc <= sc0 - 4
    @test sc <= 17
    xs = [[2*randn(2, 2) for i=1:150]..., [randn(2) for i=1:100]...]
    resg = codeg(xs...)
    resk = codek(xs...)
    @test resg ≈ resk

    # slice with fixed slices
    Random.seed!(2)
    code = EinCode(random_regular_eincode(20, 3).ixs, [3,10,2])
    code0 = optimize_tree(code, uniformsize(code, 2); nslices=5, fixed_slices=[5,3,8,1,2,4,11])
    code1 = optimize_tree(code, uniformsize(code, 2); ntrials=1, nslices=5)
    code2 = optimize_tree(code, uniformsize(code, 2); ntrials=1, nslices=5, fixed_slices=[5,3])
    code3 = optimize_tree(code, uniformsize(code, 2); ntrials=1, nslices=5, fixed_slices=[5,1,4,3,2])
    code4 = optimize_tree(code, uniformsize(code, 2); ntrials=1, fixed_slices=[5,1,4,3,2])
    xs = [[2*randn(2, 2) for i=1:30]..., [randn(2) for i=1:20]...]
    @test length(code0.slicing) == 7 && code0.slicing.legs == [5,3,8,1,2,4,11]
    @test length(code2.slicing) == 5 && code2.slicing.legs[1:2] == [5,3]
    @test length(code3.slicing) == 5 && code3.slicing.legs == [5,1,4,3,2]
    @test length(code4.slicing) == 5 && code4.slicing.legs == [5,1,4,3,2]
    @test code1(xs...) ≈ code2(xs...)
    @test code1(xs...) ≈ code3(xs...)
end