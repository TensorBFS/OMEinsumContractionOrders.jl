using OMEinsumContractionOrders, Test, Random
using OMEinsumContractionOrders: random_exprtree, ExprTree, ExprInfo, ruleset, update_tree!, tcsc, optimize_subtree!, LeafNode
using OMEinsum, LightGraphs

function random_regular_eincode(n, k)
    g = LightGraphs.random_regular_graph(n, k)
    ixs = [minmax(e.src,e.dst) for e in LightGraphs.edges(g)]
    return EinCode((ixs..., [(i,) for i in LightGraphs.vertices(g)]...), ())
end

@testset "random expr tree" begin
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
    @test tcsc(t1, log2_sizes) == (2.0, 1.0)
    @test tcsc(t2, log2_sizes) == (2.0, 1.0)
    @test tcsc(t3, log2_sizes) == (1.0, 1.0)
    @test tcsc(t4, log2_sizes) == (2.0, 1.0)
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
    Random.seed!(2)
    n = 20
    log2_sizes = rand(n+n÷2) * 2
    code = random_regular_eincode(20, 3)
    optcode = optimize_greedy(code, uniformsize(code, 2))
    tree = ExprTree(optcode)
    tc0, sc0 = timespace_complexity(tree, exp2.(log2_sizes))
    size_dict = Dict([j=>exp2(log2_sizes[j]) for j=1:length(log2_sizes)])
    tc0_, sc0_ = timespace_complexity(NestedEinsum(tree), size_dict)
    @test tc0 ≈ tc0_ && sc0 ≈ sc0_
    opt_tree = optimize_subtree!(copy(tree), 20.0, log2_sizes, 5, 2.0)
    tc1, sc1 = timespace_complexity(opt_tree, exp2.(log2_sizes))
    @test sc1 < sc0
end