using OMEinsumContractionOrders, Test, Random
using OMEinsumContractionOrders: random_exprtree, exprtree, ExprTree, ExprInfo, ruleset, update_tree!, tcsc
using OMEinsum, LightGraphs

@testset "random expr tree" begin
    Random.seed!(2)
    tree = random_exprtree([[1,2,5], [2,3], [2,4]], [5], 5)
    @test tree isa OMEinsumContractionOrders.ExprTree
    tree2 = random_exprtree(EinCode(((1,2,5), (2,3), (2,4)), (5,)))
    @test tree isa OMEinsumContractionOrders.ExprTree
    function random_regular_eincode(n, k)
        g = LightGraphs.random_regular_graph(n, k)
        ixs = [minmax(e.src,e.dst) for e in LightGraphs.edges(g)]
        return EinCode((ixs..., [(i,) for i in LightGraphs.vertices(g)]...), ())
    end
    code = random_regular_eincode(20, 3)
    optcode = optimize_greedy(code, uniformsize(code, 2))
    tree3 = exprtree(optcode)
    @test tree isa OMEinsumContractionOrders.ExprTree
end

@testset "rules" begin
    t1 = ExprTree([1,2], ExprTree([2,3], [1,4], ExprInfo([1,2])), ExprInfo([2]))
    t2 = ExprTree(ExprTree([2,3], [1,4], ExprInfo([1,2,3])), [1,2], ExprInfo([2]))
    t3 = ExprTree([2,3], [1,2], ExprInfo([2]))
    t4 = ExprTree(ExprTree([2,3], [1,4], ExprInfo([1,2])), ExprTree([5,1], [1], ExprInfo([1])), ExprInfo([2]))
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
    @test t11 == ExprTree([2,3], ExprTree([1,2], [1,4], ExprInfo([2])), ExprInfo([2]))
    t11_ = update_tree!(copy(t1), 4)
    @test t11_ == ExprTree([1,4], ExprTree([2,3], [1,2], ExprInfo([1,2])), ExprInfo([2]))
    t22 = update_tree!(copy(t2), 1)
    @test t22 == ExprTree(ExprTree([2,3], [1,2], ExprInfo([1,2])), [1,4], ExprInfo([2]))
    t22_ = update_tree!(copy(t2), 2)
    @test t22_ == ExprTree(ExprTree([1,2], [1,4], ExprInfo([2])), [2,3], ExprInfo([2]))
    t44 = update_tree!(copy(t4), 1)
    @test t44 == ExprTree(ExprTree([2,3], ExprTree([5,1], [1], ExprInfo([1])), ExprInfo([1,2])), [1,4], ExprInfo([2]))
end