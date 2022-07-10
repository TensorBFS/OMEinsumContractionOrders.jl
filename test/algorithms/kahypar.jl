using Graphs
using Test, Random
using SparseArrays
using OMEinsumContractionOrders.ContractionOrderAlgorithms
using OMEinsumContractionOrders.ContractionOrderAlgorithms: get_coarse_grained_graph, _connected_components, bipartite_sc, group_sc, coarse_graphed_optimize,
    map_tree_to_parts
using KaHyPar

@testset "graph coarse graining" begin
    Random.seed!(2)
    adj = zeros(5, 6)
    for ind in [[1,1], [1,3], [2,2], [3,2], [4,4], [4,5], [4,3], [5, 4]]
        adj[ind...] = 1
    end
    parts = [[1,3], [2], [4]]
    incidence_list = get_coarse_grained_graph(sparse(adj), parts)
    @test incidence_list.v2e[1] == [1,2,3]
    @test incidence_list.v2e[2] == [2]
    @test incidence_list.v2e[3] == [3,4,5]
    @test incidence_list.openedges == [4]

    @test length(_connected_components(adj, parts[1])) == 2

    res = coarse_grained_optimize(adj, parts, ones(6), GreedyMethod(OMEinsum.MinSpaceOut(), 10))
    @test res == OMEinsum.ContractionTree(OMEinsum.ContractionTree(1,2), 3)
    @test map_tree_to_parts(res, [[[1,2], 3], [7,6], [9, [4,1]]]) == [[[[1,2], 3], [7,6]], [9, [4,1]]]
end

@testset "kahypar" begin
    Random.seed!(2)
    function random_regular_eincode(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [minmax(e.src,e.dst) for e in Graphs.edges(g)]
        return EinCode((ixs..., [(i,) for i in Graphs.vertices(g)]...), ())
    end

    g = random_regular_graph(220, 3)
    rows = Int[]
    cols = Int[]
    for (i,edge) in enumerate(edges(g))
        push!(rows, edge.src, edge.dst)
        push!(cols, i, i)
    end
    graph = sparse(rows, cols, ones(Int, length(rows)))
    sc_target = 28.0
    log2_sizes = fill(1, size(graph, 2))
    b = KaHyParBipartite(sc_target=sc_target, imbalances=[0.0:0.02:0.8...])
    group1, group2 = bipartite_sc(b, graph, collect(1:size(graph, 1)), log2_sizes)
    @test group_sc(graph, group1, log2_sizes) <= sc_target
    @test group_sc(graph, group2, log2_sizes) <= sc_target
    sc_target = 27.0
    group11, group12 = bipartite_sc(b, graph, group1, log2_sizes)
    @test group_sc(graph, group11, log2_sizes) <= sc_target
    @test group_sc(graph, group12, log2_sizes) <= sc_target

    code = random_regular_eincode(220, 3)
    res = optimize_kahypar(code,uniformsize(code, 2); max_group_size=50, sc_target=30)
    tc, sc = OMEinsum.timespace_complexity(res, uniformsize(code, 2))
    @test sc <= 30

    # contraction test
    code = random_regular_eincode(50, 3)
    codeg = optimize_kahypar(code, uniformsize(code, 2); max_group_size=10, sc_target=12)
    codek = optimize_greedy(code, uniformsize(code, 2))
    tc, sc = OMEinsum.timespace_complexity(codek, uniformsize(code, 2))
    @test sc <= 12
    xs = [[2*randn(2, 2) for i=1:75]..., [randn(2) for i=1:50]...]
    resg = codeg(xs...)
    resk = codek(xs...)
    @test resg ≈ resk

    Random.seed!(2)
    code = random_regular_eincode(220, 3)
    codeg_auto = optimize_kahypar_auto(code, uniformsize(code, 2), greedy_method=OMEinsum.MinSpaceOut())
    tc, sc = OMEinsum.timespace_complexity(codeg_auto, uniformsize(code, 2))
    @test sc <= 30
end

@testset "regression test" begin
    code = ein"i->"
    optcode = optimize_kahypar(code, Dict('i'=>4), sc_target=10, max_group_size=10)
    @test optcode isa NestedEinsum
    x = randn(4)
    @test optcode(x) ≈ code(x)

    code = ein"i,j->"
    optcode = optimize_kahypar(code, Dict('i'=>4, 'j'=>4), sc_target=10, max_group_size=10)
    @test optcode isa NestedEinsum
    x = randn(4)
    y = randn(4)
    @test optcode(x, y) ≈ code(x, y)

    code = ein"ij,jk,kl->ijl"
    optcode = optimize_kahypar(code, Dict('i'=>4, 'j'=>4, 'k'=>4, 'l'=>4), sc_target=4, max_group_size=2)
    @test optcode isa NestedEinsum
    a, b, c = [rand(4,4) for i=1:4]
    @test optcode(a, b, c) ≈ code(a, b, c)

    code = ein"ij,jk,kl->ijl"
    optcode = optimize_kahypar(code, Dict('i'=>3, 'j'=>3, 'k'=>3, 'l'=>3), sc_target=4, max_group_size=2)
    @test optcode isa NestedEinsum
    a, b, c = [rand(3,3) for i=1:4]
    @test optcode(a, b, c) ≈ code(a, b, c)
end
