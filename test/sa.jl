using Test, OMEinsumContractionOrders, OMEinsum
using LightGraphs, Random
using SparseArrays

@testset "sa bipartition" begin
    Random.seed!(3)
    g = random_regular_graph(120, 3)
    adj, = OMEinsumContractionOrders.adjacency_matrix([(e.src,e.dst) for e in edges(g)])
    ws = fill(log2(2), ne(g))
    vertices = 1:110
    b = SABipartite(βs=0.1:0.2:20.0, niters=1000, ntrials=100, sc_target=40)
    g1, g2 = OMEinsumContractionOrders.bipartite_sc(b, adj, vertices, ws)
    @test length(g1) + length(g2) == 110
end

@testset "sa" begin
    Random.seed!(2)
    function random_regular_eincode(n, k)
        g = LightGraphs.random_regular_graph(n, k)
        ixs = [minmax(e.src,e.dst) for e in LightGraphs.edges(g)]
        return EinCode((ixs..., [(i,) for i in LightGraphs.vertices(g)]...), ())
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
    b = SABipartite(βs=0.1:0.2:20.0, niters=1000, ntrials=100, sc_target=28)
    group1, group2 = OMEinsumContractionOrders.bipartite_sc(b, graph, collect(1:size(graph, 1)), log2_sizes)
    @test OMEinsumContractionOrders.group_sc(graph, group1, log2_sizes) <= sc_target
    @test OMEinsumContractionOrders.group_sc(graph, group2, log2_sizes) <= sc_target
    sc_target = 27.0
    group11, group12 = OMEinsumContractionOrders.bipartite_sc(b, graph, group1, log2_sizes)
    @test OMEinsumContractionOrders.group_sc(graph, group11, log2_sizes) <= sc_target
    @test OMEinsumContractionOrders.group_sc(graph, group12, log2_sizes) <= sc_target

    code = random_regular_eincode(220, 3)
    res = optimize_sa(code,uniformsize(code, 2); sc_target=30)
    tc, sc = OMEinsum.timespace_complexity(res, uniformsize(code, 2))
    @test sc <= 30

    # contraction test
    code = random_regular_eincode(50, 3)
    codeg = optimize_sa(code, uniformsize(code, 2); sc_target=12)
    codek = optimize_greedy(code, uniformsize(code, 2))
    tc, sc = OMEinsum.timespace_complexity(codek, uniformsize(code, 2))
    @test sc <= 12
    xs = [[2*randn(2, 2) for i=1:75]..., [randn(2) for i=1:50]...]
    resg = codeg(xs...)
    resk = codek(xs...)
    @test resg ≈ resk
end

