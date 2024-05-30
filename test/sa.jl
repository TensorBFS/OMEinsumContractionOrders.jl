using Test, OMEinsumContractionOrders
using Graphs, Random
using SparseArrays
using OMEinsumContractionOrders: bipartite_sc, adjacency_matrix, SABipartite, group_sc, bipartite_sc, optimize_sa, optimize_greedy
using OMEinsum: decorate

@testset "sa bipartition" begin
    Random.seed!(3)
    g = random_regular_graph(120, 3)
    adj, = adjacency_matrix([(e.src,e.dst) for e in edges(g)])
    ws = fill(log2(2), ne(g))
    vertices = 1:110
    b = SABipartite(βs=0.1:0.2:20.0, niters=1000, ntrials=100, sc_target=40, initializer=:random)
    g1, g2 = bipartite_sc(b, adj, vertices, ws)
    @test length(g1) + length(g2) == 110
end

@testset "sa" begin
    Random.seed!(2)
    function random_regular_eincode(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [[minmax(e.src,e.dst)...] for e in Graphs.edges(g)]
        return OMEinsumContractionOrders.EinCode([ixs..., [[i] for i in Graphs.vertices(g)]...], Int[])
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
    βs = 0.01:0.05:10.0
    b = SABipartite(βs=βs, niters=1000, ntrials=100, sc_target=28)
    group1, group2 = bipartite_sc(b, graph, collect(1:size(graph, 1)), log2_sizes)
    @test group_sc(graph, group1, log2_sizes) <= sc_target+2
    @test group_sc(graph, group2, log2_sizes) <= sc_target+2
    sc_target = 27.0
    group11, group12 = bipartite_sc(b, graph, group1, log2_sizes)
    @test group_sc(graph, group11, log2_sizes) <= sc_target+2
    @test group_sc(graph, group12, log2_sizes) <= sc_target+2

    code = random_regular_eincode(220, 3)
    res = optimize_sa(code,uniformsize(code, 2); sc_target=30, βs=βs)
    cc = contraction_complexity(res, uniformsize(code, 2))
    @test cc.sc <= 32

    tc1, sc1, rw1 = timespacereadwrite_complexity(res, uniformsize(code, 2))
    cc = contraction_complexity(res, uniformsize(code, 2))
    @test (tc1, sc1, rw1) == (cc...,)
    println(cc)

    # contraction test
    code = random_regular_eincode(50, 3)
    codeg = optimize_sa(code, uniformsize(code, 2); sc_target=12, βs=βs, ntrials=1, initializer=:greedy)
    codek = optimize_greedy(code, uniformsize(code, 2))
    cc = contraction_complexity(codek, uniformsize(code, 2))
    @test cc.sc <= 12
    xs = [[2*randn(2, 2) for i=1:75]..., [randn(2) for i=1:50]...]
    resg = decorate(codeg)(xs...)
    resk = decorate(codek)(xs...)
    @test resg ≈ resk
end

