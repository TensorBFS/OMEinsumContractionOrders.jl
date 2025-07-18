using Graphs
using Test, Random
using SparseArrays
using OMEinsumContractionOrders
using OMEinsumContractionOrders: get_coarse_grained_graph, _connected_components, bipartite_sc, group_sc, coarse_grained_optimize,
    map_tree_to_parts, ContractionTree, optimize_greedy, optimize_kahypar
using KaHyPar
using OMEinsum: decorate

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

    res = coarse_grained_optimize(adj, parts, ones(6), GreedyMethod())
    @test res == ContractionTree(ContractionTree(1,2), 3)
    @test map_tree_to_parts(res, [[[1,2], 3], [7,6], [9, [4,1]]]) == [[[[1,2], 3], [7,6]], [9, [4,1]]]
end

@testset "kahypar" begin
    Random.seed!(2)
    function random_regular_eincode(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [[minmax(e.src,e.dst)...] for e in Graphs.edges(g)]
        return OMEinsumContractionOrders.EinCode([ixs..., [[i] for i in Graphs.vertices(g)]...], Int[])
    end

    function random_regular_open_eincode(n, k, m)
        g = Graphs.random_regular_graph(n, k)
        ixs = [[minmax(e.src,e.dst)...] for e in Graphs.edges(g)]
        iy = Int[]
        while length(iy) < m
            v = rand(1:n)
            if !(v in iy)
                push!(iy, v)
            end
        end
        sort!(iy)
        return OMEinsumContractionOrders.EinCode([ixs..., [[i] for i in Graphs.vertices(g)]...], iy)
    end

    g = random_regular_graph(220, 3)
    rows = Int[]
    cols = Int[]
    for (i,edge) in enumerate(edges(g))
        push!(rows, edge.src, edge.dst)
        push!(cols, i, i)
    end
    graph = sparse(rows, cols, ones(Int, length(rows)))
    sc_target = 30.0
    log2_sizes = fill(1, size(graph, 2))
    b = KaHyParBipartite(sc_target=sc_target, imbalances=[0.0:0.02:0.8...])
    res = bipartite_sc(b, graph, collect(1:size(graph, 1)), log2_sizes)
    @test group_sc(graph, res.part1, log2_sizes) <= sc_target
    @test group_sc(graph, res.part2, log2_sizes) <= sc_target
    sc_target = 30.0
    res = bipartite_sc(b, graph, collect(1:size(graph, 1)), log2_sizes)
    @test group_sc(graph, res.part1, log2_sizes) <= sc_target
    @test group_sc(graph, res.part2, log2_sizes) <= sc_target

    code = random_regular_eincode(220, 3)
    res = optimize_kahypar(code,uniformsize(code, 2); max_group_size=50, sc_target=30)
    cc = contraction_complexity(res, uniformsize(code, 2))
    @test cc.sc <= 30

    # contraction test
    code = random_regular_eincode(50, 3)
    codeg = optimize_kahypar(code, uniformsize(code, 2); max_group_size=10, sc_target=10)
    codet = optimize_kahypar(code, uniformsize(code, 2); max_group_size=10, sc_target=10, sub_optimizer = TreeSA())
    codek = optimize_greedy(code, uniformsize(code, 2); α=0.0, temperature=0.0)

    cc_kg = contraction_complexity(codeg, uniformsize(code, 2))
    cc_kt = contraction_complexity(codet, uniformsize(code, 2))
    cc_g = contraction_complexity(codek, uniformsize(code, 2))

    @test cc_kg.sc <= 12
    @test cc_kt.sc <= 12
    @test cc_g.sc <= 12

    xs = [[2*randn(2, 2) for i=1:75]..., [randn(2) for i=1:50]...]
    resg = decorate(codeg)(xs...)
    resk = decorate(codek)(xs...)
    @test resg ≈ resk

    Random.seed!(2)
    code = random_regular_open_eincode(50, 3, 3)
    codeg = optimize_kahypar(code, uniformsize(code, 2); max_group_size=10, sc_target=10)
    codet = optimize_kahypar(code, uniformsize(code, 2); max_group_size=10, sc_target=10, sub_optimizer = TreeSA())
    codek = optimize_greedy(code, uniformsize(code, 2); α=0.0, temperature=0.0)

    xs = [[2*randn(2, 2) for i=1:75]..., [randn(2) for i=1:50]...]
    resg = decorate(codeg)(xs...)
    rest = decorate(codet)(xs...)
    resk = decorate(codek)(xs...)
    @test rest ≈ resk
    @test resg ≈ resk
end

@testset "kahypar no error" begin
    function random_regular_open_eincode(n, k, m)
        g = Graphs.random_regular_graph(n, k)
        ixs = [[minmax(e.src,e.dst)...] for e in Graphs.edges(g)]
        iy = Int[]
        while length(iy) < m
            v = rand(1:n)
            if !(v in iy)
                push!(iy, v)
            end
        end
        sort!(iy)
        return OMEinsumContractionOrders.EinCode([ixs..., [[i] for i in Graphs.vertices(g)]...], iy)
    end


    Random.seed!(2)
    code = random_regular_open_eincode(50, 3, 3)
    codeg = optimize_kahypar(code, uniformsize(code, 2); max_group_size=10, sc_target=5)
    codek = optimize_greedy(code, uniformsize(code, 2); α=0.0, temperature=0.0)

    xs = [[2*randn(2, 2) for i=1:75]..., [randn(2) for i=1:50]...]
    resg = decorate(codeg)(xs...)
    resk = decorate(codek)(xs...)
    @test resg ≈ resk
end
