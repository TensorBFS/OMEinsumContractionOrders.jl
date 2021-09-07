using Test, OMEinsumContractionOrders
using LightGraphs, Random

@testset "sa bipartition" begin
    Random.seed!(3)
    g = random_regular_graph(120, 3)
    adj, = OMEinsumContractionOrders.adjacency_matrix([(e.src,e.dst) for e in edges(g)])
    ws = fill(log2(2), ne(g))
    vertices = 1:110
    best, tc = OMEinsumContractionOrders.divide_group_sa(adj, vertices; βs=0.1:0.2:20.0, niters=1000, ntrials=100, sc_target=40, log2_sizes=ws)
    @test best isa OMEinsumContractionOrders.PartitionState
end


@testset "sa bipartition on 220 node closed graph" begin
    Random.seed!(3)
    g = random_regular_graph(120, 3)
    adj, = OMEinsumContractionOrders.adjacency_matrix([(e.src,e.dst) for e in edges(g)])
    ws = fill(log2(2), ne(g))
    vertices = 1:120
    best, tc = OMEinsumContractionOrders.divide_group_sa(adj, vertices; βs=0.1:0.2:20.0, niters=1000, ntrials=100, sc_target=40, log2_sizes=ws)
    @test best isa OMEinsumContractionOrders.PartitionState
end
