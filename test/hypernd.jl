using Graphs
using OMEinsumContractionOrders
using OMEinsumContractionOrders: optimize_hyper_nd
using Random
using Test

import KaHyPar

@testset "chain and ring" begin
    # chain
    code = OMEinsumContractionOrders.EinCode([[1,2], [2,3], [3,4], [4,5], [5,6], [6,7], [7,8], [8,9], [9,10]], Int[1, 10])
    size_dict = Dict([i=>2 for i in 1:10])
    optcode = optimize_hyper_nd(HyperND(; width=5), code, size_dict)
    cc = contraction_complexity(optcode, size_dict)
    @test cc.sc == 2

    # ring
    ring = OMEinsumContractionOrders.EinCode([[1,2], [2,3], [3,4], [4,5], [5,6], [6,7], [7,8], [8,9], [9,10], [10,1]], Int[])
    optcode = optimize_hyper_nd(HyperND(; width=5), ring, size_dict)
    cc = contraction_complexity(optcode, size_dict)
    @test cc.sc == 2

    # petersen
    graph = smallgraph(:tutte)
    code = OMEinsumContractionOrders.EinCode([[e.src, e.dst] for e in edges(graph)], Int[])
    size_dict = Dict([i=>2 for i in 1:nv(graph)])
    optcode = optimize_hyper_nd(HyperND(; width=5), code, size_dict)
    cc = contraction_complexity(optcode, size_dict)
    @test cc.sc == 5
end
