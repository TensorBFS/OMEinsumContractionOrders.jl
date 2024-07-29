using OMEinsumContractionOrders
using OMEinsumContractionOrders: IncidenceList, optimize_exact_treewidth

using Test, Random
@testset "tree width" begin
    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], Vector{Char}())
    size_dict = Dict([c=>(1<<i) for (i,c) in enumerate(['a', 'b', 'c', 'd', 'e', 'f'])]...)
    optimizer = ExactTreewidth()
    optcode = optimize_exact_treewidth(optimizer, eincode, size_dict)
    cc = contraction_complexity(optcode, size_dict)
    # test flop
    @test cc.tc ≈ log2(flop(optcode, size_dict))
    @test 16 <= cc.tc <= log2(exp2(10)+exp2(16)+exp2(15)+exp2(9))
    @test cc.sc == 11

    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], ['a'])
    optcode2 = optimize_exact_treewidth(optimizer, eincode, size_dict)
    cc = contraction_complexity(optcode, size_dict)
    @test cc.sc == 11

    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e'], ['e'], ['f']], ['a', 'f'])
    optcode3 = optimize_exact_treewidth(optimizer, eincode, size_dict)
    cc = contraction_complexity(optcode, size_dict)
    @test cc.sc == 11
end