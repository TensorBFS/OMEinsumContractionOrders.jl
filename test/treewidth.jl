using OMEinsumContractionOrders
using OMEinsumContractionOrders: IncidenceList, optimize_exact_treewidth, getixsv
using OMEinsum: decorate

using Test, Random
@testset "tree width" begin

    optimizer = ExactTreewidth()
    size_dict = Dict([c=>(1<<i) for (i,c) in enumerate(['a', 'b', 'c', 'd', 'e', 'f'])]...)

    # eincode with no open edges
    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], Vector{Char}())
    tensors = [rand([size_dict[j] for j in ixs]...) for ixs in getixsv(eincode)]
    optcode = optimize_exact_treewidth(optimizer, eincode, size_dict)
    cc = contraction_complexity(optcode, size_dict)
    # test flop
    @test cc.tc ≈ log2(flop(optcode, size_dict))
    @test 16 <= cc.tc <= log2(exp2(10)+exp2(16)+exp2(15)+exp2(9))
    @test cc.sc == 11
    @test decorate(eincode)(tensors...) ≈ decorate(optcode)(tensors...)

    # eincode with open edges
    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], ['a'])
    tensors = [rand([size_dict[j] for j in ixs]...) for ixs in getixsv(eincode)]
    optcode = optimize_exact_treewidth(optimizer, eincode, size_dict)
    cc = contraction_complexity(optcode, size_dict)
    @test cc.sc == 11
    @test decorate(eincode)(tensors...) ≈ decorate(optcode)(tensors...)

    # disconnect contraction tree
    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e'], ['e'], ['f']], ['a', 'f'])
    tensors = [rand([size_dict[j] for j in ixs]...) for ixs in getixsv(eincode)]
    optcode = optimize_exact_treewidth(optimizer, eincode, size_dict)
    cc = contraction_complexity(optcode, size_dict)
    @test cc.sc == 7
    @test decorate(eincode)(tensors...) ≈ decorate(optcode)(tensors...)

    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e'], ['e'], ['f'], Char[]], ['a', 'f'])
    tensors = tensors ∪ fill(2.0,())
    optcode = optimize_exact_treewidth(optimizer, eincode, size_dict)
    cc = contraction_complexity(optcode, size_dict)
    @test decorate(eincode)(tensors...) ≈ decorate(optcode)(tensors...)
end