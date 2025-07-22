using OMEinsumContractionOrders
using OMEinsumContractionOrders: IncidenceList, optimize_treewidth, getixsv
using OMEinsumContractionOrders: BFS, MCS, LexBFS, RCMMD, RCMGL, MCSM, LexM, AMF, MF, MMD, BT, SafeRules
using OMEinsum: decorate
using Test, Random, JSON

@testset "tree width" begin

    optimizer = ExactTreewidth()
    size_dict = Dict([c=>(1<<i) for (i,c) in enumerate(['a', 'b', 'c', 'd', 'e', 'f'])]...)

    # eincode with no open edges
    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], Vector{Char}())
    tensors = [rand([size_dict[j] for j in ixs]...) for ixs in getixsv(eincode)]
    optcode = optimize_treewidth(optimizer, eincode, size_dict)
    cc = contraction_complexity(optcode, size_dict)
    # test flop
    @test cc.tc ≈ log2(flop(optcode, size_dict))
    @test (16 <= cc.tc <= log2(exp2(10)+exp2(16)+exp2(15)+exp2(9))) | (cc.tc ≈ log2(exp2(10)+exp2(16)+exp2(15)+exp2(9)))
    @test cc.sc == 11
    @test decorate(eincode)(tensors...) ≈ decorate(optcode)(tensors...)

    # eincode with open edges
    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], ['a'])
    tensors = [rand([size_dict[j] for j in ixs]...) for ixs in getixsv(eincode)]
    optcode = optimize_treewidth(optimizer, eincode, size_dict)
    cc = contraction_complexity(optcode, size_dict)
    @test cc.sc == 11
    @test decorate(eincode)(tensors...) ≈ decorate(optcode)(tensors...)

    # disconnect contraction tree
    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e'], ['e'], ['f']], ['a', 'f'])
    tensors = [rand([size_dict[j] for j in ixs]...) for ixs in getixsv(eincode)]
    optcode = optimize_treewidth(optimizer, eincode, size_dict)
    cc = contraction_complexity(optcode, size_dict)
    @test cc.sc == 7
    @test decorate(eincode)(tensors...) ≈ decorate(optcode)(tensors...)

    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e'], ['e'], ['f'], Char[]], ['a', 'f'])
    tensors = tensors ∪ [fill(2.0,())]
    optcode = optimize_treewidth(optimizer, eincode, size_dict)
    cc = contraction_complexity(optcode, size_dict)
    @test decorate(eincode)(tensors...) ≈ decorate(optcode)(tensors...)
end

@testset "tcsc" begin
    optimizer = ExactTreewidth()
    size_dict = Dict([c=>2 for (i,c) in enumerate(['a' + (i-1) for i in 1:10])]...)

    eincode = OMEinsumContractionOrders.EinCode([['a', 'b', 'c', 'd', 'e', 'f'], ['e', 'f', 'g', 'h', 'i', 'j'], ['a', 'b', 'c', 'd', 'g', 'h', 'i', 'j']], Vector{Char}())
    ixs = getixsv(eincode)
    tensors = [rand([size_dict[j] for j in ix]...) for ix in ixs]
    optcode = optimize_treewidth(optimizer, eincode, size_dict)
    incidence_list = IncidenceList(Dict([i=>ix for (i,ix) in enumerate(ixs)] ∪ [(length(ixs) + 1 => eincode.iy)]))
    lg = OMEinsumContractionOrders.il2lg(incidence_list, collect(keys(size_dict)))
    tw = OMEinsumContractionOrders.TreeWidthSolver.exact_treewidth(lg)
    cc = contraction_complexity(optcode, size_dict)
    @test isapprox(cc.tc, tw + 1, atol=0.9)
end

@testset "treewidth" begin
    for alg in [SafeRules(BT()), MF(), MMD(), AMF(), LexM(), LexBFS(), BFS(), MCS(), RCMMD(), RCMGL(), MCSM()]
        optimizer = Treewidth(alg=alg)
        size_dict = Dict([c=>(1<<i) for (i,c) in enumerate(['a', 'b', 'c', 'd', 'e', 'f'])]...)

        # eincode with no open edges
        eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], Vector{Char}())
        tensors = [rand([size_dict[j] for j in ixs]...) for ixs in getixsv(eincode)]
        optcode = optimize_treewidth(optimizer, eincode, size_dict)
        cc = contraction_complexity(optcode, size_dict)
        # test flop
        @test cc.tc ≈ log2(flop(optcode, size_dict))
        @test (16 <= cc.tc <= log2(exp2(10)+exp2(16)+exp2(15)+exp2(9))) | (cc.tc ≈ log2(exp2(10)+exp2(16)+exp2(15)+exp2(9)))
        @test cc.sc == 11
        @test decorate(eincode)(tensors...) ≈ decorate(optcode)(tensors...)
    end
end

@testset "fix issue 87" begin
    file = "rnd_mixed_07.json"
    dict = JSON.parsefile(joinpath(@__DIR__, file))

    einsum = dict["einsum"]

    ixs = Vector{Vector{Int}}(einsum["ixs"])
    iy = Vector{Int}(einsum["iy"])

    size_dict = Dict(
        parse(Int, key) => val
        for (key, val) in dict["size"]
    )

    code = OMEinsumContractionOrders.EinCode(ixs, iy)
    optimizer = Treewidth(; alg=AMF())
    optcode = optimize_code(code, size_dict, optimizer);
    @test optcode isa OMEinsumContractionOrders.NestedEinsum
end