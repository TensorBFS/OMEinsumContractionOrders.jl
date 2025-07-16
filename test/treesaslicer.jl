using OMEinsumContractionOrders, Test, Random
using OMEinsumContractionOrders: random_exprtree, ExprTree, ExprInfo,
    ruleset, update_tree!, tcscrw, optimize_subtree!, optimize_tree_sa!, labels, tree_timespace_complexity, fast_log2sumexp2,
    ExprTree, _label_dict, Slicer, slice_tree, SlicedEinsum,
    optimize_greedy, optimize_tree, optimize_hyper_nd
using Graphs
using OMEinsum: decorate
using KaHyPar

@testset "slicing" begin
    function random_regular_eincode(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [[minmax(e.src,e.dst)...] for e in Graphs.edges(g)]
        return OMEinsumContractionOrders.EinCode([ixs..., [[i] for i in Graphs.vertices(g)]...], Int[])
    end
    # contraction test
    Random.seed!(2)
    code = random_regular_eincode(100, 3)

    for code0 in [optimize_greedy(code, uniformsize(code, 2)),
            optimize_tree(code, uniformsize(code, 2); initializer=:greedy, greedy_method=GreedyMethod(nrepeat=1), sc_target=10, βs=0.1:0.05:20.0, ntrials=2, niters=10, sc_weight=1.0, rw_weight=1.0),
            optimize_hyper_nd(HyperND(), code, uniformsize(code, 2))]

        codet = slice_tree(code0, uniformsize(code, 2); sc_target = 10)
        @test codet isa SlicedEinsum

        cc0 = contraction_complexity(code0, uniformsize(code, 2))
        cc = contraction_complexity(codet, uniformsize(code, 2))

        @show cc.tc, cc.sc, cc0.tc, cc0.sc
        @show length(codet.slicing)
        @test cc.sc <= 10

        xs = [[2*randn(2, 2) for i=1:150]..., [randn(2) for i=1:100]...]
        res0 = decorate(code0)(xs...)
        rest = decorate(codet)(xs...)
        @test res0 ≈ rest
    end

    # with open edges
    Random.seed!(2)
    code = OMEinsumContractionOrders.EinCode(random_regular_eincode(100, 3).ixs, [3,81,2])

    for code0 in [optimize_greedy(code, uniformsize(code, 2)), optimize_tree(code, uniformsize(code, 2)), optimize_hyper_nd(HyperND(), code, uniformsize(code, 2))]
        codet = slice_tree(code0, uniformsize(code, 2); sc_target = 10, ntrials = 10)
        @test codet isa SlicedEinsum

        cc0 = contraction_complexity(code0, uniformsize(code, 2))
        cc = contraction_complexity(codet, uniformsize(code, 2))

        @show cc.tc, cc.sc, cc0.tc, cc0.sc
        @test cc.sc <= 10

        xs = [[2*randn(2, 2) for i=1:150]..., [randn(2) for i=1:100]...]
        res0 = decorate(code0)(xs...)
        rest = decorate(codet)(xs...)
        @test res0 ≈ rest
    end

    # slice with fixed slices
    Random.seed!(2)
    code = OMEinsumContractionOrders.EinCode(random_regular_eincode(100, 3).ixs, [3,10,2])
    for code0 in [optimize_greedy(code, uniformsize(code, 2)), optimize_tree(code, uniformsize(code, 2)), optimize_hyper_nd(HyperND(), code, uniformsize(code, 2))]
        code1 = slice_tree(code0, uniformsize(code, 2); sc_target = 10, fixed_slices=[10, 89, 26, 3, 51])
        code2 = slice_tree(code0, uniformsize(code, 2); sc_target = 10, fixed_slices=[10, 89])

        @test length(code1.slicing) >= 5 && Set([10, 89, 26, 3, 51]) ⊆ Set(code1.slicing)
        @test length(code2.slicing) >= 2 && Set([10, 89]) ⊆ Set(code2.slicing)
        @test contraction_complexity(code1, uniformsize(code, 2)).sc <= 10
        @test contraction_complexity(code2, uniformsize(code, 2)).sc <= 10

        xs = [[2*randn(2, 2) for i=1:150]..., [randn(2) for i=1:100]...]
        res0 = decorate(code0)(xs...)
        res1 = decorate(code1)(xs...)
        res2 = decorate(code2)(xs...)
        @test res0 ≈ res1
        @test res0 ≈ res2
    end

    function random_regular_eincode_char(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [[minmax('0' + e.src, '0'+e.dst)...] for e in Graphs.edges(g)]
        return OMEinsumContractionOrders.EinCode([ixs..., [['0'+i] for i in Graphs.vertices(g)]...], Char[])
    end
    code = OMEinsumContractionOrders.EinCode(random_regular_eincode_char(20, 3).ixs, ['3','8','2'])
    code0 = optimize_tree(code, uniformsize(code, 2))
    code1 = slice_tree(code0, uniformsize(code, 2); sc_target = 5, fixed_slices=['7'])
    @test eltype(code1.eins.eins.iy) == Char

    xs = [[2*randn(2, 2) for i=1:30]..., [randn(2) for i=1:20]...]
    res0 = decorate(code0)(xs...)
    res1 = decorate(code1)(xs...)
    @test res0 ≈ res1
end
