using OMEinsum, OMEinsumContractionOrders
using Test, Random, Graphs
using KaHyPar
using OMEinsum: NestedEinsum, DynamicEinCode

@testset "interface" begin
    function random_regular_eincode(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [minmax(e.src,e.dst) for e in Graphs.edges(g)]
        return EinCode((ixs..., [(i,) for i in Graphs.vertices(g)]...), ())
    end
    Random.seed!(2)
    g = random_regular_graph(100, 3)
    code = random_regular_eincode(100, 3)
    xs = [[randn(2,2) for i=1:150]..., [randn(2) for i=1:100]...]

    results = Float64[]
    for optimizer in [TreeSA(ntrials=1), TreeSA(ntrials=1, nslices=5), GreedyMethod(), KaHyParBipartite(sc_target=18), SABipartite(sc_target=18, ntrials=1)]
        for simplifier in (nothing, MergeVectors(), MergeGreedy())
            @info "optimizer = $(optimizer), simplifier = $(simplifier)"
            res = optimize_code(code,uniformsize(code, 2), optimizer, simplifier)
            tc, sc = OMEinsum.timespace_complexity(res, uniformsize(code, 2))
            @test sc <= 18
            push!(results, res(xs...)[])
        end
    end
    for i=1:length(results)-1
        @test results[i] ≈ results[i+1]
    end
end

@testset "corner case: smaller contraction orders" begin
    code = ein"i->"
    sizes = uniformsize(code, 2)
    ne = NestedEinsum((1,), code)
    dne = NestedEinsum((1,), DynamicEinCode(code))
    @test optimize_code(code, sizes, GreedyMethod()) == ne
    @test optimize_code(code, sizes, TreeSA()) == SlicedEinsum(Slicing(Char[]), dne)
    @test optimize_code(code, sizes, TreeSA(nslices=2)) == SlicedEinsum(Slicing(Char[]), dne)
    @test optimize_code(code, sizes, KaHyParBipartite(sc_target=25)) == dne
    @test optimize_code(code, sizes, SABipartite(sc_target=25)) == dne
end

@testset "peak memory" begin
    Random.seed!(2)
    code = ein"(ab,a),ac->bc"
    @test peak_memory(code, uniformsize(code, 5)) == 75

    function random_regular_eincode(n, k)
        g = Graphs.random_regular_graph(n, k)
        ixs = [minmax(e.src,e.dst) for e in Graphs.edges(g)]
        return EinCode((ixs..., [(i,) for i in Graphs.vertices(g)]...), ())
    end
    code = random_regular_eincode(50, 3)
    @test peak_memory(code, uniformsize(code, 5)) == (25 * 75 + 5 * 50)
    code1 = optimize_code(code, uniformsize(code, 5), GreedyMethod())
    pm1 = peak_memory(code1, uniformsize(code, 5))
    tc1, sc1, rw1 = timespacereadwrite_complexity(code1, uniformsize(code, 5))
    code2 = optimize_code(code, uniformsize(code, 5), TreeSA(ntrials=1, nslices=5))
    pm2 = peak_memory(code2, uniformsize(code, 5))
    tc2, sc2, rw2 = timespacereadwrite_complexity(code2, uniformsize(code, 5))
    @test 5 * 2^sc1 > pm1 > 2^sc1
    @test 5 * 2^sc2 > pm2 > 2^sc2
end