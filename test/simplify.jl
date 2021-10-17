using OMEinsumContractionOrders, OMEinsum
using Test, Random

@testset "simplify vectors" begin
    tn = ein"ab,bc,cd,de,a,b,c,d,f,f,e->ab"
    xs = [randn(fill(2, length(ix))...) for ix in OMEinsum.getixs(tn)]
    simplifier, tn2 = merge_vectors(tn)
    @test tn2 == ein"ab,bc,cd,de,f->ab"
    @test tn(xs...) ≈ tn2(simplifier(xs...)...)
    tn3 = optimize_greedy(tn2, uniformsize(tn2, 2))
    tn4 = embed_simplifier(tn3, simplifier)
    @test tn(xs...) ≈ tn4(xs...)

    tn = OMEinsum.DynamicEinCode(tn)
    xs = [randn(fill(2, length(ix))...) for ix in OMEinsum.getixs(tn)]
    simplifier, tn2 = merge_vectors(tn)
    @test tn2 == OMEinsum.DynamicEinCode(ein"ab,bc,cd,de,f->ab")
    @test tn(xs...) ≈ tn2(simplifier(xs...)...)
    tn3 = optimize_greedy(tn2, uniformsize(tn2, 2))
    tn4 = embed_simplifier(tn3, simplifier)
    @test tn(xs...) ≈ tn4(xs...)
end
