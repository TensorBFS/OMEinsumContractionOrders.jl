using OMEinsumContractionOrders, OMEinsum
using Test, Random

@testset "simplify vectors" begin
    tn = ein"ab,bc,cd,de,a,b,c,d,f,f->ab"
    xs = [randn(fill(2, length(ix))...) for ix in OMEinsum.getixs(tn)]
    simplifier, tn2 = merge_vectors(tn)
    @test tn2 == ein"ab,bc,cd,de,f->ab"
    @test tn(xs...) ≈ tn2(simplifier(xs...)...)

    tn = OMEinsum.DynamicEinCode(tn)
    xs = [randn(fill(2, length(ix))...) for ix in OMEinsum.getixs(tn)]
    simplifier, tn2 = merge_vectors(tn)
    @test tn2 == OMEinsum.DynamicEinCode(ein"ab,bc,cd,de,f->ab")
    @test tn(xs...) ≈ tn2(simplifier(xs...)...)
end