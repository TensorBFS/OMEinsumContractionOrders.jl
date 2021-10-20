using OMEinsumContractionOrders, OMEinsum
using Test, Random
using OMEinsum: getixs, getiy, isleaf

function have_identity(ne::NestedEinsum)
    if isleaf(ne)
        return false
    elseif length(getixs(ne.eins)) == 1 && getixs(ne.eins)[1] == getiy(ne.eins)
        return true
    else
        return any(have_identity, ne.args)
    end
end

@testset "simplify vectors" begin
    tn = ein"ab,bc,cd,de,a,b,c,d,f,f,e->ab"
    xs = [randn(fill(2, length(ix))...) for ix in OMEinsum.getixs(tn)]
    simplifier, tn2 = merge_vectors(tn)
    @test tn2 == ein"ab,bc,cd,de,f->ab"
    @test tn(xs...) ≈ tn2(simplifier(xs...)...)
    tn3 = optimize_greedy(tn2, uniformsize(tn2, 2))
    tn4 = embed_simplifier(tn3, simplifier)
    @test !have_identity(tn4)
    @test tn(xs...) ≈ tn4(xs...)

    tn = OMEinsum.DynamicEinCode(tn)
    xs = [randn(fill(2, length(ix))...) for ix in OMEinsum.getixs(tn)]
    simplifier, tn2 = merge_vectors(tn)
    @test tn2 == OMEinsum.DynamicEinCode(ein"ab,bc,cd,de,f->ab")
    @test tn(xs...) ≈ tn2(simplifier(xs...)...)
    tn3 = optimize_greedy(tn2, uniformsize(tn2, 2))
    tn4 = embed_simplifier(tn3, simplifier)
    @test !have_identity(tn4)
    @test tn(xs...) ≈ tn4(xs...)
end

@testset "simplify greedy" begin
    tn = ein"ab,bc,cd,de,a,b,c,d,f,f,e->ab"
    size_dict = uniformsize(tn, 2)
    xs = [randn(fill(2, length(ix))...) for ix in OMEinsum.getixs(tn)]
    simplifier, tn2 = merge_greedy(tn, size_dict)
    @test tn2 == ein"ab,->ab"
    @test tn(xs...) ≈ tn2(simplifier(xs...)...)
    tn3 = optimize_greedy(tn2, uniformsize(tn2, 2))
    tn4 = embed_simplifier(tn3, simplifier)
    @test !have_identity(tn4)
    @test tn(xs...) ≈ tn4(xs...)

    tn = OMEinsum.DynamicEinCode(tn)
    xs = [randn(fill(2, length(ix))...) for ix in OMEinsum.getixs(tn)]
    simplifier, tn2 = merge_greedy(tn, size_dict)
    @test tn2 == OMEinsum.DynamicEinCode(ein"ab,->ab")
    @test tn(xs...) ≈ tn2(simplifier(xs...)...)
    tn3 = optimize_greedy(tn2, uniformsize(tn2, 2))
    tn4 = embed_simplifier(tn3, simplifier)
    @test !have_identity(tn4)
    @test tn(xs...) ≈ tn4(xs...)
end
