using OMEinsumContractionOrders
using Test, Random
using OMEinsumContractionOrders: merge_vectors, merge_greedy, optimize_greedy, embed_simplifier
using OMEinsum: decorate, rawcode, @ein_str
import OMEinsum

function apply_simplifier(s, xs)
    map(s.operations) do op
        return decorate(op)(xs...)
    end
end

function have_identity(ne::OMEinsumContractionOrders.NestedEinsum)
    if OMEinsumContractionOrders.isleaf(ne)
        return false
    elseif length(OMEinsumContractionOrders.getixsv(ne.eins)) == 1 && OMEinsumContractionOrders.getixsv(ne.eins)[1] == OMEinsumContractionOrders.getiyv(ne.eins)
        return true
    else
        return any(have_identity, ne.args)
    end
end

@testset "simplify vectors" begin
    tn = rawcode(ein"ab,bc,cd,de,a,b,c,d,f,f,e->ab")
    xs = [randn(fill(2, length(ix))...) for ix in OMEinsumContractionOrders.getixsv(tn)]
    simplifier, tn2 = merge_vectors(tn)
    @test tn2 == rawcode(ein"ab,bc,cd,de,f->ab")
    @test decorate(tn)(xs...) ≈ decorate(tn2)(apply_simplifier(simplifier, xs)...)
    tn3 = optimize_greedy(tn2, uniformsize(tn2, 2))
    tn4 = embed_simplifier(tn3, simplifier)
    @test !have_identity(tn4)
    @test decorate(tn)(xs...) ≈ decorate(tn4)(xs...)

    tn3 = optimize_greedy(tn2, uniformsize(tn2, 2))
    tn4 = embed_simplifier(tn3, simplifier)
    @test !have_identity(tn4)
    @test decorate(tn)(xs...) ≈ decorate(tn4)(xs...)
end

@testset "simplify greedy" begin
    tn = rawcode(ein"ab,bc,cd,de,a,b,c,d,f,f,e->ab")
    size_dict = uniformsize(tn, 2)
    xs = [randn(fill(2, length(ix))...) for ix in OMEinsumContractionOrders.getixsv(tn)]
    simplifier, tn2 = merge_greedy(tn, size_dict)
    @test tn2 == rawcode(ein"ab,->ab") || tn2 == rawcode(ein"ba,->ab")
    @test decorate(tn)(xs...) ≈ decorate(tn2)(apply_simplifier(simplifier, xs)...)
    tn3 = optimize_greedy(tn2, uniformsize(tn2, 2))
    tn4 = embed_simplifier(tn3, simplifier)
    @test !have_identity(tn4)
    @test decorate(tn)(xs...) ≈ decorate(tn4)(xs...)

    tn3 = optimize_greedy(tn2, uniformsize(tn2, 2))
    tn4 = embed_simplifier(tn3, simplifier)
    @test !have_identity(tn4)
    @test decorate(tn)(xs...) ≈ decorate(tn4)(xs...)
end

@testset "optimize permute" begin
    code = ein"lcij,lcjk->kcli" |> rawcode
    @test OMEinsumContractionOrders.optimize_output_permute(OMEinsumContractionOrders.getixsv(code), OMEinsumContractionOrders.getiyv(code)) == collect("iklc")
    code = ein"lcij,lcjk,c->kcli" |> rawcode
    @test OMEinsumContractionOrders.optimize_output_permute(OMEinsumContractionOrders.getixsv(code), OMEinsumContractionOrders.getiyv(code)) == collect("kcli")
    code = ein"lcij,lcjk->" |> rawcode
    @test OMEinsumContractionOrders.optimize_output_permute(OMEinsumContractionOrders.getixsv(code), OMEinsumContractionOrders.getiyv(code)) == Char[]
end

@testset "optimize permute" begin
    code = ein"(lcij,lcjk),lcik->" |> rawcode
    res = OMEinsumContractionOrders.NestedEinsum([OMEinsumContractionOrders.NestedEinsum([OMEinsumContractionOrders.NestedEinsum{Char}(1),
        OMEinsumContractionOrders.NestedEinsum{Char}(2)], ein"lcij,lcjk->iklc" |> rawcode), OMEinsumContractionOrders.NestedEinsum{Char}(3)], ein"iklc,lcik->" |> rawcode)
    @test optimize_permute(code) == res
    code = OMEinsumContractionOrders.SlicedEinsum(['c'], ein"(lcij,lcjk),lcik->" |> rawcode)
    @test optimize_permute(code) == OMEinsumContractionOrders.SlicedEinsum(['c'], res)
end