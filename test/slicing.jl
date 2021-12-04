using OMEinsum, OMEinsumContractionOrders
using Test, Random

@testset "slicer" begin
    log2_sizes = [1, 2,3 ,4.0]
    s = Slicer(log2_sizes, 3)
    push!(s, 1)
    @test_throws AssertionError push!(s, 1)
    push!(s, 2)
    push!(s, 3)
    @test_throws AssertionError push!(s, 4)
    replace!(s, 1=>4)
    @test s.log2_sizes == [1, 0.0, 0.0, 0.0]
    @test s.legs == Dict(2=>2.0, 3=>3.0, 4=>4.0)
    @test_throws AssertionError replace!(s, 1=>4)
end


@testset "SlicedEinsum" begin
    se = SlicedEinsum(Slicing(['i', 'l']), ein"(ij,jk),(kl,lm)->im")
    @test OMEinsum.flatten(se) == OMEinsum.flatten(se.eins)
    @test OMEinsum.labeltype(se) == Char
    xs = (randn(2,3), randn(3,4), randn(4,5), randn(5,6))
    size_info = Dict{Char,Int}()
    @test OMEinsum.get_size_dict!(se, xs, size_info) == Dict('i'=>2, 'j'=>3, 'k'=>4, 'l'=>5, 'm'=>6)
    @test OMEinsum.getixsv(se) == [['i','j'],['j','k'],['k','l'],['l','m']]
    @test OMEinsum.getiyv(se) == ['i','m']
    @test se(xs...) ≈ se.eins(xs...)
end