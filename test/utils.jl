using OMEinsumContractionOrders
using OMEinsumContractionOrders: log2sumexp2, _log2_size_dict, convert_label
using Test

@testset "log2sumexp2" begin
    @test log2sumexp2([1.0, 2.0, 3.0]) ≈ log2(exp2(1.0) + exp2(2.0) + exp2(3.0))
    @test log2sumexp2([10.0, 10.0]) ≈ 11.0
    @test log2sumexp2([0.0]) ≈ 0.0
end

@testset "_log2_size_dict" begin
    size_dict = Dict('a' => 2, 'b' => 4, 'c' => 8)
    log2_dict = _log2_size_dict(size_dict)
    @test log2_dict['a'] ≈ 1.0
    @test log2_dict['b'] ≈ 2.0
    @test log2_dict['c'] ≈ 3.0
end

@testset "convert_label" begin
    # Single leaf
    ne = OMEinsumContractionOrders.NestedEinsum{Char}(1)
    labelmap = Dict{Char, Int}()
    result = convert_label(ne, labelmap)
    @test result.tensorindex == 1
    
    # Simple contraction
    ne = OMEinsumContractionOrders.NestedEinsum([OMEinsumContractionOrders.NestedEinsum{Char}(1), OMEinsumContractionOrders.NestedEinsum{Char}(2)], 
                      OMEinsumContractionOrders.EinCode([['a','b'], ['b','c']], ['a','c']))
    labelmap = Dict('a' => 1, 'b' => 2, 'c' => 3)
    result = convert_label(ne, labelmap)
    @test result.eins.ixs == [[1,2], [2,3]]
    @test result.eins.iy == [1,3]
end

