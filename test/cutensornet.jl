using Test
using CUDA
using OMEinsumContractionOrders
using OMEinsumContractionOrders: EinCode
using cuTensorNet: CuTensorNetwork

@testset "Fake types" begin
    elty = Float32
    felty = Int32
    bound_dim = 8
    code = EinCode(
        [['a', 'b'], ['b', 'c']],
        ['a', 'c']
    )
    d = Dict(i => bound_dim for i = ('a', 'b', 'c'))
    fxs = [CUDA.rand(felty, bound_dim, bound_dim), CUDA.rand(felty, bound_dim, bound_dim)]
    cte = optimize_code(code, d, CuTensorOptimizer(2^28); elty)
    @test collect(cte(fxs...)) != collect(fxs[1] * fxs[2])
end

@testset "Matrix MUL via TN (elty = $elty, dim = $bound_dim)" for elty in 
        (Float16, Float32, Float64), 
        bound_dim in [2^k for k in 2:5]
    code = EinCode(
        [['a', 'b'], ['b', 'c']],
        ['a', 'c']
    )
    d = Dict(i => bound_dim for i = ('a', 'b', 'c'))
    xs = [CUDA.rand(elty, bound_dim, bound_dim), CUDA.rand(elty, bound_dim, bound_dim)]
    cte = optimize_code(code, d, CuTensorOptimizer(2^28); elty)
    @test flop(code, d) == flop(cte)
    @test collect(cte(xs...)) ≈ collect(xs[1] * xs[2])
end
