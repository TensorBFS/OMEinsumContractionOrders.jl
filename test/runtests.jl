using OMEinsumContractionOrders
using Test

@testset "kahypar" begin
    include("kahypar.jl")
end

@testset "sa" begin
    include("sa.jl")
end
