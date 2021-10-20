using OMEinsumContractionOrders
using Test

@testset "kahypar" begin
    include("kahypar.jl")
end

@testset "sa" begin
    include("sa.jl")
end

@testset "treesa" begin
    include("treesa.jl")
end

@testset "simplify" begin
    include("simplify.jl")
end

@testset "interfaces" begin
    include("interfaces.jl")
end
