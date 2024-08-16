using OMEinsumContractionOrders
using Test

@testset "Core" begin
    include("Core.jl")
end

@testset "greedy" begin
    include("greedy.jl")
end

@testset "sa" begin
    include("sa.jl")
end

if isdefined(Base, :get_extension)
    @testset "kahypar" begin
        include("kahypar.jl")
    end
end

@testset "treesa" begin
    include("treesa.jl")
end

@testset "treewidth" begin
    include("treewidth.jl")
end

@testset "simplify" begin
    include("simplify.jl")
end

@testset "interfaces" begin
    include("interfaces.jl")
end

@testset "json" begin
    include("json.jl")
end

# testing the extension `LuxorTensorPlot` for visualization
if isdefined(Base, :get_extension)
    @testset "visualization" begin
        include("visualization.jl")
    end
end