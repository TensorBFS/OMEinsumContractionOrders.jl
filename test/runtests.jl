using OMEinsumContractionOrders
using Test
using Documenter
using Aqua

@testset "Aqua" begin
    # Supressor does not pass the stale_deps test
    Aqua.test_all(OMEinsumContractionOrders; stale_deps = false)
end

@testset "Core" begin
    include("Core.jl")
end

@testset "utils" begin
    include("utils.jl")
end

@testset "greedy" begin
    include("greedy.jl")
end

@testset "sabipartite" begin
    include("sabipartite.jl")
end

if isdefined(Base, :get_extension)
    @testset "kahypar" begin
        include("kahypar.jl")
    end
end

@testset "treesa" begin
    include("treesa.jl")
end

@testset "treesaslicer" begin
    include("treesaslicer.jl")
end

@testset "treewidth" begin
    include("treewidth.jl")
end

@testset "hypernd" begin
    include("hypernd.jl")
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

DocMeta.setdocmeta!(OMEinsumContractionOrders, :DocTestSetup, :(using OMEinsumContractionOrders); recursive=true)
Documenter.doctest(OMEinsumContractionOrders; manual=false, fix=false)
