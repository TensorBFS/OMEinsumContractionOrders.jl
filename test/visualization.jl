using OMEinsum
using OMEinsumContractionOrders: ein2hypergraph, ein2elimination
using Test, OMEinsumContractionOrders

# tests before the extension loaded
@testset "luxor tensor plot dependency check" begin
    @test_throws ArgumentError begin
        eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], ['a'])
        ein2hypergraph(eincode)
    end

    @test_throws ArgumentError begin
        eincode = OMEinsum.rawcode(ein"((ij, jk), kl), lm -> im")
        ein2elimination(eincode)
    end

    @test_throws ArgumentError begin
        eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], Vector{Char}())
        viz_eins(eincode)
    end

    @test_throws ArgumentError begin
        eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], Vector{Char}())
        nested_code = optimize_code(eincode, uniformsize(eincode, 2), GreedyMethod())
        viz_contraction(nested_code, pathname = "")
    end
end

using LuxorGraphPlot
using LuxorGraphPlot.Luxor

@testset "eincode to hypergraph" begin
    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], ['a'])
    g1 = ein2hypergraph(eincode)

    nested_code = optimize_code(eincode, uniformsize(eincode, 2), GreedyMethod())
    g2 = ein2hypergraph(nested_code)

    sliced_code = optimize_code(eincode, uniformsize(eincode, 2), TreeSA(nslices = 1))
    g3 = ein2hypergraph(sliced_code)

    @test g1 == g2 == g3
    @test size(g1.adjacency_matrix, 1) == 5
    @test size(g1.adjacency_matrix, 2) == 6
end

@testset "eincode to elimination order" begin
    eincode = OMEinsum.rawcode(ein"((ij, jk), kl), lm -> im")
    elimination_order = ein2elimination(eincode)
    @test elimination_order == ['j', 'k', 'l']
end

@testset "visualize eincode" begin
    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], Vector{Char}())
    t = viz_eins(eincode)
    @test t isa Luxor.Drawing

    nested_code = optimize_code(eincode, uniformsize(eincode, 2), GreedyMethod())
    t = viz_eins(nested_code)
    @test t isa Luxor.Drawing

    sliced_code = optimize_code(eincode, uniformsize(eincode, 2), TreeSA())
    t = viz_eins(sliced_code)
    @test t isa Luxor.Drawing

    open_eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], ['a'])
    t = viz_eins(open_eincode)
    @test t isa Luxor.Drawing

    # filename and location specified
    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], Vector{Char}())
    filename = tempname() * ".png"
    viz_eins(eincode; filename, locs=vcat([(randn() * 60, 0.0) for i=1:5], [(randn() * 60, 320.0) for i=1:6]))
    @test isfile(filename)
end

@testset "visualize contraction" begin
    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], Vector{Char}())
    nested_code = optimize_code(eincode, uniformsize(eincode, 2), GreedyMethod())
    t_mp4 = viz_contraction(nested_code)
    tempmp4 = tempname() * ".mp4"
    tempgif = tempname() * ".gif"
    t_mp4_2 = viz_contraction(nested_code, filename = tempmp4)
    @test t_mp4 isa String
    @test t_mp4_2 isa String
    t_gif = viz_contraction(nested_code, filename = tempgif)
    @test t_gif isa String

    @test_throws AssertionError begin
        viz_contraction(nested_code, filename = "test.avi")
    end

    sliced_code = optimize_code(eincode, uniformsize(eincode, 2), TreeSA())
    t_mp4 = viz_contraction(sliced_code)
    t_mp4_2 = viz_contraction(sliced_code, filename = tempmp4)
    @test t_mp4 isa String
    @test t_mp4_2 isa String
    t_gif = viz_contraction(sliced_code, filename = tempgif)
    @test t_gif isa String

    sliced_code2 =  optimize_code(eincode, uniformsize(eincode, 2), TreeSA(nslices = 1))
    t_mp4 = viz_contraction(sliced_code2)
    t_mp4_2 = viz_contraction(sliced_code2, filename = tempmp4)
    @test t_mp4 isa String
    @test t_mp4_2 isa String
    t_gif = viz_contraction(sliced_code2, filename = tempgif)
    @test t_gif isa String
end