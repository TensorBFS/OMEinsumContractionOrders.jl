using LuxorGraphPlot
using LuxorGraphPlot.Luxor
using OMEinsumContractionOrders: ein2hypergraph, ein2elimination

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
    @test typeof(t) == Luxor.Drawing

    nested_code = optimize_code(eincode, uniformsize(eincode, 2), GreedyMethod())
    t = viz_eins(nested_code)
    @test typeof(t) == Luxor.Drawing

    sliced_code = optimize_code(eincode, uniformsize(eincode, 2), TreeSA())
    t = viz_eins(sliced_code)
    @test typeof(t) == Luxor.Drawing

    open_eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], ['a'])
    t = viz_eins(open_eincode)
    @test typeof(t) == Luxor.Drawing
end

@testset "visualize contraction" begin
    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], Vector{Char}())
    nested_code = optimize_code(eincode, uniformsize(eincode, 2), GreedyMethod())
    t_mp4 = viz_contraction(nested_code, pathname = "")
    @test typeof(t_mp4) == String
    t_gif = viz_contraction(nested_code, pathname = "", create_gif = true)
    @test typeof(t_gif) == String

    sliced_code = optimize_code(eincode, uniformsize(eincode, 2), TreeSA())
    t_mp4 = viz_contraction(sliced_code, pathname = "")
    @test typeof(t_mp4) == String
    t_gif = viz_contraction(sliced_code, pathname = "", create_gif = true)
    @test typeof(t_gif) == String

    sliced_code2 =  optimize_code(eincode, uniformsize(eincode, 2), TreeSA(nslices = 1))
    t_mp4 = viz_contraction(sliced_code2, pathname = "")
    @test typeof(t_mp4) == String
    t_gif = viz_contraction(sliced_code2, pathname = "", create_gif = true)
    @test typeof(t_gif) == String
end