using Test, OMEinsumContractionOrders

@testset "save load" begin
    for code in [
        OMEinsumContractionOrders.EinCode([[1,2], [2,3], [3,4]], [1,4]),
        OMEinsumContractionOrders.EinCode([['a','b'], ['b','c'], ['c','d']], ['a','d'])
    ]
        code0 = optimize_code(code, uniformsize(code, 2), GreedyMethod())
        code1 = slice_code(code0, uniformsize(code, 2), TreeSASlicer(sc_target=2))
        for optcode in [code0, code1]
            filename = tempname()
            OMEinsumContractionOrders.writejson(filename, optcode)
            code2 = OMEinsumContractionOrders.readjson(filename)
            @test optcode == code2
        end
    end
end
