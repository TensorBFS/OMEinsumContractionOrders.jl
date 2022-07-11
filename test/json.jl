using Test, OMEinsumContractionOrders

@testset "save load" begin
    for code in [
        OMEinsumContractionOrders.EinCode([[1,2], [2,3], [3,4]], [1,4]),
        OMEinsumContractionOrders.EinCode([['a','b'], ['b','c'], ['c','d']], ['a','d'])
    ]
        for optcode in [optimize_code(code, uniformsize(code, 2), GreedyMethod()),
            optimize_code(code, uniformsize(code, 2), TreeSA(nslices=1))]
            filename = tempname()
            writejson(filename, optcode)
            code2 = readjson(filename)
            @test optcode == code2
        end
    end
end