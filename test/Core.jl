using OMEinsumContractionOrders, OMEinsum
using OMEinsumContractionOrders: pivot_tree, path_to_tensor

using Test

@testset "tree reformulate" begin
    eincode = ein"((ik, jkl), ij), (lm, m) -> "
    code = OMEinsum.rawcode(eincode)

    size_dict = Dict([c=>(1<<i) for (i,c) in enumerate(['i', 'j', 'k', 'l', 'm'])]...)
    tensor_labels = [['i', 'k'], ['j', 'k', 'l'], ['i', 'j'], ['l', 'm'], ['m']]
    size_tensors = [log2(prod(size_dict[l] for l in tensor)) for tensor in tensor_labels]
    tensors = [rand([size_dict[j] for j in tensor_labels[i]]...) for i in 1:5]
    for tensor_index in 1:5

        path = path_to_tensor(code, tensor_index)
        tensor = reduce((x, y) -> x.args[y], path, init = code)
        @test tensor.tensorindex == tensor_index

        new_code = pivot_tree(code, tensor_index)
        @test contraction_complexity(new_code, size_dict).sc == max(contraction_complexity(code, size_dict).sc, size_tensors[tensor_index])

        closed_code = OMEinsumContractionOrders.NestedEinsum([new_code, tensor], OMEinsumContractionOrders.EinCode([OMEinsumContractionOrders.getiyv(new_code), tensor_labels[tensor_index]], Char[]))
        new_eincode = OMEinsum.decorate(closed_code)
        @test eincode(tensors...) ≈ new_eincode(tensors...)
    end
end

@testset "contraction complexity for empty tensor network" begin
    code = OMEinsumContractionOrders.EinCode(Vector{Char}[], Char[])
    @test contraction_complexity(code, Dict{Char, Int}()).sc == 0

    code = OMEinsumContractionOrders.NestedEinsum(OMEinsumContractionOrders.NestedEinsum{Char}[], OMEinsumContractionOrders.EinCode(Vector{Char}[], Char[]))
    @test contraction_complexity(code, Dict{Char, Int}()).sc == 0
end
