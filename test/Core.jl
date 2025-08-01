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

@testset "score" begin
    score = ScoreFunction(tc_weight=1.0, sc_weight=1.0, rw_weight=0.1, sc_target=20.0)
    @test score(10, 10, 10) ≈ 1024.0 + 102.4
    @test score(10, 30, 10) ≈ 1024.0 + exp2(30) - exp2(20) + 102.4
end

@testset "fix_binary_tree" begin
    # case 1: leaf and unary
    leaf1 = OMEinsumContractionOrders.NestedEinsum{Int}(2)
    leaf2 = OMEinsumContractionOrders.NestedEinsum{Int}(2)
    unary1 = OMEinsumContractionOrders.NestedEinsum([leaf1], OMEinsumContractionOrders.EinCode([[1, 2]], [2, 1]))
    unary2 = OMEinsumContractionOrders.NestedEinsum([leaf2], OMEinsumContractionOrders.EinCode([[1, 2]], [2, 1]))
    nested_unary = OMEinsumContractionOrders.NestedEinsum([unary1], OMEinsumContractionOrders.EinCode([[1, 2]], [1]))
    @test OMEinsumContractionOrders.fix_binary_tree(leaf1) == leaf1
    @test OMEinsumContractionOrders.fix_binary_tree(unary1) == unary1
    expected_fix1 = OMEinsumContractionOrders.NestedEinsum([leaf1], OMEinsumContractionOrders.EinCode([[1, 2]], [1]))
    @test OMEinsumContractionOrders.fix_binary_tree(nested_unary) == expected_fix1

    # case 2: binary tree
    binary1 = OMEinsumContractionOrders.NestedEinsum(
        [leaf1, leaf2],
        OMEinsumContractionOrders.EinCode([[1, 2], [1, 2]], [1])
    )
    @test OMEinsumContractionOrders.fix_binary_tree(binary1) == binary1
    @test OMEinsumContractionOrders.is_binary(binary1)

    binary2 = OMEinsumContractionOrders.NestedEinsum(
        [unary1, nested_unary],
        OMEinsumContractionOrders.EinCode([[1, 2], [1, 2]], [1])
    )
    expected_fix2 = OMEinsumContractionOrders.NestedEinsum([leaf1, leaf1], OMEinsumContractionOrders.EinCode([[1, 2], [1, 2]], [1]))
    @test OMEinsumContractionOrders.fix_binary_tree(binary2) == expected_fix2
    @test !OMEinsumContractionOrders.is_binary(binary2)
    @test OMEinsumContractionOrders.is_binary(expected_fix2)
    
    # case 3: unary with binary
    unary3 = OMEinsumContractionOrders.NestedEinsum([binary1], OMEinsumContractionOrders.EinCode([[1]], Int[]))
    expected_fix3 = OMEinsumContractionOrders.NestedEinsum(binary1.args, OMEinsumContractionOrders.EinCode([[1, 2], [1, 2]], Int[]))
    @test OMEinsumContractionOrders.fix_binary_tree(unary3) == expected_fix3
    @test !OMEinsumContractionOrders.is_binary(unary3)
    @test OMEinsumContractionOrders.is_binary(expected_fix3)

    unary4 = OMEinsumContractionOrders.NestedEinsum([binary2], OMEinsumContractionOrders.EinCode([[1]], Int[]))
    expected_fix4 = OMEinsumContractionOrders.NestedEinsum(expected_fix2.args, OMEinsumContractionOrders.EinCode([[1, 2], [1, 2]], Int[]))
    @test OMEinsumContractionOrders.fix_binary_tree(unary4) == expected_fix4
    @test !OMEinsumContractionOrders.is_binary(unary4)
    @test OMEinsumContractionOrders.is_binary(expected_fix4)
end