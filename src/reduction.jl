struct Reduction{O <: CodeOptimizer} <: CodeOptimizer
    inner::O
end

struct ReductionEA{O <: CodeOptimizer} <: CliqueTrees.EliminationAlgorithm
    inner::O
end

#=
function CliqueTrees.permutation(weights::AbstractVector, graph, alg::ReductionEA)
    # reduce graph
    width = CliqueTrees.lowerbound(weights, graph)
    weights, graph, stack, index, width = CliqueTrees.saferules(weights, graph, width)
    
    # feed reduced graph back to OMEinsumContractionOrders
    code = EinCode(maximal_cliques(CliqueTrees.Graph(graph)), Int[])
    sizes = Dict{Int, Int}(v => round(Int, 2^weights[v]) for v in CliqueTrees.vertices(graph))
    opt = optimize_code(code, sizes, alg.inner).eins

    # compute ordering
    for v in eincode2order(opt)
        append!(stack, CliqueTrees.neighbors(index, v))
    end
    
    return stack, invperm(stack)
end
=#

function CliqueTrees.permutation(weights::AbstractVector, graph, alg::ReductionEA)
    # reduce graph
    kernel, stack, label, width = CliqueTrees.pr4(graph, CliqueTrees.lowerbound(graph))
    
    # feed reduced graph back to OMEinsumContractionOrders
    code = EinCode(maximal_cliques(kernel), Int[])
    sizes = Dict{Int, Int}(v => round(Int, 2^weights[label[v]]) for v in CliqueTrees.vertices(kernel))
    opt = optimize_code(code, sizes, alg.inner).eins
    
    # compute ordering
    append!(stack, label[eincode2order(opt)])
    return stack, invperm(stack)
end

function eincode2order(code::NestedEinsum{L}) where {L}
    elimination_order = Vector{L}()
    isleaf(code) && return elimination_order
    
    for node in PostOrderDFS(code)
        if !(node isa LeafString)
            for id in setdiff(vcat(getixsv(node.eins)...), getiyv(node.eins))
                push!(elimination_order, id)
            end
        end
    end
    
    return elimination_order
end
