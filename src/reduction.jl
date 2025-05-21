struct Reduction{S, O <: CodeOptimizer} <: CodeOptimizer
    inner::O
end

function Reduction{S}(inner::O) where {O <: CodeOptimizer}
    return Reduction{S, O}(inner)
end

function Reduction(inner::CodeOptimizer)
    return Reduction{3}(inner)
end

struct ReductionEA{S, O <: CodeOptimizer} <: CliqueTrees.EliminationAlgorithm
    inner::O
end

function ReductionEA{S}(inner::O) where {O <: CodeOptimizer}
    return ReductionEA{S, O}(inner)
end

# rules applied:
#  - islet
#  - twig
#  - series
#  - triangle
#  - buddy
#  - cube
function CliqueTrees.permutation(weights::AbstractVector, graph, alg::ReductionEA{1})
    # reduce graph
    width = CliqueTrees.lowerbound(weights, graph)
    graph, stack, label, width = CliqueTrees.pr3(weights, graph, width)
    
    # feed reduced graph back to OMEinsumContractionOrders
    code = EinCode(maximal_cliques(graph), Int[])
    sizes = Dict{Int, Int}(v => round(Int, 2^weights[label[v]]) for v in CliqueTrees.vertices(graph))
    opt = optimize_code(code, sizes, alg.inner).eins
    
    # compute ordering
    append!(stack, label[eincode2order(opt)])
    return stack, invperm(stack)
end

# rules applied:
#  - islet
#  - twig
#  - series
#  - triangle
#  - buddy
#  - cube
#  - simplicial
#  - almost-simplicial
function CliqueTrees.permutation(weights::AbstractVector, graph, alg::ReductionEA{2})
    # reduce graph
    width = CliqueTrees.lowerbound(weights, graph)
    graph, stack, label, width = CliqueTrees.pr4(graph, weights)
    
    # feed reduced graph back to OMEinsumContractionOrders
    code = EinCode(maximal_cliques(graph), Int[])
    sizes = Dict{Int, Int}(v => round(Int, 2^weights[label[v]]) for v in CliqueTrees.vertices(graph))
    opt = optimize_code(code, sizes, alg.inner).eins
    
    # compute ordering
    append!(stack, label[eincode2order(opt)])
    return stack, invperm(stack)
end

# rules applied:
#  - islet
#  - twig
#  - series
#  - triangle
#  - buddy
#  - cube
#  - simplicial
#  - almost-simplicial
#  - indistinguishable neighbors
function CliqueTrees.permutation(weights::AbstractVector, graph, alg::ReductionEA{3})
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
