module OMEinsumContractionOrders

using Requires
using SparseArrays, Suppressor
using OMEinsum
using OMEinsum: NestedEinsum, getixsv, getiyv, DynamicEinCode, StaticEinCode, isleaf, MinSpaceDiff, MinSpaceOut
export MinSpaceDiff, MinSpaceOut

using Requires
function __init__()
    @require KaHyPar="2a6221f6-aa48-11e9-3542-2d9e0ef01880" begin
        using .KaHyPar
        @info "`OMEinsumContractionOrders` loads `KaHyPar` module successfully."
    end
end

abstract type CodeOptimizer end

"""
    GreedyMethod{MT}
    GreedyMethod(; method=MinSpaceOut(), nrepeat=10)

The fast but poor greedy optimizer. Input arguments are

* `method` is `MinSpaceDiff()` or `MinSpaceOut`.
    * `MinSpaceOut` choose one of the contraction that produces a minimum output tensor size,
    * `MinSpaceDiff` choose one of the contraction that decrease the space most.
* `nrepeat` is the number of repeatition, returns the best contraction order.
"""
Base.@kwdef struct GreedyMethod{MT} <: CodeOptimizer
    method::MT = MinSpaceOut()
    nrepeat::Int = 10
end

include("kahypar.jl")
include("sa.jl")
include("treesa.jl")
include("simplify.jl")
include("interfaces.jl")

end
