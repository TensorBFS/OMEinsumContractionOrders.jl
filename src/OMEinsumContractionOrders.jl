module OMEinsumContractionOrders

using Requires
using SparseArrays, Suppressor
using OMEinsum
using OMEinsum: NestedEinsum, getixsv, getiyv, DynamicEinCode, StaticEinCode, isleaf, MinSpaceDiff, MinSpaceOut
export MinSpaceDiff, MinSpaceOut
using JSON

using Requires
function __init__()
    @require KaHyPar="2a6221f6-aa48-11e9-3542-2d9e0ef01880" begin
        using .KaHyPar
        @info "`OMEinsumContractionOrders` loads `KaHyPar` module successfully."
    end
end

include("slicing.jl")
include("algorithms/algorithms.jl")
include("interfaces.jl")
include("json.jl")

end
