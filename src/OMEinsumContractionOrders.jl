module OMEinsumContractionOrders

using Requires
using SparseArrays, Suppressor
using OMEinsum
using OMEinsum: NestedEinsum, getixsv, getiyv, DynamicEinCode, StaticEinCode, isleaf

using Requires
function __init__()
    @require KaHyPar="2a6221f6-aa48-11e9-3542-2d9e0ef01880" begin
        using .KaHyPar
        @info "`OMEinsumContractionOrders` loads `KaHyPar` module successfully."
    end
end

include("kahypar.jl")
include("sa.jl")
include("treesa.jl")
include("simplify.jl")

end
