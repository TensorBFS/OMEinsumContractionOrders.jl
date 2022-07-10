module ContractionOrderAlgorithms

using SparseArrays
using Base: RefValue
using BetterExp
using Base.Threads
using Suppressor: @suppress

using Requires
function __init__()
    @require KaHyPar="2a6221f6-aa48-11e9-3542-2d9e0ef01880" begin
        using .KaHyPar
        @info "`OMEinsumContractionOrders` loads `KaHyPar` module successfully."
    end
end

include("Core.jl")
include("utils.jl")

# greedy method
include("incidencelist.jl")
include("greedy.jl")

# bipartition based methods
include("sa.jl")
include("kahypar.jl")

# local search method
include("treesa.jl")

# simplification passes
include("simplify.jl")

# interfaces
include("complexity.jl")
include("interfaces.jl")

end