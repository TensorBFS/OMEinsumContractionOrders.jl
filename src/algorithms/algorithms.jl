module ContractionOrderAlgorithms

using SparseArrays
using Base: RefValue
using BetterExp
using Base.Threads
using Suppressor: @suppress

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