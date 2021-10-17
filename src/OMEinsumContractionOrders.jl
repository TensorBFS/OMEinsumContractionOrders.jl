module OMEinsumContractionOrders

using KaHyPar
using SparseArrays, Suppressor
using OMEinsum
using OMEinsum: NestedEinsum, getixs, getiy

include("kahypar.jl")
include("sa.jl")
include("treesa.jl")
include("simplify.jl")

end
