using OMEinsumContractionOrders
using OMEinsumContractionOrders: EinCode
using CUDA
using cuTensorNet: CuTensorNetwork, cuTensorNet
using cuTensorNet: handle, cutensornetContractionOptimizerInfoGetAttribute, 
CUTENSORNET_CONTRACTION_OPTIMIZER_INFO_NUM_SLICED_MODES,
CUTENSORNET_CONTRACTION_OPTIMIZER_INFO_PATH, cutensornetContractionPath_t,
cutensornetNodePair_t

elty = Float32
bound_dim = 8
nmat = 5

ids = ['a'+i for i = 0:nmat]
code = EinCode(
    [[ids[i], ids[i+1]] for i = 1:nmat],
    [ids[1], ids[nmat+1]]
)
code.ixs[[1, 5]] = code.ixs[[5, 1]] 
code
d = Dict(ids[i] => nmat+3-i for i in 1:(nmat+1))

cte = optimize_code(code, d, CuTensorOptimizer(2^28); elty)
info = cte.info

arr = Vector{cuTensorNet.cutensornetNodePair_t}(undef, nmat - 1)
buf = Ref(cutensornetContractionPath_t(Int32(0), pointer(arr)))
cutensornetContractionOptimizerInfoGetAttribute(
    handle(), info, CUTENSORNET_CONTRACTION_OPTIMIZER_INFO_PATH,
    buf, sizeof(buf)
)
arr
flop(cte)
