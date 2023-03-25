using CUDA
using cuTensorNet
using cuTensorNet: OptimizerConfig, CuTensorNetwork, 
    CuTensorNetworkContractionOptimizerInfo, NoAutoTune, AutoTune,
    cutensornetGraphAlgo_t, cutensornetMemoryModel_t, cutensornetOptimizerCost_t,
    CUTENSORNET_GRAPH_ALGO_KWAY, CUTENSORNET_MEMORY_MODEL_CUTENSOR, CUTENSORNET_OPTIMIZER_COST_FLOPS,
    rehearse_contraction, perform_contraction!

struct CuTensorOptimizer <: CodeOptimizer
    config::OptimizerConfig
    sc_target::Int
end

CuTensorOptimizer(sc_target = 2^28; 
    num_graph_partitions::Int32 = Int32(8),
    graph_cutoff_size::Int32 = Int32(8),
    graph_algorithm::cutensornetGraphAlgo_t = CUTENSORNET_GRAPH_ALGO_KWAY,
    graph_imbalance_factor::Int32 = Int32(200),
    graph_num_iterations::Int32 = Int32(60),
    graph_num_cuts::Int32 = Int32(10),
    reconfig_num_iterations::Int32 = Int32(500),
    reconfig_num_leaves::Int32 = Int32(8),
    slicer_disable_slicing::Int32 = Int32(0),
    slicer_memory_model::cutensornetMemoryModel_t = CUTENSORNET_MEMORY_MODEL_CUTENSOR,
    slicer_memory_factor::Int32 = Int32(80),
    slicer_min_slices::Int32 = Int32(1),
    slicer_slice_factor::Int32 = Int32(2),
    hyper_num_samples::Int32 = Int32(0),
    simplification_disable_dr::Int32 = Int32(0),
    hyper_num_threads::Int32 = Int32(Threads.nthreads()),
    seed::Int32 = Int32(0),
    cost_function_objective::cutensornetOptimizerCost_t = CUTENSORNET_OPTIMIZER_COST_FLOPS,
) = CuTensorOptimizer(
    OptimizerConfig(;
        num_graph_partitions,
        graph_cutoff_size,
        graph_algorithm,
        graph_imbalance_factor,
        graph_num_iterations,
        graph_num_cuts,
        reconfig_num_iterations,
        reconfig_num_leaves,
        slicer_disable_slicing,
        slicer_memory_model,
        slicer_memory_factor,
        slicer_min_slices,
        slicer_slice_factor,
        hyper_num_samples,
        simplification_disable_dr,
        hyper_num_threads,
        seed,
        cost_function_objective
    ),
    sc_target
)

struct CuTensorEinsum <: AbstractEinsum
    info::CuTensorNetworkContractionOptimizerInfo
    ctn::CuTensorNetwork
end

function cuTensorNet.CuTensorNetwork(elty::DataType, code::EinCode{LT}, size_dict::Dict{LT, <:Integer};
        input_strides = [C_NULL for _ = 1:length(code.ixs)], 
        input_qualifiers = Int32[0 for _ = 1:length(code.ixs)],
        output_strides = C_NULL) where LT
    ixs = code.ixs
    iy = code.iy
    idx_map = Dict{LT, Int32}()
    curr_id = one(Int32)
    for ix in ixs
        for i in ix
            !haskey(idx_map, i) && (idx_map[i] = curr_id; curr_id += 1)
        end
    end
    for i in iy
        !haskey(idx_map, i) && (idx_map[i] = curr_id; curr_id += 1)
    end
    code_mapped = EinCode{Int32}([[idx_map[i] for i in ix] for ix in ixs], [idx_map[i] for i in iy])
    size_dict_mapped = Dict{Int32, Int}(idx_map[id] => v for (id, v) in size_dict)
    return CuTensorNetwork(elty, code_mapped, size_dict_mapped; 
        input_strides, input_qualifiers, output_strides)
end

function cuTensorNet.CuTensorNetwork(elty::DataType, code::EinCode{<:Integer}, size_dict::Dict{<:Integer, <:Integer}; 
        input_strides = [C_NULL for _ = 1:length(code.ixs)], 
        input_qualifiers = Int32[0 for _ = 1:length(code.ixs)],
        output_strides = C_NULL)
    input_modes = [Int32.(ix) for ix in code.ixs]
    input_extents = [[Int(size_dict[i]) for i in ix] for ix in input_modes]
    output_modes = Int32.(code.iy)
    output_extents = [Int(size_dict[i]) for i in output_modes]
    return cuTensorNet.CuTensorNetwork(elty, 
        input_modes, input_extents, input_strides, input_qualifiers, 
        output_modes, output_extents, output_strides)
end

function optimize_code(code::EinCode, size_dict::Dict, optimizer::CuTensorOptimizer = CuTensorOptimizer(); elty::DataType = Float32)
    return optimize_cutensornet(code, size_dict, optimizer; elty)
end

function optimize_cutensornet(code::EinCode, size_dict::Dict, optimizer::CuTensorOptimizer = CuTensorOptimizer(); elty::DataType = Float32)
    ctn = CuTensorNetwork(elty, code, size_dict)
    config = optimizer.config
    sc_target = optimizer.sc_target
    info = rehearse_contraction(ctn, sc_target, config)
    return CuTensorEinsum(info, ctn)
end

tensor_eltype(::CuTensorNetwork{T}) where T = T

function (code::CuTensorEinsum)(xs::CuArray...; tuning::Union{NoAutoTune, AutoTune} = AutoTune(), kwargs...)
    ctn = code.ctn
    info = code.info
    elty = tensor_eltype(ctn)
    if !all(eltype(x) === elty for x in xs)
        @warn "Applying on fake types"
        elty = eltype(xs[1])
        ctn = CuTensorNetwork{elty}(
            ctn.desc,
            ctn.input_modes,
            ctn.input_extents,
            ctn.input_strides,
            ctn.input_qualifiers,
            collect(xs),
            ctn.output_modes,
            ctn.output_extents,
            ctn.output_strides,
            CUDA.zeros(elty, ctn.output_extents...)
        )
    else
        ctn.input_arrs = collect(xs)
        ctn.output_arr = CUDA.zeros(elty, ctn.output_extents...)
    end
    perform_contraction!(ctn, info, tuning; kwargs...)
    return ctn.output_arr
end
