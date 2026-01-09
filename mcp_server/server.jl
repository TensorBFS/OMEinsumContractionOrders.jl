#!/usr/bin/env julia
"""
Pure Julia MCP Server for Tensor Network Contraction Order Optimization.

This server provides tools for optimizing tensor network contraction orders
using OMEinsumContractionOrders.jl directly, with no subprocess overhead.

Usage:
    julia --project=. server.jl

For Cursor MCP configuration (~/.cursor/mcp.json):
{
  "mcpServers": {
    "tensor-network-optimizer": {
      "command": "julia",
      "args": ["--project=/path/to/mcp_server", "/path/to/mcp_server/server.jl"]
    }
  }
}
"""

using ModelContextProtocol
using OMEinsumContractionOrders
using JSON3

# ============================================================================
# Helper Functions
# ============================================================================

"""
Parse label type from inputs - determines if labels are Char or Int.
"""
function parse_label_type(inputs::Vector, output::Vector)
    all_labels = vcat(inputs..., output)
    isempty(all_labels) && return Int
    first_label = first(all_labels)
    if first_label isa String && length(first_label) == 1
        return Char
    else
        return Int
    end
end

"""
Convert labels to the target type.
"""
function convert_labels(inputs::Vector, output::Vector, ::Type{LT}) where LT
    conv = if LT == Char
        l -> l isa String ? only(l) : Char(l)
    elseif LT == Int
        l -> l isa String ? parse(Int, l) : Int(l)
    else
        l -> string(l)
    end
    converted_inputs = [[conv(l) for l in ix] for ix in inputs]
    converted_output = [conv(l) for l in output]
    return converted_inputs, converted_output
end

"""
Parse size dictionary with proper label types.
"""
function parse_size_dict(size_dict::Dict, ::Type{LT}) where LT
    conv = if LT == Char
        k -> k isa String && length(k) == 1 ? only(k) : Char(k)
    elseif LT == Int
        k -> k isa String ? parse(Int, k) : Int(k)
    else
        k -> string(k)
    end
    return Dict{LT,Int}(conv(k) => Int(v) for (k, v) in size_dict)
end

"""
Create optimizer from type string and parameters.
"""
function create_optimizer(optimizer_type::String, params::Dict)
    opt = lowercase(optimizer_type)
    
    if opt == "greedy"
        α = get(params, "alpha", get(params, "α", 0.0))
        temperature = get(params, "temperature", 0.0)
        return GreedyMethod(; α=α, temperature=temperature)
        
    elseif opt == "treesa"
        sc_target = get(params, "sc_target", 20.0)
        ntrials = get(params, "ntrials", 10)
        niters = get(params, "niters", 50)
        score = ScoreFunction(sc_target=sc_target)
        return TreeSA(; ntrials=ntrials, niters=niters, score=score)
        
    elseif opt == "exacttreewidth" || opt == "exact_treewidth"
        return ExactTreewidth()
        
    elseif opt == "treewidth"
        return Treewidth()
        
    elseif opt == "sabipartite"
        sc_target = get(params, "sc_target", 25)
        ntrials = get(params, "ntrials", 50)
        niters = get(params, "niters", 1000)
        return SABipartite(; sc_target=sc_target, ntrials=ntrials, niters=niters)
        
    else
        error("Unknown optimizer: $optimizer_type. Options: greedy, treesa, treewidth, exacttreewidth, sabipartite")
    end
end

"""
Convert NestedEinsum/SlicedEinsum to dictionary using built-in JSON support.
"""
to_result_dict(ne) = OMEinsumContractionOrders._todict(ne)

# ============================================================================
# MCP Tool Handlers
# ============================================================================

"""
Handler for optimize_contraction_order tool.
"""
function handle_optimize(params::Dict)
    try
        # Extract parameters
        inputs_raw = params["inputs"]
        output_raw = get(params, "output", [])
        size_dict_raw = params["size_dict"]
        optimizer_type = get(params, "optimizer", "greedy")
        optimizer_params = get(params, "optimizer_params", Dict())
        do_slicing = get(params, "slicing", false)
        slicing_params = get(params, "slicing_params", Dict())
        
        # Convert types
        LT = parse_label_type(inputs_raw, output_raw)
        inputs, output = convert_labels(inputs_raw, output_raw, LT)
        size_dict = parse_size_dict(size_dict_raw, LT)
        
        # Create einsum code and optimizer
        code = OMEinsumContractionOrders.EinCode(inputs, output)
        optimizer = create_optimizer(optimizer_type, optimizer_params)
        
        # Run optimization
        optcode = optimize_code(code, size_dict, optimizer)
        
        # Apply slicing if requested
        if do_slicing
            sc_target = get(slicing_params, "sc_target", 20.0)
            slicer = TreeSASlicer(score=ScoreFunction(sc_target=sc_target))
            optcode = slice_code(optcode, size_dict, slicer)
        end
        
        # Compute complexity
        complexity = contraction_complexity(optcode, size_dict)
        
        # Return native Dict - MCP handles JSON serialization
        return Dict(
            "success" => true,
            "contraction_order" => to_result_dict(optcode),
            "complexity" => Dict(
                "time_complexity" => complexity.tc,
                "space_complexity" => complexity.sc,
                "read_write_complexity" => complexity.rwc
            ),
            "error" => nothing
        )
        
    catch e
        return Dict(
            "success" => false,
            "contraction_order" => nothing,
            "complexity" => nothing,
            "error" => sprint(showerror, e, catch_backtrace())
        )
    end
end

"""
Handler for analyze_einsum tool.
"""
function handle_analyze_einsum(params::Dict)
    try
        einsum_string = params["einsum_string"]
        sizes = get(params, "sizes", 2)
        
        # Parse einsum string
        if occursin("->", einsum_string)
            parts = split(einsum_string, "->")
            inputs_str = strip(parts[1])
            output_str = strip(parts[2])
        else
            inputs_str = einsum_string
            output_str = ""
        end
        
        # Parse input tensors
        input_strs = split(inputs_str, ",")
        inputs = [[c for c in strip(s)] for s in input_strs]
        output = [c for c in output_str]
        
        # Build size dictionary
        all_indices = Set{Char}()
        for inp in inputs
            union!(all_indices, inp)
        end
        union!(all_indices, output)
        
        size_dict = if sizes isa Integer
            Dict{Char,Int}(idx => sizes for idx in all_indices)
        else
            Dict{Char,Int}(only(string(k)) => Int(v) for (k, v) in sizes)
        end
        
        # Run optimization with TreeSA
        code = OMEinsumContractionOrders.EinCode(inputs, output)
        optimizer = TreeSA(ntrials=5, niters=30)
        optcode = optimize_code(code, size_dict, optimizer)
        complexity = contraction_complexity(optcode, size_dict)
        
        # Return native Dict - MCP handles JSON serialization
        return Dict(
            "success" => true,
            "contraction_order" => to_result_dict(optcode),
            "complexity" => Dict(
                "time_complexity" => complexity.tc,
                "space_complexity" => complexity.sc,
                "read_write_complexity" => complexity.rwc
            ),
            "einsum_parsed" => Dict(
                "inputs" => inputs,
                "output" => output,
                "indices" => collect(all_indices)
            ),
            "error" => nothing
        )
        
    catch e
        return Dict(
            "success" => false,
            "error" => sprint(showerror, e, catch_backtrace())
        )
    end
end

"""
Handler for compare_optimizers tool.
"""
function handle_compare_optimizers(params::Dict)
    try
        inputs_raw = params["inputs"]
        output_raw = get(params, "output", [])
        size_dict_raw = params["size_dict"]
        optimizers = get(params, "optimizers", ["greedy", "treesa"])
        
        # Convert types
        LT = parse_label_type(inputs_raw, output_raw)
        inputs, output = convert_labels(inputs_raw, output_raw, LT)
        size_dict = parse_size_dict(size_dict_raw, LT)
        code = OMEinsumContractionOrders.EinCode(inputs, output)
        
        results = Dict{String,Any}()
        best_tc = Inf
        best_opt = nothing
        
        for opt_name in optimizers
            try
                optimizer = create_optimizer(opt_name, Dict())
                optcode = optimize_code(code, size_dict, optimizer)
                complexity = contraction_complexity(optcode, size_dict)
                
                results[opt_name] = Dict(
                    "success" => true,
                    "complexity" => Dict(
                        "time_complexity" => complexity.tc,
                        "space_complexity" => complexity.sc,
                        "read_write_complexity" => complexity.rwc
                    ),
                    "error" => nothing
                )
                
                if complexity.tc < best_tc
                    best_tc = complexity.tc
                    best_opt = opt_name
                end
            catch e
                results[opt_name] = Dict(
                    "success" => false,
                    "complexity" => nothing,
                    "error" => sprint(showerror, e)
                )
            end
        end
        
        # Return native Dict - MCP handles JSON serialization
        return Dict(
            "results" => results,
            "best_optimizer" => best_opt,
            "best_time_complexity" => best_tc == Inf ? nothing : best_tc
        )
        
    catch e
        return Dict(
            "success" => false,
            "error" => sprint(showerror, e, catch_backtrace())
        )
    end
end

# ============================================================================
# MCP Server Definition
# ============================================================================

server = mcp_server(
    name = "tensor-network-optimizer",
    version = "1.0.0",
    tools = [
        MCPTool(
            name = "optimize_contraction_order",
            description = """Optimize the contraction order of a tensor network.
            
Finds an efficient order to contract tensors, minimizing time and space complexity.

Parameters:
- inputs: List of tensor index lists, e.g., [["i","j"], ["j","k"], ["k","l"]]
- size_dict: Dictionary mapping indices to sizes, e.g., {"i": 100, "j": 50}
- output: Output tensor indices (default: [] for scalar)
- optimizer: Algorithm - "greedy" (default), "treesa", "treewidth", "exacttreewidth", "sabipartite"
- optimizer_params: Optional parameters for the optimizer
- slicing: Enable slicing for memory reduction (default: false)

Returns JSON with: success, contraction_order, complexity (tc/sc/rwc in log2 scale), error""",
            parameters = [
                ToolParameter(name="inputs", type="array", description="List of tensor index lists", required=true),
                ToolParameter(name="size_dict", type="object", description="Index to size mapping", required=true),
                ToolParameter(name="output", type="array", description="Output tensor indices", required=false),
                ToolParameter(name="optimizer", type="string", description="Optimizer: greedy, treesa, treewidth, exacttreewidth, sabipartite", required=false),
                ToolParameter(name="optimizer_params", type="object", description="Optimizer parameters", required=false),
                ToolParameter(name="slicing", type="boolean", description="Enable slicing", required=false),
                ToolParameter(name="slicing_params", type="object", description="Slicing parameters", required=false)
            ],
            handler = handle_optimize
        ),
        MCPTool(
            name = "analyze_einsum",
            description = """Parse and optimize an einsum expression.

Parameters:
- einsum_string: Einstein notation, e.g., "ij,jk,kl->il"
- sizes: Size dict {"i": 10, "j": 20} or uniform integer size

Returns JSON with optimized contraction order and complexity.""",
            parameters = [
                ToolParameter(name="einsum_string", type="string", description="Einsum notation like 'ij,jk->ik'", required=true),
                ToolParameter(name="sizes", type="object", description="Size dict or uniform integer", required=false)
            ],
            handler = handle_analyze_einsum
        ),
        MCPTool(
            name = "compare_optimizers",
            description = """Compare different optimization algorithms on the same tensor network.

Parameters:
- inputs, size_dict, output: Same as optimize_contraction_order
- optimizers: List of optimizer names (default: ["greedy", "treesa"])

Returns comparison results with best optimizer identified.""",
            parameters = [
                ToolParameter(name="inputs", type="array", description="List of tensor index lists", required=true),
                ToolParameter(name="size_dict", type="object", description="Index to size mapping", required=true),
                ToolParameter(name="output", type="array", description="Output tensor indices", required=false),
                ToolParameter(name="optimizers", type="array", description="List of optimizer names to compare", required=false)
            ],
            handler = handle_compare_optimizers
        )
    ],
    resources = [
        MCPResource(
            uri = "optimizer://docs",
            name = "Optimizer Documentation",
            description = "Documentation about available optimizers",
            mime_type = "text/markdown",
            data_provider = function()
                """
# Tensor Network Contraction Order Optimizers

## Available Optimizers

### GreedyMethod (optimizer="greedy")
Fast greedy algorithm. Parameters: alpha (α), temperature.

### TreeSA (optimizer="treesa")  
Tree Simulated Annealing. Better quality for complex networks.
Parameters: sc_target, ntrials, niters.

### Treewidth (optimizer="treewidth")
Heuristic treewidth-based optimizer.

### ExactTreewidth (optimizer="exacttreewidth")
Optimal for small networks (<20 tensors).

### SABipartite (optimizer="sabipartite")
SA with bipartite decomposition.
Parameters: sc_target, ntrials, niters.

## Complexity Metrics (log2 scale)
- time_complexity: Total FLOPs
- space_complexity: Largest intermediate tensor
- read_write_complexity: Memory operations
"""
            end
        )
    ]
)

# ============================================================================
# Main Entry Point
# ============================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    start!(server)
end
