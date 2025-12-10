#!/usr/bin/env julia
"""
Julia service script for tensor network contraction order optimization.
This script reads JSON input from stdin and writes JSON output to stdout.

Input JSON format:
{
    "inputs": [[1, 2], [2, 3], [3, 4]],  // input tensor indices
    "output": [1, 4],                      // output tensor indices  
    "size_dict": {"1": 2, "2": 2, "3": 2, "4": 2},  // label -> size mapping
    "optimizer": "greedy",                 // optimizer type: "greedy", "treesa", "exacttreewidth"
    "optimizer_params": {}                 // optional optimizer parameters
}

Output JSON format:
{
    "success": true,
    "contraction_order": {...},           // nested einsum tree in JSON
    "complexity": {
        "time_complexity": 10.5,          // log2 of FLOPs
        "space_complexity": 5.0,          // log2 of max tensor size
        "read_write_complexity": 8.0      // log2 of read-write operations
    },
    "error": null
}
"""

using Pkg
# Activate the project environment
project_dir = dirname(@__DIR__)
Pkg.activate(project_dir)

using OMEinsumContractionOrders
using JSON

"""
    parse_label_type(inputs, output)

Determine the label type (Int or Char) from input/output specifications.
"""
function parse_label_type(inputs, output)
    # Check if labels are strings of length 1 (chars) or integers
    all_labels = vcat(inputs..., output)
    if all(l -> l isa String && length(l) == 1, all_labels)
        return Char
    elseif all(l -> l isa Integer || (l isa String && tryparse(Int, l) !== nothing), all_labels)
        return Int
    else
        return String
    end
end

"""
    convert_labels(inputs, output, LT)

Convert labels to the specified type.
"""
function convert_labels(inputs, output, ::Type{LT}) where LT
    if LT == Char
        conv = l -> l isa String ? only(l) : Char(l)
    elseif LT == Int
        conv = l -> l isa String ? parse(Int, l) : Int(l)
    else
        conv = l -> string(l)
    end
    return [[conv(l) for l in ix] for ix in inputs], [conv(l) for l in output]
end

"""
    parse_size_dict(size_dict_json, ::Type{LT}) where LT

Parse size dictionary with proper label types.
"""
function parse_size_dict(size_dict_json, ::Type{LT}) where LT
    if LT == Char
        conv = k -> k isa String ? (length(k) == 1 ? only(k) : error("Invalid char key: $k")) : Char(k)
    elseif LT == Int
        conv = k -> k isa String ? parse(Int, k) : Int(k)
    else
        conv = k -> string(k)
    end
    return Dict{LT, Int}(conv(k) => Int(v) for (k, v) in size_dict_json)
end

"""
    create_optimizer(optimizer_type, params)

Create an optimizer instance based on type string and parameters.
"""
function create_optimizer(optimizer_type::String, params::Dict)
    optimizer_type = lowercase(optimizer_type)
    
    if optimizer_type == "greedy"
        α = get(params, "alpha", get(params, "α", 0.0))
        temperature = get(params, "temperature", 0.0)
        return GreedyMethod(; α=α, temperature=temperature)
        
    elseif optimizer_type == "treesa"
        sc_target = get(params, "sc_target", 20.0)
        βs = get(params, "betas", get(params, "βs", 0.01:0.05:15.0))
        if βs isa Vector
            βs = βs[1]:βs[2]:βs[3]  # Convert [start, step, stop] to range
        end
        ntrials = get(params, "ntrials", 10)
        niters = get(params, "niters", 50)
        score = ScoreFunction(sc_target=sc_target)
        return TreeSA(; βs=βs, ntrials=ntrials, niters=niters, score=score)
        
    elseif optimizer_type == "exacttreewidth" || optimizer_type == "exact_treewidth"
        return ExactTreewidth()
        
    elseif optimizer_type == "treewidth"
        return Treewidth()
        
    elseif optimizer_type == "sabipartite"
        sc_target = get(params, "sc_target", 25)
        ntrials = get(params, "ntrials", 50)
        βs = get(params, "betas", get(params, "βs", 0.1:0.2:15.0))
        if βs isa Vector
            βs = βs[1]:βs[2]:βs[3]
        end
        niters = get(params, "niters", 1000)
        return SABipartite(; sc_target=sc_target, ntrials=ntrials, βs=βs, niters=niters)
        
    else
        error("Unknown optimizer type: $optimizer_type. Supported: greedy, treesa, exacttreewidth, treewidth, sabipartite")
    end
end

"""
    nested_einsum_to_dict(ne)

Convert a NestedEinsum or SlicedEinsum to a JSON-serializable dictionary.
"""
function nested_einsum_to_dict(ne)
    return OMEinsumContractionOrders._todict(ne)
end

"""
    optimize_contraction(request::Dict)

Main optimization function that processes a request dictionary.
"""
function optimize_contraction(request::Dict)
    try
        # Parse inputs
        inputs_raw = request["inputs"]
        output_raw = get(request, "output", [])
        size_dict_raw = request["size_dict"]
        optimizer_type = get(request, "optimizer", "greedy")
        optimizer_params = get(request, "optimizer_params", Dict())
        do_slicing = get(request, "slicing", false)
        slicing_params = get(request, "slicing_params", Dict())
        
        # Determine and convert label types
        LT = parse_label_type(inputs_raw, output_raw)
        inputs, output = convert_labels(inputs_raw, output_raw, LT)
        size_dict = parse_size_dict(size_dict_raw, LT)
        
        # Create einsum code
        code = OMEinsumContractionOrders.EinCode(inputs, output)
        
        # Create optimizer
        optimizer = create_optimizer(optimizer_type, optimizer_params)
        
        # Run optimization
        optcode = optimize_code(code, size_dict, optimizer)
        
        # Optionally apply slicing
        if do_slicing
            sc_target = get(slicing_params, "sc_target", 20.0)
            slicer = TreeSASlicer(score=ScoreFunction(sc_target=sc_target))
            optcode = slice_code(optcode, size_dict, slicer)
        end
        
        # Compute complexity
        complexity = contraction_complexity(optcode, size_dict)
        
        # Prepare response
        result = Dict(
            "success" => true,
            "contraction_order" => nested_einsum_to_dict(optcode),
            "complexity" => Dict(
                "time_complexity" => complexity.tc,
                "space_complexity" => complexity.sc,
                "read_write_complexity" => complexity.rwc
            ),
            "error" => nothing
        )
        
        return result
        
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
    to_dict(obj)

Recursively convert JSON objects to Dict.
"""
to_dict(obj::AbstractDict) = Dict{String, Any}(k => to_dict(v) for (k, v) in obj)
to_dict(obj::AbstractVector) = [to_dict(v) for v in obj]
to_dict(obj) = obj

"""
    process_request(json_str::String)

Process a JSON request string and return a JSON response string.
"""
function process_request(json_str::String)
    request = to_dict(JSON.parse(json_str))
    result = optimize_contraction(request)
    return JSON.json(result)
end

# Main entry point - read from stdin, write to stdout
function main()
    if !isinteractive()
        input = read(stdin, String)
        if !isempty(strip(input))
            output = process_request(input)
            println(output)
        end
    end
end

# Run main if this is the entry script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
