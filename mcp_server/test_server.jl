#!/usr/bin/env julia
"""
Test the Julia MCP server locally
"""

include("server.jl")

println("="^60)
println("Testing Tensor Network Optimizer MCP Server")
println("="^60)

# Test 1: optimize_contraction_order
println("\n1. Testing optimize_contraction_order...")
params1 = Dict(
    "inputs" => [["i", "j"], ["j", "k"], ["k", "l"]],
    "output" => ["i", "l"],
    "size_dict" => Dict("i" => 100, "j" => 50, "k" => 200, "l" => 100),
    "optimizer" => "greedy"
)

result1 = handle_optimize(params1)
if result1["success"]
    println("✓ optimize_contraction_order works!")
    println("  Time complexity: 2^$(round(result1["complexity"]["time_complexity"], digits=2))")
    println("  Space complexity: 2^$(round(result1["complexity"]["space_complexity"], digits=2))")
else
    println("✗ Failed: $(result1["error"])")
end

# Test 2: analyze_einsum
println("\n2. Testing analyze_einsum...")
params2 = Dict(
    "einsum_string" => "ij,jk,kl->il",
    "sizes" => 10
)

result2 = handle_analyze_einsum(params2)
if result2["success"]
    println("✓ analyze_einsum works!")
    println("  Parsed: $(result2["einsum_parsed"]["inputs"]) -> $(result2["einsum_parsed"]["output"])")
    println("  Time complexity: 2^$(round(result2["complexity"]["time_complexity"], digits=2))")
else
    println("✗ Failed: $(result2["error"])")
end

# Test 3: compare_optimizers
println("\n3. Testing compare_optimizers...")
params3 = Dict(
    "inputs" => [[1, 2], [2, 3], [3, 4]],
    "output" => [1, 4],
    "size_dict" => Dict("1" => 10, "2" => 20, "3" => 30, "4" => 40),
    "optimizers" => ["greedy", "treewidth"]
)

result3 = handle_compare_optimizers(params3)
if result3["best_optimizer"] !== nothing
    println("✓ compare_optimizers works!")
    println("  Best optimizer: $(result3["best_optimizer"])")
    println("  Best time complexity: 2^$(round(result3["best_time_complexity"], digits=2))")
    for (opt, res) in result3["results"]
        if res["success"]
            tc = res["complexity"]["time_complexity"]
            println("    $opt: 2^$(round(tc, digits=2))")
        end
    end
else
    println("✗ Failed")
end

println("\n" * "="^60)
println("All tests completed!")
println("="^60)
