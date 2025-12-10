#!/usr/bin/env python3
"""
MCP Server for Tensor Network Contraction Order Optimization.

This server provides tools for optimizing tensor network contraction orders
using the Julia package OMEinsumContractionOrders.
"""

import json
import subprocess
import sys
import os
from pathlib import Path
from typing import Any

from mcp.server.fastmcp import FastMCP

# Initialize MCP server
mcp = FastMCP("tensor-network-optimizer")

# Get the directory containing this script
SCRIPT_DIR = Path(__file__).parent.resolve()
JULIA_SERVICE = SCRIPT_DIR / "optimize_service.jl"
PROJECT_DIR = SCRIPT_DIR.parent


def run_julia_optimization(request: dict) -> dict:
    """
    Run the Julia optimization service with the given request.
    
    Args:
        request: Dictionary containing optimization parameters
        
    Returns:
        Dictionary with optimization results
    """
    try:
        # Prepare the Julia command
        julia_cmd = [
            "julia",
            "--project=" + str(PROJECT_DIR),
            str(JULIA_SERVICE)
        ]
        
        # Run Julia with the request as JSON input
        result = subprocess.run(
            julia_cmd,
            input=json.dumps(request),
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        if result.returncode != 0:
            return {
                "success": False,
                "error": f"Julia process failed: {result.stderr}",
                "contraction_order": None,
                "complexity": None
            }
        
        # Parse the output
        output = result.stdout.strip()
        if not output:
            return {
                "success": False,
                "error": "Empty output from Julia process",
                "contraction_order": None,
                "complexity": None
            }
            
        return json.loads(output)
        
    except subprocess.TimeoutExpired:
        return {
            "success": False,
            "error": "Optimization timed out after 5 minutes",
            "contraction_order": None,
            "complexity": None
        }
    except json.JSONDecodeError as e:
        return {
            "success": False,
            "error": f"Failed to parse Julia output: {e}",
            "contraction_order": None,
            "complexity": None
        }
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "contraction_order": None,
            "complexity": None
        }


@mcp.tool()
def optimize_contraction_order(
    inputs: list[list[Any]],
    size_dict: dict[str, int],
    output: list[Any] | None = None,
    optimizer: str = "greedy",
    optimizer_params: dict | None = None,
    slicing: bool = False,
    slicing_params: dict | None = None
) -> dict:
    """
    Optimize the contraction order of a tensor network.
    
    This tool finds an efficient order to contract tensors in a tensor network,
    minimizing time and space complexity.
    
    Args:
        inputs: List of tensor index lists. Each inner list contains the indices
                (labels) for one input tensor. Can use integers or single characters.
                Example: [[1, 2], [2, 3], [3, 4]] or [["i", "j"], ["j", "k"]]
        
        size_dict: Dictionary mapping index labels to their dimensions.
                   Example: {"1": 2, "2": 3, "3": 4} or {"i": 10, "j": 20}
        
        output: List of output tensor indices. Default is [] (scalar output).
                Example: [1, 4] means the result tensor has indices 1 and 4.
        
        optimizer: Optimization algorithm to use. Options:
                   - "greedy": Fast greedy algorithm (default)
                   - "treesa": Tree simulated annealing (better quality, slower)
                   - "exacttreewidth": Exact treewidth algorithm (optimal for small networks)
                   - "treewidth": Heuristic treewidth algorithm
                   - "sabipartite": Simulated annealing bipartite method
        
        optimizer_params: Optional parameters for the optimizer:
                         - For "greedy": {"alpha": 0.0, "temperature": 0.0}
                         - For "treesa": {"sc_target": 20, "ntrials": 10, "niters": 50}
                         - For "sabipartite": {"sc_target": 25, "ntrials": 50, "niters": 1000}
        
        slicing: Whether to apply slicing to reduce space complexity. Default False.
        
        slicing_params: Parameters for slicing: {"sc_target": 20.0}
    
    Returns:
        Dictionary containing:
        - success: Whether optimization succeeded
        - contraction_order: The optimized contraction tree in JSON format
        - complexity: Time/space/read-write complexity (log2 scale)
        - error: Error message if failed
    
    Example:
        # Matrix chain multiplication: A(i,j) * B(j,k) * C(k,l) -> D(i,l)
        result = optimize_contraction_order(
            inputs=[["i", "j"], ["j", "k"], ["k", "l"]],
            size_dict={"i": 100, "j": 50, "k": 200, "l": 100},
            output=["i", "l"],
            optimizer="treesa"
        )
    """
    request = {
        "inputs": inputs,
        "output": output if output is not None else [],
        "size_dict": size_dict,
        "optimizer": optimizer,
        "optimizer_params": optimizer_params or {},
        "slicing": slicing,
        "slicing_params": slicing_params or {}
    }
    
    return run_julia_optimization(request)


@mcp.tool()
def analyze_einsum(
    einsum_string: str,
    sizes: dict[str, int] | int = 2
) -> dict:
    """
    Analyze an einsum expression and compute optimal contraction complexity.
    
    Args:
        einsum_string: Einstein summation notation string.
                      Format: "ij,jk,kl->il" (like numpy.einsum)
                      
        sizes: Either a dictionary mapping indices to sizes,
               or a single integer for uniform sizes.
               Example: {"i": 10, "j": 20} or just 2 for all size 2
    
    Returns:
        Dictionary with optimization results and complexity analysis.
    
    Example:
        # Analyze matrix multiplication
        result = analyze_einsum("ij,jk->ik", {"i": 100, "j": 50, "k": 100})
        
        # Analyze with uniform size
        result = analyze_einsum("abc,bcd,cde->ae", 10)
    """
    try:
        # Parse einsum string
        if "->" in einsum_string:
            inputs_str, output_str = einsum_string.split("->")
        else:
            inputs_str = einsum_string
            output_str = ""
        
        # Parse input tensors (comma-separated)
        input_strs = inputs_str.split(",")
        inputs = [[c for c in s.strip()] for s in input_strs]
        output = [c for c in output_str.strip()]
        
        # Build size dictionary
        all_indices = set()
        for inp in inputs:
            all_indices.update(inp)
        all_indices.update(output)
        
        if isinstance(sizes, int):
            size_dict = {idx: sizes for idx in all_indices}
        else:
            size_dict = sizes
            # Check all indices have sizes
            missing = all_indices - set(size_dict.keys())
            if missing:
                return {
                    "success": False,
                    "error": f"Missing sizes for indices: {missing}",
                    "contraction_order": None,
                    "complexity": None
                }
        
        # Run optimization
        request = {
            "inputs": inputs,
            "output": output,
            "size_dict": size_dict,
            "optimizer": "treesa",
            "optimizer_params": {"ntrials": 5, "niters": 30}
        }
        
        result = run_julia_optimization(request)
        result["einsum_parsed"] = {
            "inputs": inputs,
            "output": output,
            "indices": list(all_indices)
        }
        
        return result
        
    except Exception as e:
        return {
            "success": False,
            "error": f"Failed to parse einsum string: {e}",
            "contraction_order": None,
            "complexity": None
        }


@mcp.tool()
def compare_optimizers(
    inputs: list[list[Any]],
    size_dict: dict[str, int],
    output: list[Any] | None = None,
    optimizers: list[str] | None = None
) -> dict:
    """
    Compare different optimization algorithms on the same tensor network.
    
    Args:
        inputs: List of tensor index lists (same as optimize_contraction_order)
        size_dict: Dictionary mapping index labels to dimensions
        output: Output tensor indices (default: [] for scalar)
        optimizers: List of optimizer names to compare.
                   Default: ["greedy", "treesa"]
    
    Returns:
        Dictionary with results from each optimizer for comparison.
    
    Example:
        result = compare_optimizers(
            inputs=[[1, 2], [2, 3], [3, 4], [4, 5]],
            size_dict={"1": 10, "2": 20, "3": 30, "4": 40, "5": 50},
            output=[1, 5],
            optimizers=["greedy", "treesa", "treewidth"]
        )
    """
    if optimizers is None:
        optimizers = ["greedy", "treesa"]
    
    results = {}
    
    for opt_name in optimizers:
        request = {
            "inputs": inputs,
            "output": output if output is not None else [],
            "size_dict": size_dict,
            "optimizer": opt_name,
            "optimizer_params": {}
        }
        
        result = run_julia_optimization(request)
        results[opt_name] = {
            "success": result.get("success", False),
            "complexity": result.get("complexity"),
            "error": result.get("error")
        }
    
    # Find best optimizer
    best_tc = float('inf')
    best_opt = None
    for opt_name, res in results.items():
        if res["success"] and res["complexity"]:
            tc = res["complexity"].get("time_complexity", float('inf'))
            if tc < best_tc:
                best_tc = tc
                best_opt = opt_name
    
    return {
        "results": results,
        "best_optimizer": best_opt,
        "best_time_complexity": best_tc if best_opt else None
    }


@mcp.resource("optimizer://docs")
def get_optimizer_docs() -> str:
    """Get documentation about available optimizers."""
    return """
# Tensor Network Contraction Order Optimizers

## Available Optimizers

### GreedyMethod (optimizer="greedy")
A fast greedy algorithm that iteratively contracts the pair of tensors with
minimum cost. Good for quick results but may not find optimal solutions.

Parameters:
- alpha (α): Controls the trade-off in the loss function (default: 0.0)
- temperature: For stochastic selection, 0 means deterministic (default: 0.0)

### TreeSA (optimizer="treesa")
Tree Simulated Annealing - uses simulated annealing to optimize the contraction
tree structure. Generally produces better results than greedy for complex networks.

Parameters:
- sc_target: Target space complexity (log2), optimization focuses on time after this
- ntrials: Number of independent trials (default: 10)
- niters: Iterations per trial (default: 50)
- betas: Temperature schedule [start, step, stop]

### ExactTreewidth (optimizer="exacttreewidth")
Finds the optimal contraction order using exact treewidth algorithms.
Only practical for small networks (< ~20 tensors).

### Treewidth (optimizer="treewidth") 
Heuristic treewidth-based optimizer using various algorithms from graph theory.
Good balance between quality and speed.

### SABipartite (optimizer="sabipartite")
Simulated annealing with bipartite recursive decomposition.
Good for networks with hierarchical structure.

Parameters:
- sc_target: Target space complexity (default: 25)
- ntrials: Number of trials (default: 50)
- niters: Iterations per trial (default: 1000)

## Complexity Metrics

Results include three complexity measures (all in log2 scale):
- time_complexity: Total number of floating point operations
- space_complexity: Size of the largest intermediate tensor
- read_write_complexity: Total memory read/write operations

## Example Usage

For a tensor network contraction like:
A[i,j] * B[j,k] * C[k,l] * D[l,m] -> E[i,m]

Use:
```
inputs = [["i","j"], ["j","k"], ["k","l"], ["l","m"]]
output = ["i", "m"]
size_dict = {"i": 100, "j": 50, "k": 200, "l": 100, "m": 150}
```
"""


def main():
    """Run the MCP server."""
    mcp.run()


if __name__ == "__main__":
    main()
