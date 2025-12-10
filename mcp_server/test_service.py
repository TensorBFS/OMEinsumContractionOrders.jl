#!/usr/bin/env python3
"""Test script for the tensor network optimizer MCP server."""

import json
import subprocess
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.resolve()
JULIA_SERVICE = SCRIPT_DIR / "optimize_service.jl"
PROJECT_DIR = SCRIPT_DIR.parent


def test_julia_service():
    """Test the Julia optimization service directly."""
    print("Testing Julia optimization service...")
    
    # Test case: Simple matrix chain multiplication
    request = {
        "inputs": [["i", "j"], ["j", "k"], ["k", "l"]],
        "output": ["i", "l"],
        "size_dict": {"i": 10, "j": 20, "k": 30, "l": 40},
        "optimizer": "greedy"
    }
    
    julia_cmd = [
        "julia",
        "--project=" + str(PROJECT_DIR),
        str(JULIA_SERVICE)
    ]
    
    print(f"Running: {' '.join(julia_cmd)}")
    print(f"Input: {json.dumps(request, indent=2)}")
    
    result = subprocess.run(
        julia_cmd,
        input=json.dumps(request),
        capture_output=True,
        text=True,
        timeout=120
    )
    
    print(f"\nReturn code: {result.returncode}")
    
    if result.stderr:
        print(f"\nStderr:\n{result.stderr}")
    
    if result.stdout:
        print(f"\nStdout:\n{result.stdout}")
        try:
            output = json.loads(result.stdout.strip())
            print(f"\nParsed output:\n{json.dumps(output, indent=2)}")
        except json.JSONDecodeError as e:
            print(f"\nFailed to parse JSON: {e}")
    
    return result.returncode == 0


def test_mcp_server():
    """Test the MCP server tools (requires mcp package)."""
    try:
        from server import optimize_contraction_order, analyze_einsum, compare_optimizers
        
        print("\n" + "="*60)
        print("Testing MCP server tools...")
        print("="*60)
        
        # Test 1: optimize_contraction_order
        print("\n1. Testing optimize_contraction_order...")
        result = optimize_contraction_order(
            inputs=[["i", "j"], ["j", "k"], ["k", "l"]],
            size_dict={"i": 10, "j": 20, "k": 30, "l": 40},
            output=["i", "l"],
            optimizer="greedy"
        )
        print(f"Success: {result.get('success')}")
        if result.get('complexity'):
            print(f"Time complexity: 2^{result['complexity']['time_complexity']:.2f}")
            print(f"Space complexity: 2^{result['complexity']['space_complexity']:.2f}")
        if result.get('error'):
            print(f"Error: {result['error']}")
        
        # Test 2: analyze_einsum
        print("\n2. Testing analyze_einsum...")
        result = analyze_einsum("ij,jk,kl->il", {"i": 10, "j": 20, "k": 30, "l": 40})
        print(f"Success: {result.get('success')}")
        if result.get('complexity'):
            print(f"Time complexity: 2^{result['complexity']['time_complexity']:.2f}")
        if result.get('error'):
            print(f"Error: {result['error']}")
        
        # Test 3: compare_optimizers
        print("\n3. Testing compare_optimizers...")
        result = compare_optimizers(
            inputs=[[1, 2], [2, 3], [3, 4]],
            size_dict={"1": 10, "2": 20, "3": 30, "4": 40},
            output=[1, 4],
            optimizers=["greedy"]
        )
        print(f"Best optimizer: {result.get('best_optimizer')}")
        for opt, res in result.get('results', {}).items():
            print(f"  {opt}: TC=2^{res['complexity']['time_complexity']:.2f}" if res.get('complexity') else f"  {opt}: failed")
        
        print("\nAll MCP server tests completed!")
        return True
        
    except ImportError as e:
        print(f"\nCould not import MCP server (mcp package may not be installed): {e}")
        print("Install with: pip install mcp")
        return False
    except Exception as e:
        print(f"\nMCP server test failed: {e}")
        return False


if __name__ == "__main__":
    print("="*60)
    print("Tensor Network Optimizer - Test Suite")
    print("="*60)
    
    julia_ok = test_julia_service()
    
    if julia_ok:
        print("\n✓ Julia service test passed!")
        test_mcp_server()
    else:
        print("\n✗ Julia service test failed!")
        print("\nTroubleshooting:")
        print("1. Make sure Julia is installed and in your PATH")
        print("2. Run: julia --project=/path/to/OMEinsumContractionOrders -e 'using Pkg; Pkg.instantiate()'")
