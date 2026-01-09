# Tensor Network Optimizer MCP Server

An MCP (Model Context Protocol) server that provides tensor network contraction order optimization services powered by [OMEinsumContractionOrders.jl](https://github.com/TensorBFS/OMEinsumContractionOrders.jl).

**Two implementations available:**
- **Julia (recommended)** - Pure Julia using ModelContextProtocol.jl, faster after startup
- **Python** - Uses subprocess to call Julia, easier Python ecosystem integration

## Prerequisites

- **Julia** (v1.9 or later) installed and available in PATH

## Quick Start (Julia Version - Recommended)

### 1. Install Dependencies

```bash
cd /path/to/OMEinsumContractionOrders/mcp_server
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### 2. Configure Cursor

Add to `~/.cursor/mcp.json`:

```json
{
  "mcpServers": {
    "tensor-network-optimizer": {
      "command": "julia",
      "args": [
        "--project=/path/to/OMEinsumContractionOrders/mcp_server",
        "/path/to/OMEinsumContractionOrders/mcp_server/server.jl"
      ]
    }
  }
}
```

### 3. Restart Cursor

The MCP server will start automatically.

### Manual Testing

```bash
cd /path/to/OMEinsumContractionOrders/mcp_server
julia --project=. server.jl
```

---

## Alternative: Python Version

If you prefer Python or need Python ecosystem integration:

### 1. Install Dependencies

```bash
cd /path/to/OMEinsumContractionOrders/mcp_server
uv sync
```

### 2. Configure Cursor

```json
{
  "mcpServers": {
    "tensor-network-optimizer": {
      "command": "uv",
      "args": [
        "--directory",
        "/path/to/OMEinsumContractionOrders/mcp_server",
        "run",
        "server.py"
      ]
    }
  }
}
```

---

## Available Tools

### 1. `optimize_contraction_order`

Optimize the contraction order of a tensor network.

**Parameters:**
- `inputs`: List of tensor index lists (e.g., `[["i","j"], ["j","k"], ["k","l"]]`)
- `size_dict`: Dictionary mapping indices to sizes (e.g., `{"i": 100, "j": 50}`)
- `output`: Output tensor indices (default: `[]` for scalar)
- `optimizer`: Algorithm - "greedy", "treesa", "exacttreewidth", "treewidth", "sabipartite"
- `optimizer_params`: Optional algorithm parameters
- `slicing`: Enable slicing for memory reduction (default: false)

**Example:**
```
optimize_contraction_order(
    inputs=[["i", "j"], ["j", "k"], ["k", "l"]],
    size_dict={"i": 100, "j": 50, "k": 200, "l": 100},
    output=["i", "l"],
    optimizer="treesa"
)
```

### 2. `analyze_einsum`

Parse and optimize an einsum expression.

**Parameters:**
- `einsum_string`: Einstein notation (e.g., `"ij,jk->ik"`)
- `sizes`: Size dictionary or uniform size integer

**Example:**
```
analyze_einsum("abc,bcd,cde->ae", 10)
```

### 3. `compare_optimizers`

Compare different optimization algorithms on the same network.

**Parameters:**
- `inputs`, `size_dict`, `output`: Same as `optimize_contraction_order`
- `optimizers`: List of optimizer names to compare

**Example:**
```
compare_optimizers(
    inputs=[[1, 2], [2, 3], [3, 4]],
    size_dict={"1": 10, "2": 20, "3": 30, "4": 40},
    output=[1, 4],
    optimizers=["greedy", "treesa", "treewidth"]
)
```

## Optimizers

| Optimizer | Description | Best For |
|-----------|-------------|----------|
| `greedy` | Fast greedy algorithm | Quick results, simple networks |
| `treesa` | Tree simulated annealing | Complex networks, better quality |
| `treewidth` | Heuristic treewidth | Good balance of speed/quality |
| `exacttreewidth` | Exact optimal | Small networks (<20 tensors) |
| `sabipartite` | SA bipartite | Hierarchical structures |

## Complexity Metrics

Results include three complexity measures (all in log₂ scale):
- **time_complexity**: Total FLOPs for the contraction
- **space_complexity**: Size of the largest intermediate tensor
- **read_write_complexity**: Total memory operations

## Troubleshooting

### Julia not found
Make sure Julia is installed and in your PATH:
```bash
which julia
julia --version
```

### Package not found
Run package instantiation:
```bash
cd /path/to/OMEinsumContractionOrders/mcp_server
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### First call is slow
The Julia version needs to compile on first use. Subsequent calls are fast.
Consider using `--startup-file=no` for faster startup:
```json
{
  "command": "julia",
  "args": ["--startup-file=no", "--project=...", "server.jl"]
}
```
