# Tensor Network Optimizer MCP Server

An MCP (Model Context Protocol) server that provides tensor network contraction order optimization services powered by [OMEinsumContractionOrders.jl](https://github.com/TensorBFS/OMEinsumContractionOrders.jl).

## Prerequisites

1. **Julia** (v1.9 or later) installed and available in PATH
2. **uv** - Python package manager ([install uv](https://docs.astral.sh/uv/getting-started/installation/))
3. **OMEinsumContractionOrders.jl** package installed in Julia

### Installing Julia Dependencies

```bash
cd /path/to/OMEinsumContractionOrders
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Installing Python Dependencies

```bash
cd mcp_server
uv sync
```

## Usage in Cursor

Add the following to your Cursor MCP settings (`~/.cursor/mcp.json`):

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

**Note:** Replace `/path/to/OMEinsumContractionOrders` with the actual path to your installation.

### Running Manually for Testing

```bash
cd /path/to/OMEinsumContractionOrders/mcp_server
uv run server.py
```

Or run the test script:

```bash
uv run test_service.py
```

## Available Tools

### 1. `optimize_contraction_order`

Optimize the contraction order of a tensor network.

**Parameters:**
- `inputs`: List of tensor index lists (e.g., `[[1, 2], [2, 3], [3, 4]]`)
- `size_dict`: Dictionary mapping indices to sizes (e.g., `{"1": 10, "2": 20}`)
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

## Resources

### `optimizer://docs`

Access documentation about available optimizers and their parameters.

## Complexity Metrics

Results include three complexity measures (all in log2 scale):
- **time_complexity**: Total FLOPs for the contraction
- **space_complexity**: Size of the largest intermediate tensor
- **read_write_complexity**: Total memory operations

## Troubleshooting

### Julia not found
Make sure Julia is installed and in your PATH. You can specify the full path to Julia by modifying the `run_julia_optimization` function in `server.py`.

### Package not found
Run the Julia package instantiation:
```bash
julia --project=/path/to/OMEinsumContractionOrders -e 'using Pkg; Pkg.instantiate()'
```

### Timeout errors
For large networks, increase the timeout in `server.py` (default: 300 seconds).
