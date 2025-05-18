# Examples

## Setup

Navigate to the repo root and run:
```bash
$ make init  # or make update
```
It will install the dependencies for your examples environment.

## Run examples

1. [nqueens.jl](nqueens.jl) - N-Queens problem.
    ```bash
    $ julia --project=examples examples/nqueens.jl
    ```
    Source: https://github.com/nzy1997/TensorNQueens.jl/
2. [circuit.jl](circuit.jl) - Quantum circuit simulation.
    ```bash
    $ julia --project=examples examples/circuit.jl
    ```
    Source: https://github.com/CodingThrust/AMAT5315-2025Spring-Homeworks/tree/main/hw8
3. [slicing_multigpu.jl](slicing_multigpu.jl) - Multi-GPU simulation of a quantum circuit (parallel in slices).
    ```bash
    $ julia -e 'using Pkg; Pkg.add("CUDA")'
    $ julia --project=examples examples/slicing_multigpu.jl
    ```
    Here, we first install `CUDA` package in the global environment, and then run the example.