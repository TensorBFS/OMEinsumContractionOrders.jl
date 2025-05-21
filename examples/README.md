# Examples

## Setup

Navigate to the repo root and run:
```bash
$ make init  # or make update
```
It will install the dependencies for your examples environment.

## Run examples

### Contraction order optimization in practical cases

1. [nqueens.jl](nqueens.jl) - N-Queens problem counting.
    ```bash
    $ julia --project=examples examples/nqueens.jl
    ```
    Source: https://github.com/nzy1997/TensorNQueens.jl/
2. [circuit.jl](circuit.jl) - Quantum circuit simulation.
    ```bash
    $ julia --project=examples examples/circuit.jl
    ```
    Source: https://github.com/CodingThrust/AMAT5315-2025Spring-Homeworks/tree/main/hw8
3. [mis.jl](mis.jl) - Maximum Independent Set problem.
    ```bash
    $ julia --project=examples examples/mis.jl
    ```
4. [inference.jl](inference.jl) - Inference in a Bayesian network, the relational dataset that difficult to solve with TreeSA.
    ```bash
    $ julia --project=examples examples/inference.jl
    ```
    Source: https://github.com/TensorBFS/TensorInference.jl/issues/15
5. [qec.jl](qec.jl) - Quantum error correction.
    ```bash
    $ julia --project=examples examples/qec.jl
    ```
    TBW.

### Multi-GPU simulation of a quantum circuit

Additionally, we have a [slicing_multigpu.jl](slicing_multigpu.jl) example for multi-GPU simulation of a quantum circuit (parallel in slices).
```bash
$ julia -e 'using Pkg; Pkg.add("CUDA")'
$ julia --project=examples examples/slicing_multigpu.jl
```
Here, we first install `CUDA` package in the global environment, and then run the example.