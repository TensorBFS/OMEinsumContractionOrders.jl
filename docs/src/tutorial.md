# Tutorial

## Example 1: Optimize contraction order without backend specified
`OMEinsumContractionOrders` can be used as a standalone package to optimize the contraction order of an einsum notation.
The first step is to construct an [`OMEinsumContractionOrders.EinCode`](@ref) object, which is a data type to represent the [einsum notation](https://numpy.org/doc/stable/reference/generated/numpy.einsum.html).
```@repl tutorial
using OMEinsumContractionOrders, Graphs, KaHyPar
function random_regular_eincode(n, k; optimize=nothing, seed)
    g = Graphs.random_regular_graph(n, k; seed)
    ixs = [[minmax(e.src,e.dst)...] for e in Graphs.edges(g)]  # input indices
    iy = Int[]  # output indices (scalar output)
    return OMEinsumContractionOrders.EinCode(ixs, iy)
end
code = random_regular_eincode(100, 3; seed=42);
```

Here, we define an einsum notation with 3-regular graph topology. The vertices correspond to indices. On each edge, we specify a rank 2 tensor that associates with the vertices it connects. The output is a scalar.

One can use [`contraction_complexity`](@ref) function to get the time, space and rewrite cost for contracting this einsum notation. This function takes two arguments: the einsum notation and a dictionary to specify the size of the variables. Here, we use the [`uniformsize`](@ref) function to specify that all variables have the same size $2$. This function returns a dictionary that maps each variable to $2$: `Dict(i=>2 for i in uniquelabels(code))`, where [`OMEinsumContractionOrders.uniquelabels`](@ref) returns the unique labels in the einsum notation.

```@repl tutorial
size_dict = uniformsize(code, 2)
contraction_complexity(code, size_dict)
```

Since we do not specify a contraction order, the direct contraction corresponds to brute force enumeration and costs $2^{\text{number of vertices}}$ operations. No space is required to store the intermediate contraction result and the space complexity is $0$. The read-write complexity corresponds to how many element-wise read and write operations are required to perform contraction.

The order of contraction is optimized by the [`optimize_code`](@ref) function. It takes three arguments: `code`, `size_dict`, and `optimizer`. The `optimizer` argument is the optimizer to be used. The available optimizers are listed in the [optimizers](optimizers.md) page.

```@repl tutorial
score = ScoreFunction(sc_target=10, sc_weight=3.0)
optcode_tree = optimize_code(code, uniformsize(code, 2),
	TreeSA(; βs=0.1:0.1:10, ntrials=2, niters=20, score);
	slicer=TreeSASlicer(; score)
    );
contraction_complexity(optcode_tree, size_dict)
optcode_tree.slicing
```

The `optimize_code` function returns the optimized contraction order. The optimized contraction order is a [`OMEinsumContractionOrders.NestedEinsum`](@ref) object, which is a data type to represent the nested einsum notation.

## Example 2: Use it in `OMEinsum`

`OMEinsumContractionOrders` is shipped with [`OMEinsum`](https://github.com/under-Peter/OMEinsum.jl) package, which is a powerful package for tensor network contraction (or einsum contraction). You can use it to optimize the contraction order of an `OMEinsum` notation.

```@repl tutorial
using OMEinsum

code = ein"ij, jk, kl, il->"

optimized_code = optimize_code(code, uniformsize(code, 2), TreeSA())
```

## Example 3: Visualization

The visualization is provided by the `LuxorTensorPlot` extension. To use it, just load the extension by `using LuxorGraphPlot`.
```julia
pkg> add OMEinsumContractionOrders, LuxorGraphPlot

julia> using OMEinsumContractionOrders, LuxorGraphPlot
```

The extension provides the following two functions, [`viz_eins`](@ref) and [`viz_contraction`](@ref), where the former will plot the tensor network as a graph, and the latter will generate a video or gif of the contraction process.

Here is an example:
```julia-repl tutorial
julia> using OMEinsumContractionOrders, LuxorGraphPlot

julia> eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], ['a'])
ab, acd, bcef, e, df -> a

julia> viz_eins(eincode, filename = "eins.png")

julia> nested_eins = optimize_code(eincode, uniformsize(eincode, 2), GreedyMethod())
ab, ab -> a
├─ ab
└─ acf, bcf -> ab
   ├─ acd, df -> acf
   │  ├─ acd
   │  └─ df
   └─ bcef, e -> bcf
      ├─ bcef
      └─ e


julia> viz_contraction(nested_code)
[ Info: Generating frames, 7 frames in total
[ Info: Creating video at: /var/folders/3y/xl2h1bxj4ql27p01nl5hrrnc0000gn/T/jl_SiSvrH/contraction.mp4
"/var/folders/3y/xl2h1bxj4ql27p01nl5hrrnc0000gn/T/jl_SiSvrH/contraction.mp4"
```

The resulting image and video will be saved in the current working directory, and the image is shown below:
```@raw html
<div style="text-align:center">
	<img src="/assets/eins.png" alt="Image" width="40%" />
</div>
```
The large white nodes represent the tensors, and the small colored nodes represent the indices, red for closed indices and green for open indices.
