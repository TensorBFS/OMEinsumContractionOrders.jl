# Reference

## Data structures and interfaces
```@autodocs
Modules = [OMEinsumContractionOrders]
Pages = ["Core.jl"]
```

## Time and space complexity
```@autodocs
Modules = [OMEinsumContractionOrders]
Pages = ["complexity.jl"]
```

## Contraction order optimizers
```@autodocs
Modules = [OMEinsumContractionOrders]
Pages = ["interfaces.jl", "greedy.jl", "treesa.jl", "treewidth.jl", "hypernd.jl", "kahypar.jl", "sabipartite.jl"]
```

## Preprocessing code
Some optimizers (e.g. [`TreeSA`](@ref)) are too costly to run. We provide a preprocessing step to reduce the time of calling optimizers.
```@autodocs
Modules = [OMEinsumContractionOrders]
Pages = ["simplify.jl"]
```

## Dump and load contraction orders
```@autodocs
Modules = [OMEinsumContractionOrders]
Pages = ["json.jl"]
```

## Visualization
Requires `using LuxorGraphPlot` to load the extension.

```@autodocs
Modules = [OMEinsumContractionOrders]
Pages = ["visualization.jl"]
```