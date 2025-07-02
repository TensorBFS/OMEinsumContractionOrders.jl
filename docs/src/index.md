```@meta
CurrentModule = OMEinsumContractionOrders
```

# OMEinsumContractionOrders

This is the documentation for [OMEinsumContractionOrders](https://github.com/TensorBFS/OMEinsumContractionOrders.jl),
a Julia package for the optimization of the contraction order of tensor networks.

**Installation** guide is available in [README.md](https://github.com/TensorBFS/OMEinsumContractionOrders.jl). You can also access its features in [OMEinsum](https://github.com/under-Peter/OMEinsum.jl), which uses it as the default contraction order optimizer.

## Related Packages

### Derived packages
- [OMEinsum.jl](https://github.com/under-Peter/OMEinsum.jl): a Julia package for tensor network contraction.
- [TensorInference.jl](https://github.com/TensorBFS/TensorInference.jl): a Julia package for exact probabilistic inference based on the tensor network representation.
- [Yao.jl](https://github.com/QuantumBFS/Yao.jl): a Julia package for quantum computing. Its tensor network based simulation backend use OMEinsumContractionOrders.jl as the default contraction order optimizer.
- [GenericTensorNetworks.jl](https://github.com/TensorBFS/GenericTensorNetworks.jl): a Julia package for computing solution space properties of computational hard problems. It is based on the generic tensor network method and its contraction order optimization backend is OMEinsumContractionOrders.jl.
- [ITensorNetworks.jl](https://github.com/mtfishman/ITensorNetworks.jl): a Julia package for physics simulation, which uses OMEinsumContractionOrders.jl as alternative contraction order optimization backend.

### Similar packages
- [Cotengra](https://cotengra.readthedocs.io/en/latest/) : a python library for contracting tensor networks or einsum expressions involving large numbers of tensors.

