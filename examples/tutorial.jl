### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ c4f8bf2b-d07b-478d-bcd9-ed61ccd5abe3
using Pkg; Pkg.activate()

# ╔═╡ ee4e8d2a-04db-425a-af21-f5bc5ce6e26b
using OMEinsum

# ╔═╡ 8d093496-2fe3-419b-a88b-f4a7b7dadad8
using JunctionTrees

# ╔═╡ 744be89f-da21-455e-81b2-77278bd1ddaa
using Zygote, LinearAlgebra

# ╔═╡ 0d247952-0e51-11ed-3550-279e22cfc49b
md"# OMEinsum deep dive"

# ╔═╡ 9927eea0-64be-45cc-83e9-0e25eea1499d
md"## Specify the einsum contraction"

# ╔═╡ 553aa4bb-966c-4eba-94f5-d191fc99befe
tensor_network = ein"ab,bc,cd,def,f->"

# ╔═╡ 3fc61408-18a4-4679-8880-7fc7fff3acb6
md"One can extract inputs and output labels using `getixsv` and `getiy`"

# ╔═╡ 91660727-6a55-4e38-939d-3d27ec10cc3e
getixsv(tensor_network)

# ╔═╡ b93688af-3bb0-4a9b-81f9-125897310c49
getiyv(tensor_network)

# ╔═╡ 6197f512-0ab6-4975-b3d0-83eb72490d85
md"One can also construct `EinCode` programatically as"

# ╔═╡ e71fd0a9-e4fe-48a8-b685-c0971c7f8918
EinCode(getixsv(tensor_network), getiyv(tensor_network))

# ╔═╡ d255114a-fec8-406c-b710-20f53c682a3f
md"The label type does not has to be `Char`, it can be any type."

# ╔═╡ 64b77f52-d476-4f57-b9ad-493a031bf296
md"**Example: loading factor graph from an `uai` file**"

# ╔═╡ 55a597c3-1a74-45a0-9237-178990a7eb9d
md"In the following, we use tensor network for probabilistic modeling. The first step is loading a factor graph form an file of [UAI format](https://mroavi.github.io/JunctionTrees.jl/stable/file_formats/uai/#UAI-model-file-format). The loaded model is decribed in detail in the [documentation of JunctionTrees.jl](https://mroavi.github.io/JunctionTrees.jl/stable/usage/)."

# ╔═╡ 42ded005-d55a-423d-b642-2f142c579036
uai_folder = joinpath(pkgdir(JunctionTrees), "docs", "src", "problems", "asia")

# ╔═╡ 917d33c0-c3c2-4865-8842-60c7248f0a67
# outputs are (number of variables,
#				cardinalities of each variable,
# 				number of cliques,
# 				the factors (labelled tensors) for the cliques)
nvars, cards, nclique, factors = JunctionTrees.read_uai_file(joinpath(uai_folder, "asia.uai"))

# ╔═╡ 3a5b2c4e-aa6b-4ac8-8f43-eb7e217a85cb
# outputs are observables and its values
obsvars, obsvals = [], []#JunctionTrees.read_uai_evid_file(joinpath(uai_folder, "asia.uai.evid"))

# ╔═╡ 5af37871-593c-43f2-acd1-45eb682ca0ac
md"The first 8 inputs labels are for unity tensors."

# ╔═╡ 896d755e-6a99-4f8d-a358-aa7714a036d5
factor_graph = EinCode([[[i] for i in 1:nvars]..., [[factor.vars...] for factor in factors]...], Int[])  # labels for edge tensors

# ╔═╡ da889e9e-4fb7-473b-908c-177e8e11c83d
fixedvertices=Dict(zip(obsvars, obsvals .- 1))

# ╔═╡ 07cb8edf-f009-428e-a0e3-4e9b30974ce4
size_dict = Dict([v=>(haskey(fixedvertices, v) ? 1 : 2) for v in 1:nvars])

# ╔═╡ 3505a5b2-6c20-4643-9937-4c1af27b7398
md"## Einsum with a specific contraction order"

# ╔═╡ 441fcd01-9b84-4389-adf7-d5eef09c0ebd
md"**Example: optimize the contraction order for inference**"

# ╔═╡ 58a3f8ab-0ecf-4024-9b35-d76ebaba2b32
optimized_factor_graph = optimize_code(factor_graph, size_dict, TreeSA(ntrials=1))

# ╔═╡ 7f728828-2710-46ca-ae54-a8dc512ecce7
md"prepare input tensors"

# ╔═╡ 7a105d4c-fd9d-4756-b6e2-93362b326219
clip_size(labels, fixedvertices) = [haskey(fixedvertices, v) ? (fixedvertices[v]+1:fixedvertices[v]+1) : Colon() for v in labels]

# ╔═╡ 65be0aaa-73bf-4120-a1d2-284c2cc7ad1d
tensors = [
	[ones(haskey(fixedvertices, i) ? 1 : 2) for i=1:length(cards)]...,  # unity tensors
	[factor.vals[clip_size(factor.vars, fixedvertices)...]
	for factor in factors]...
]

# ╔═╡ 6fcc8ec7-da6d-4198-a30f-6d792f2682a7
size.(tensors)

# ╔═╡ b1083a2c-3283-45f6-aaf7-36651643d8f8
getixsv(optimized_factor_graph)

# ╔═╡ f6414390-820f-435f-9d18-2bf2b499f268
optimized_factor_graph(tensors...)

# ╔═╡ 90f8f78a-f323-466c-b501-d95db4db76c7
marginals_raw = Zygote.gradient((args...)->optimized_factor_graph(args...)[], tensors...)

# ╔═╡ 0e89da78-08c0-49fb-8ec3-559345b32b0c
marginals_raw ./ sum.(marginals_raw)

# ╔═╡ Cell order:
# ╟─0d247952-0e51-11ed-3550-279e22cfc49b
# ╠═c4f8bf2b-d07b-478d-bcd9-ed61ccd5abe3
# ╠═ee4e8d2a-04db-425a-af21-f5bc5ce6e26b
# ╠═8d093496-2fe3-419b-a88b-f4a7b7dadad8
# ╟─9927eea0-64be-45cc-83e9-0e25eea1499d
# ╠═553aa4bb-966c-4eba-94f5-d191fc99befe
# ╟─3fc61408-18a4-4679-8880-7fc7fff3acb6
# ╠═91660727-6a55-4e38-939d-3d27ec10cc3e
# ╠═b93688af-3bb0-4a9b-81f9-125897310c49
# ╟─6197f512-0ab6-4975-b3d0-83eb72490d85
# ╠═e71fd0a9-e4fe-48a8-b685-c0971c7f8918
# ╟─d255114a-fec8-406c-b710-20f53c682a3f
# ╟─64b77f52-d476-4f57-b9ad-493a031bf296
# ╟─55a597c3-1a74-45a0-9237-178990a7eb9d
# ╠═42ded005-d55a-423d-b642-2f142c579036
# ╠═917d33c0-c3c2-4865-8842-60c7248f0a67
# ╠═3a5b2c4e-aa6b-4ac8-8f43-eb7e217a85cb
# ╟─5af37871-593c-43f2-acd1-45eb682ca0ac
# ╠═896d755e-6a99-4f8d-a358-aa7714a036d5
# ╠═da889e9e-4fb7-473b-908c-177e8e11c83d
# ╠═07cb8edf-f009-428e-a0e3-4e9b30974ce4
# ╟─3505a5b2-6c20-4643-9937-4c1af27b7398
# ╟─441fcd01-9b84-4389-adf7-d5eef09c0ebd
# ╠═58a3f8ab-0ecf-4024-9b35-d76ebaba2b32
# ╟─7f728828-2710-46ca-ae54-a8dc512ecce7
# ╠═7a105d4c-fd9d-4756-b6e2-93362b326219
# ╠═65be0aaa-73bf-4120-a1d2-284c2cc7ad1d
# ╠═6fcc8ec7-da6d-4198-a30f-6d792f2682a7
# ╠═b1083a2c-3283-45f6-aaf7-36651643d8f8
# ╠═f6414390-820f-435f-9d18-2bf2b499f268
# ╠═744be89f-da21-455e-81b2-77278bd1ddaa
# ╠═90f8f78a-f323-466c-b501-d95db4db76c7
# ╠═0e89da78-08c0-49fb-8ec3-559345b32b0c
