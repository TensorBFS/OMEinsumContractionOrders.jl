using GenericTensorNetworks, Graphs, Metis, OMEinsumContractionOrders

@info "Random 3-regular graph of size 200"
rg3 = Graphs.random_regular_graph(200, 3; seed=42)

problem = IndependentSet(rg3)
optimizer = Treewidth(; alg=ND(MF(), METISND(; ufactor=150); limit=200, level=6))

@info "Using optimizer: $(optimizer)"
net = GenericTensorNetwork(problem; optimizer)
@info "The generated tensor network is: $(net)"

@info "Random diagonal coupled graph (King's subgraph, or KSG) of size 40"
ksg = GenericTensorNetworks.random_diagonal_coupled_graph(40, 40, 0.8)

problem = IndependentSet(ksg)
net = GenericTensorNetwork(problem; optimizer)
@info "The generated tensor network is: $(net)"