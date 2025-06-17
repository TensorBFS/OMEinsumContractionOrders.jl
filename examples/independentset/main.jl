using GenericTensorNetworks, GenericTensorNetworks.Graphs, OMEinsumContractionOrders

function main(optimizer)
    @info "Running independent set with optimizer: $(optimizer)"
    return [single_run(optimizer, case) for case in [:rg3, :ksg]]
end

function single_run(optimizer, case)
    if case == :rg3
        @info "Random 3-regular graph of size 200"
        graph = Graphs.random_regular_graph(200, 3; seed=42)
    elseif case == :ksg
        @info "Random diagonal coupled graph (King's subgraph, or KSG) of size 40"
        graph = GenericTensorNetworks.random_diagonal_coupled_graph(40, 40, 0.8)
    else
        error("Invalid case: $case")
    end
    time_elapsed = @elapsed net = GenericTensorNetwork(IndependentSet(graph); optimizer)
    @info "Contraction complexity: $(contraction_complexity(net)), time cost: $(time_elapsed)s"
    return contraction_complexity(net)
end