module KaHyParExt

using OMEinsumContractionOrders: KaHyParBipartite, SparseMatrixCSC, group_sc, induced_subhypergraph, convert2int
import KaHyPar
import OMEinsumContractionOrders: bipartite_sc
using Suppressor: @suppress

function bipartite_sc(bipartiter::KaHyParBipartite, adj::SparseMatrixCSC, vertices, log2_sizes)
    n_v = length(vertices)
    subgraph, remaining_edges = induced_subhypergraph(adj, vertices)
    hypergraph = KaHyPar.HyperGraph(subgraph, ones(n_v), convert2int(log2_sizes[remaining_edges]))
    local parts
    min_sc = 999999
    for imbalance in bipartiter.imbalances
        parts = @suppress KaHyPar.partition(hypergraph, 2; imbalance=imbalance, configuration=:edge_cut)
        part0 = vertices[parts .== 0]
        part1 = vertices[parts .== 1]
        sc0, sc1 = group_sc(adj, part0, log2_sizes), group_sc(adj, part1, log2_sizes)
        sc = max(sc0, sc1)
        min_sc = min(sc, min_sc)
        @debug "imbalance $imbalance: sc = $sc, group = ($(length(part0)), $(length(part1)))"
        if sc <= bipartiter.sc_target
            return part0, part1
        end
    end
    error("fail to find a valid partition for `sc_target = $(bipartiter.sc_target)`, got minimum value `$min_sc` (imbalances = $(bipartiter.imbalances))")
end

@info "`OMEinsumContractionOrders` loads `KaHyParExt` extension successfully."
end