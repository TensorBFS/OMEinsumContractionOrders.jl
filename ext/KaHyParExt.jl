module KaHyParExt

using OMEinsumContractionOrders: KaHyParBipartite, SparseMatrixCSC, group_sc, induced_subhypergraph, convert2int
import KaHyPar
import OMEinsumContractionOrders: bipartite_sc, BipartiteResult
using Suppressor: @suppress

function bipartite_sc(bipartiter::KaHyParBipartite, adj::SparseMatrixCSC, vertices, log2_sizes)
    n_v = length(vertices)
    subgraph, remaining_edges = induced_subhypergraph(adj, vertices)
    hypergraph = KaHyPar.HyperGraph(subgraph, ones(n_v), convert2int(log2_sizes[remaining_edges]))
    local best_part0, best_part1
    best_sc = typemax(Int)
    for imbalance in bipartiter.imbalances
        parts = @suppress KaHyPar.partition(hypergraph, 2; imbalance=imbalance, configuration=:edge_cut)
        part0 = vertices[parts .== 0]
        part1 = vertices[parts .== 1]
        sc0, sc1 = group_sc(adj, part0, log2_sizes), group_sc(adj, part1, log2_sizes)
        sc = max(sc0, sc1)
        if sc <= best_sc
            best_sc = sc
            best_part0, best_part1 = part0, part1
        end
        @debug "imbalance $imbalance: sc = $sc, group = ($(length(part0)), $(length(part1)))"
        if best_sc <= bipartiter.sc_target  # first encountered valid partition
            return BipartiteResult(best_part0, best_part1, best_sc, true)
        end
    end
    return BipartiteResult(best_part0, best_part1, best_sc, false)
end

@debug "`OMEinsumContractionOrders` loads `KaHyParExt` extension successfully."
end