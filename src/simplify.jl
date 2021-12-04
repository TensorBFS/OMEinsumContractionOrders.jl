using OMEinsum: StaticEinCode, DynamicEinCode
using OMEinsum.ContractionOrder: IncidenceList, ContractionTree
export merge_vectors, embed_simplifier
export merge_greedy

struct NetworkSimplifier
    operations::Vector{NestedEinsum}
end

function merge_vectors(code::EinCode)
    ET = code isa DynamicEinCode ? typeof(code) : StaticEinCode
    ixs = OMEinsum.getixs(code)
    mask = trues(length(ixs))
    ops = [NestedEinsum{ET}(i) for i=1:length(ixs)]
    for i in 1:length(ixs)
        if length(ixs[i]) == 1
            for j in 1:length(ixs)
                if i!=j && mask[j] && ixs[i][1] ∈ ixs[j]  # merge i to j
                    mask[i] = false
                    ops[j] = NestedEinsum([ops[i], ops[j]],
                        _similar(code, (ixs[i], ixs[j]), ixs[j]))
                    break
                end
            end
        end
    end
    newcode = _similar(code, ixs[mask], OMEinsum.getiy(code))
    return NetworkSimplifier(ops[mask]), newcode
end
_similar(::DynamicEinCode, ixs, iy) = DynamicEinCode(collect(collect.(ixs)), collect(iy))
_similar(::StaticEinCode, ixs, iy) = StaticEinCode{(Tuple.(ixs)...,), (iy...,)}()

function merge_greedy(code::EinCode, size_dict; threshhold=-1e-12)
    ixs, iy, L = getixsv(code), getiyv(code), OMEinsum.labeltype(code)
    ET = code isa DynamicEinCode ? typeof(code) : StaticEinCode
    log2_edge_sizes = Dict{L,Float64}()
    for (k, v) in size_dict
        log2_edge_sizes[k] = log2(v)
    end
    incidence_list = IncidenceList(Dict([i=>ixs[i] for i=1:length(ixs)]); openedges=iy)
    n = ContractionOrder.nv(incidence_list)
    if n == 0
        return nothing
    elseif n == 1
        return collect(ContractionOrder.vertices(incidence_list))[1]
    end
    tree = Dict{Int,NestedEinsum}([v=>NestedEinsum{ET}(v) for v in ContractionOrder.vertices(incidence_list)])
    cost_values = ContractionOrder.evaluate_costs(MinSpaceDiff(), incidence_list, log2_edge_sizes)
    while true
        if length(cost_values) == 0
            return _buildsimplifier(code, tree, incidence_list)
        end
        v, pair = findmin(cost_values)
        if v <= threshhold
            _, _, c = ContractionOrder.contract_pair!(incidence_list, pair..., log2_edge_sizes)
            tree[pair[1]] = NestedEinsum((tree[pair[1]], tree[pair[2]]), _similar(code, c.first, c.second))
            if ContractionOrder.nv(incidence_list) <= 1
                return _buildsimplifier(code, tree, incidence_list)
            end
            ContractionOrder.update_costs!(cost_values, pair..., MinSpaceDiff(), incidence_list, log2_edge_sizes)
        else
            return _buildsimplifier(code, tree, incidence_list)
        end
    end
end
function _buildsimplifier(code, tree, incidence_list)
    vertices = sort!(collect(keys(incidence_list.v2e)))
    ixs = [incidence_list.v2e[v] for v in vertices]
    iy = incidence_list.openedges
    NetworkSimplifier([tree[v] for v in vertices]), _similar(code, ixs, iy)
end

function apply_simplifier(s::NetworkSimplifier, xs)
    map(s.operations) do op
        return op(xs...)
    end
end
(s::NetworkSimplifier)(xs...) = apply_simplifier(s, xs)

function embed_simplifier(code::NestedEinsum, simplifier)
    if OMEinsum.isleaf(code)
        op = simplifier.operations[code.tensorindex]
        return op
    else
        return NestedEinsum(map(code.args) do arg
            embed_simplifier(arg, simplifier)
        end, code.eins)
    end
end

embed_simplifier(code::SlicedEinsum, simplifier) = SlicedEinsum(code.slicing, embed_simplifier(code.eins, simplifier))