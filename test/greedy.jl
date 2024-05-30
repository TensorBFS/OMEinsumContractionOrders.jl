using OMEinsumContractionOrders
using OMEinsumContractionOrders: analyze_contraction, contract_pair!, evaluate_costs, contract_tree!, log2sumexp2, parse_tree
using OMEinsumContractionOrders: IncidenceList, analyze_contraction, LegInfo, tree_greedy, parse_eincode, optimize_greedy
using TropicalNumbers

using Test, Random
@testset "analyze contraction" begin
    incidence_list = IncidenceList(Dict('A' => ['a', 'b', 'k', 'o', 'f'], 'B'=>['a', 'c', 'd', 'm', 'f'], 'C'=>['b', 'c', 'e', 'f'], 'D'=>['e'], 'E'=>['d', 'f']), openedges=['c', 'f', 'o'])
    info = analyze_contraction(incidence_list, 'A', 'B')
    @test Set(info.l1) == Set(['k'])
    @test Set(info.l2) == Set(['m'])
    @test Set(info.l12) == Set(['a'])
    @test Set(info.l01) == Set(['b','o'])
    @test Set(info.l02) == Set(['c', 'd'])
    @test Set(info.l012) == Set(['f'])
end

@testset "tree greedy" begin
    Random.seed!(2)
    incidence_list = IncidenceList(Dict('A' => ['a', 'b'], 'B'=>['a', 'c', 'd'], 'C'=>['b', 'c', 'e', 'f'], 'D'=>['e'], 'E'=>['d', 'f']))
    log2_edge_sizes = Dict([c=>i for (i,c) in enumerate(['a', 'b', 'c', 'd', 'e', 'f'])]...)
    edge_sizes = Dict([c=>(1<<i) for (i,c) in enumerate(['a', 'b', 'c', 'd', 'e', 'f'])]...)
    il = copy(incidence_list)
    contract_pair!(il, 'A', 'B', log2_edge_sizes)
    target = IncidenceList(Dict('A' => ['b', 'c', 'd'], 'C'=>['b', 'c', 'e', 'f'], 'D'=>['e'], 'E'=>['d', 'f']))
    @test il.v2e == target.v2e
    @test length(target.e2v) == length(il.e2v)
    for (k,v) in il.e2v
        @test sort(target.e2v[k]) == sort(v)
    end
    costs = evaluate_costs(0.0, incidence_list, log2_edge_sizes)
    @test costs == Dict(('A', 'B')=>(2^9), ('A', 'C')=>2^15, ('B','C')=>2^18, ('B','E')=>2^10, ('C','D')=>2^11, ('C', 'E')=>2^14)
    tree, log2_tcs, log2_scs = tree_greedy(incidence_list, log2_edge_sizes)
    tcs_, scs_ = [], []
    contract_tree!(copy(incidence_list), tree, log2_edge_sizes, tcs_, scs_)
    @test all((log2sumexp2(tcs_), maximum(scs_)) .<= (log2(exp2(10)+exp2(16)+exp2(15)+exp2(9)), 11))
    vertices = ['A', 'B', 'C', 'D', 'E']
    optcode1 = parse_eincode(incidence_list, tree, vertices=vertices)
    @test optcode1 isa OMEinsumContractionOrders.NestedEinsum
    tree2 = parse_tree(optcode1, vertices)
    @test tree2 == tree

    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], Vector{Char}())
    size_dict = Dict([c=>(1<<i) for (i,c) in enumerate(['a', 'b', 'c', 'd', 'e', 'f'])]...)
    Random.seed!(2)
    optcode2 = optimize_greedy(eincode, size_dict) 
    cc = contraction_complexity(optcode2, edge_sizes)
    # test flop
    @test cc.tc ≈ log2(flop(optcode2, edge_sizes))
    @test flop(OMEinsumContractionOrders.EinCode([['i']], Vector{Char}()), Dict('i'=>4)) == 4
    @test 16 <= cc.tc <= log2(exp2(10)+exp2(16)+exp2(15)+exp2(9))
    @test cc.sc == 11
    @test optcode1 == optcode2
end

@testset "fullerene" begin
    Random.seed!(123)
    function fullerene()
        φ = (1+√5)/2
        res = NTuple{3,Float64}[]
        for (x, y, z) in ((0.0, 1.0, 3φ), (1.0, 2 + φ, 2φ), (φ, 2.0, 2φ + 1.0))
            for (α, β, γ) in ((x,y,z), (y,z,x), (z,x,y))
                for loc in ((α,β,γ), (α,β,-γ), (α,-β,γ), (α,-β,-γ), (-α,β,γ), (-α,β,-γ), (-α,-β,γ), (-α,-β,-γ))
                    if loc ∉ res
                        push!(res, loc)
                    end
                end
            end
        end
        return res
    end

    # flatten nested einsum
    function _flatten(code::OMEinsumContractionOrders.NestedEinsum, iy=nothing)
        isleaf(code) && return [tensorindex(code)=>iy]
        sibs = siblings(code)
        ixs = []
        for i=1:length(sibs)
            append!(ixs, _flatten(sibs[i], (rootcode(code).ixs)[i]))
        end
        return ixs
    end

    flatten(code::OMEinsumContractionOrders.EinCode) = code
    function flatten(code::OMEinsumContractionOrders.NestedEinsum{LT}) where LT
        ixd = Dict(_flatten(code))
        OMEinsumContractionOrders.EinCode([ixd[i] for i=1:length(ixd)], collect((code.eins).iy))
    end

    isleaf(ne::OMEinsumContractionOrders.NestedEinsum) = ne.tensorindex != -1
    siblings(ne::OMEinsumContractionOrders.NestedEinsum) = ne.args
    tensorindex(ne::OMEinsumContractionOrders.NestedEinsum) = ne.tensorindex
    rootcode(ne::OMEinsumContractionOrders.NestedEinsum) = ne.eins

    c60_xy = fullerene()
    c60_edges = [[i,j] for (i,(i2,j2,k2)) in enumerate(c60_xy), (j,(i1,j1,k1)) in enumerate(c60_xy) if i<j && (i2-i1)^2+(j2-j1)^2+(k2-k1)^2 < 5.0]
    code = OMEinsumContractionOrders.EinCode(vcat(c60_edges, [[i] for i=1:60]), Vector{Int}())
    size_dict = Dict([i=>2 for i in 1:60])
    log2_edge_sizes = Dict([i=>1 for i in 1:60])
    edge_sizes = Dict([i=>2 for i in 1:60])
    cc = contraction_complexity(code, edge_sizes)
    @test cc.tc == 60
    @test cc.sc == 0
    optcode = optimize_greedy(code, size_dict)
    cc2 = contraction_complexity(optcode, edge_sizes)
    @test cc2.sc == 10
    @test flatten(optcode) == code
    @test flatten(code) == code

    optcode_hyper = optimize_greedy(code, size_dict, α = 0.0, temperature = 100.0, nrepeat = 20)
    cc3 = contraction_complexity(optcode_hyper, edge_sizes)
    @test cc3.sc <= 12
    @test flatten(optcode_hyper) == code
end