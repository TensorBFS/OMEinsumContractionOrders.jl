using OMEinsumContractionOrders
using OMEinsumContractionOrders: analyze_contraction, contract_pair!, evaluate_costs, contract_tree!, log2sumexp2, parse_tree
using OMEinsumContractionOrders: IncidenceList, analyze_contraction, LegInfo, tree_greedy, parse_eincode, optimize_greedy
using Graphs

using Test, Random

@testset "analyze contraction" begin
    incidence_list = IncidenceList(Dict(1 => [1, 2, 11, 15, 6], 2=>[1, 3, 4, 13, 6], 3=>[2, 3, 5, 6], 4=>[5], 5=>[4, 6]), openedges=[3, 6, 15])
    info = analyze_contraction(incidence_list, 1, 2)
    @test Set(info.l1) == Set([11])
    @test Set(info.l2) == Set([13])
    @test Set(info.l12) == Set([1])
    @test Set(info.l01) == Set([2,15])
    @test Set(info.l02) == Set([3, 4])
    @test Set(info.l012) == Set([6])
end

@testset "parse eincode" begin
    incidence_list = IncidenceList(Dict(1 => [1, 2], 2=>[1, 3, 4], 3=>[2, 3, 5, 6], 4=>[5], 5=>[4, 6]))
    tree = OMEinsumContractionOrders.ContractionTree(OMEinsumContractionOrders.ContractionTree(1, 2), OMEinsumContractionOrders.ContractionTree(3, 4))
    vertices = [1, 2, 3, 4, 5, 6]
    ti, optcode = OMEinsumContractionOrders.parse_eincode!(incidence_list, tree, vertices, Dict([v=>0 for v in vertices]))
    @test ti == 1
    @test optcode isa OMEinsumContractionOrders.NestedEinsum
end

@testset "tree greedy" begin
    Random.seed!(2)
    incidence_list = IncidenceList(Dict(1 => [1, 2], 2=>[1, 3, 4], 3=>[2, 3, 5, 6], 4=>[5], 5=>[4, 6]))
    log2_edge_sizes = Dict([i=>i for i in 1:6])
    edge_sizes = Dict([i=>(1<<i) for i in 1:6])
    il = copy(incidence_list)
    contract_pair!(il, 1, 2, log2_edge_sizes)
    target = IncidenceList(Dict(1 => [2, 3, 4], 3=>[2, 3, 5, 6], 4=>[5], 5=>[4, 6]))
    @test il.v2e == target.v2e
    @test length(target.e2v) == length(il.e2v)
    for (k,v) in il.e2v
        @test sort(target.e2v[k]) == sort(v)
    end
    costs, cost_graph = evaluate_costs(0.0, incidence_list, log2_edge_sizes)
    @test costs == Dict((1, 2)=>(2^9), (1, 3)=>2^15, (2,3)=>2^18, (2,5)=>2^10, (3,4)=>2^11, (3, 5)=>2^14)
    tree, log2_tcs, log2_scs = tree_greedy(incidence_list, log2_edge_sizes)
    tcs_, scs_ = [], []
    contract_tree!(copy(incidence_list), tree, log2_edge_sizes, tcs_, scs_)
    @test all((log2sumexp2(tcs_), maximum(scs_)) .<= (log2(exp2(10)+exp2(16)+exp2(15)+exp2(9)), 11))
    vertices = [1, 2, 3, 4, 5, 6]
    optcode1 = parse_eincode(incidence_list, tree, vertices=vertices)
    @test optcode1 isa OMEinsumContractionOrders.NestedEinsum
    tree2 = parse_tree(optcode1, vertices)
    @test tree2 == tree

    eincode = OMEinsumContractionOrders.EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], Vector{Char}())
    size_dict = Dict([c=>(1<<i) for (i,c) in enumerate(['a', 'b', 'c', 'd', 'e', 'f'])]...)
    Random.seed!(2)
    optcode2 = optimize_greedy(eincode, size_dict; α=0.0, temperature=0.0)
    cc = contraction_complexity(optcode2, size_dict)
    # test flop
    @test cc.tc ≈ log2(flop(optcode2, size_dict))
    @test flop(OMEinsumContractionOrders.EinCode([['i']], Vector{Char}()), Dict('i'=>4)) == 4
    @test 16 <= cc.tc <= 17
    @test cc.sc == 11
    @test optcode1 == OMEinsumContractionOrders.convert_label(optcode2, Dict([c=>i for (i,c) in enumerate(['a', 'b', 'c', 'd', 'e', 'f'])]))
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
    optcode = optimize_greedy(code, size_dict; α=0.0, temperature=0.0)
    cc2 = contraction_complexity(optcode, edge_sizes)
    @test cc2.sc == 10
    @test flatten(optcode) == code
    @test flatten(code) == code

    cc3 = Inf
    local optcode_hyper
    for irepeat = 1:20
        optcode_hyper = optimize_greedy(code, size_dict; α=0.0, temperature=100.0)
        cc3 = contraction_complexity(optcode_hyper, edge_sizes)
        if cc3.sc < cc3.sc
            optcode_hyper = optcode_hyper
            cc3 = cc3
        end
    end
    @test cc3.sc <= 12
    @test flatten(optcode_hyper) == code
end

@testset "chain and ring" begin
    # chain
    code = OMEinsumContractionOrders.EinCode([[1,2], [2,3], [3,4], [4,5], [5,6], [6,7], [7,8], [8,9], [9,10]], Int[1, 10])
    size_dict = Dict([i=>2 for i in 1:10])
    optcode = optimize_greedy(code, size_dict; α=0.0, temperature=0.0)
    cc = contraction_complexity(optcode, size_dict)
    @test cc.sc == 2

    # ring
    ring = OMEinsumContractionOrders.EinCode([[1,2], [2,3], [3,4], [4,5], [5,6], [6,7], [7,8], [8,9], [9,10], [10,1]], Int[])
    optcode = optimize_greedy(ring, size_dict; α=0.0, temperature=0.0)
    cc = contraction_complexity(optcode, size_dict)
    @test cc.sc == 2

    # petersen
    sc_list = Float64[]
    for i=1:20
        graph = smallgraph(:tutte)
        code = OMEinsumContractionOrders.EinCode([[e.src, e.dst] for e in edges(graph)], Int[])
        size_dict = Dict([i=>2 for i in 1:nv(graph)])
        optcode = optimize_greedy(code, size_dict; α=0.0, temperature=100.0)
        cc = contraction_complexity(optcode, size_dict)
        push!(sc_list, cc.sc)
    end
    @test minimum(sc_list) == 5
end

@testset "compute_contraction_dims" begin
    # Test that compute_contraction_dims matches analyze_contraction
    using OMEinsumContractionOrders: compute_contraction_dims
    
    # Create test incidence list
    incidence_list = IncidenceList(
        Dict(1 => [1, 2, 11, 15, 6], 2 => [1, 3, 4, 13, 6], 3 => [2, 3, 5, 6], 4 => [5], 5 => [4, 6]), 
        openedges=[3, 6, 15]
    )
    log2_edge_sizes = Dict(i => Float64(i % 3 + 1) for i in [1,2,3,4,5,6,11,13,15])
    
    # Test various vertex pairs
    for (vi, vj) in [(1, 2), (2, 3), (1, 3)]
        if vi in keys(incidence_list.v2e) && vj in keys(incidence_list.v2e)
            # Get dimensions using the new function
            D1, D2, D12, D01, D02, D012 = compute_contraction_dims(incidence_list, log2_edge_sizes, vi, vj)
            
            # Compare with analyze_contraction
            legs = analyze_contraction(incidence_list, vi, vj)
            log2dim(legs_list) = isempty(legs_list) ? 0.0 : sum(l->log2_edge_sizes[l], legs_list)
            D1_ref = log2dim(legs.l1)
            D2_ref = log2dim(legs.l2)
            D12_ref = log2dim(legs.l12)
            D01_ref = log2dim(legs.l01)
            D02_ref = log2dim(legs.l02)
            D012_ref = log2dim(legs.l012)
            
            @test D1 ≈ D1_ref
            @test D2 ≈ D2_ref
            @test D12 ≈ D12_ref
            @test D01 ≈ D01_ref
            @test D02 ≈ D02_ref
            @test D012 ≈ D012_ref
        end
    end
end

@testset "greedy_loss optimization" begin
    # Original implementation using analyze_contraction (for comparison)
    function greedy_loss_with_vectors(α, incidence_list, log2_edge_sizes, vi, vj)
        log2dim(legs) = isempty(legs) ? 0 : sum(l->log2_edge_sizes[l], legs)
        legs = analyze_contraction(incidence_list, vi, vj)
        D1, D2, D12, D01, D02, D012 = log2dim.(getfield.(Ref(legs), 1:6))
        return exp2(D01+D02+D012) - α * (exp2(D01+D12+D012) + exp2(D02+D12+D012))
    end
    
    # Create a simple tensor network
    code = OMEinsumContractionOrders.EinCode([[1,2], [2,3], [3,4], [4,5], [5,1]], Int[])
    ixs = OMEinsumContractionOrders.getixsv(code)
    iy = OMEinsumContractionOrders.getiyv(code)
    size_dict = Dict([i=>2 for i in 1:5])
    log2_edge_sizes = Dict([i=>log2(size_dict[i]) for i in keys(size_dict)])
    incidence_list = IncidenceList(Dict([i=>ixs[i] for i=1:length(ixs)]); openedges=iy)
    
    # Compare optimized vs original implementation
    for α in [0.0, 0.5, 1.0]
        for vi in keys(incidence_list.v2e), vj in keys(incidence_list.v2e)
            if vi < vj
                loss_new = OMEinsumContractionOrders.greedy_loss(α, incidence_list, log2_edge_sizes, vi, vj)
                loss_old = greedy_loss_with_vectors(α, incidence_list, log2_edge_sizes, vi, vj)
                @test loss_new ≈ loss_old
            end
        end
    end
end
