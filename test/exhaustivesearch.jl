using OMEinsumContractionOrders
using OMEinsumContractionOrders: EinCode, NestedEinsum, ContractionTree, IncidenceList,
    getixsv, getiyv, parse_eincode!, optimize_exhaustive
using OMEinsum: decorate
using Test, Random

# Independent reference: enumerate every full binary contraction tree over the
# tensor indices, score each with the trusted `flop`, and return the minimum.
# This validates that the DP picks an optimal *topology* without sharing code
# with the DP itself.
function all_binary_trees(leaves::Vector{Int})
    length(leaves) == 1 && return Any[leaves[1]]
    res = Any[]
    first = leaves[1]
    rest = leaves[2:end]
    m = length(rest)
    for mask in 1:(2^m - 1)  # nonempty subset of `rest` forms the right subtree
        right = [rest[i] for i in 1:m if (mask >> (i - 1)) & 1 == 1]
        left = [first; [rest[i] for i in 1:m if (mask >> (i - 1)) & 1 == 0]]
        for tl in all_binary_trees(left), tr in all_binary_trees(right)
            push!(res, ContractionTree(tl, tr))
        end
    end
    return res
end

function brute_force_min_flop(code::EinCode, size_dict)
    ixs, iy = getixsv(code), getiyv(code)
    n = length(ixs)
    best = Inf
    for tree in all_binary_trees(collect(1:n))
        il = IncidenceList(Dict(i => ixs[i] for i in 1:n); openedges = iy)
        ne = parse_eincode!(il, tree, collect(1:n), size_dict)[2]
        best = min(best, Float64(flop(ne, size_dict)))
    end
    return best
end

@testset "exhaustive search: optimality" begin
    optimizer = ExhaustiveSearch()
    size_dict = Dict([c => (1 << i) for (i, c) in enumerate(['a', 'b', 'c', 'd', 'e', 'f'])]...)

    # fully-contracted network (all indices appear exactly twice, empty output)
    eincode = EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], Char[])
    optcode = optimize_code(eincode, size_dict, optimizer)
    cc = contraction_complexity(optcode, size_dict)
    @test cc.tc ≈ log2(flop(optcode, size_dict))
    @test flop(optcode, size_dict) == brute_force_min_flop(eincode, size_dict)
    # never worse than greedy
    greedycode = optimize_code(eincode, size_dict, GreedyMethod())
    @test flop(optcode, size_dict) <= flop(greedycode, size_dict)

    # matrix chain with open (once-appearing) output indices
    sd2 = Dict('a' => 2, 'b' => 3, 'c' => 4, 'd' => 5, 'e' => 6)
    chain = EinCode([['a', 'b'], ['b', 'c'], ['c', 'd'], ['d', 'e']], ['a', 'e'])
    optchain = optimize_code(chain, sd2, optimizer)
    @test flop(optchain, sd2) == brute_force_min_flop(chain, sd2)
end

@testset "exhaustive search: correctness" begin
    optimizer = ExhaustiveSearch()
    size_dict = Dict([c => (1 << i) for (i, c) in enumerate(['a', 'b', 'c', 'd', 'e', 'f'])]...)

    eincode = EinCode([['a', 'b'], ['a', 'c', 'd'], ['b', 'c', 'e', 'f'], ['e'], ['d', 'f']], Char[])
    tensors = [rand([size_dict[j] for j in ixs]...) for ixs in getixsv(eincode)]
    optcode = optimize_code(eincode, size_dict, optimizer)
    @test decorate(eincode)(tensors...) ≈ decorate(optcode)(tensors...)

    # output preserved at the root
    @test sort(getiyv(optcode)) == sort(getiyv(eincode))
end

@testset "exhaustive search: open indices and disconnected" begin
    optimizer = ExhaustiveSearch()
    sd = Dict(c => 2 + i for (i, c) in enumerate(['a', 'b', 'c', 'd', 'e', 'f']))

    # two disconnected chains, combined via outer product
    eincode = EinCode([['a', 'b'], ['b', 'c'], ['d', 'e'], ['e', 'f']], ['a', 'c', 'd', 'f'])
    tensors = [rand([sd[j] for j in ixs]...) for ixs in getixsv(eincode)]
    optcode = optimize_code(eincode, sd, optimizer)
    @test sort(getiyv(optcode)) == sort(getiyv(eincode))
    @test decorate(eincode)(tensors...) ≈ decorate(optcode)(tensors...)
    @test flop(optcode, sd) == brute_force_min_flop(eincode, sd)

    # trivial small cases
    one_t = EinCode([['a', 'b']], ['a', 'b'])
    @test optimize_code(one_t, sd, optimizer) isa NestedEinsum
    two_t = EinCode([['a', 'b'], ['b', 'c']], ['a', 'c'])
    t2 = [rand(sd['a'], sd['b']), rand(sd['b'], sd['c'])]
    @test decorate(two_t)(t2...) ≈ decorate(optimize_code(two_t, sd, optimizer))(t2...)
end

@testset "exhaustive search: hyperedges and batch outputs" begin
    optimizer = ExhaustiveSearch()
    sd = Dict(c => 2 + i for (i, c) in enumerate(['a', 'b', 'c', 'd', 'e']))

    # hyperedge: index 'a' is shared by three tensors and summed out
    hyper = EinCode([['a', 'b'], ['a', 'c'], ['a', 'd']], ['b', 'c', 'd'])
    th = [rand([sd[j] for j in ixs]...) for ixs in getixsv(hyper)]
    oc = optimize_code(hyper, sd, optimizer)
    @test sort(getiyv(oc)) == sort(getiyv(hyper))
    @test decorate(hyper)(th...) ≈ decorate(oc)(th...)
    @test flop(oc, sd) == brute_force_min_flop(hyper, sd)

    # batch/diagonal output: index 'a' is shared by two tensors and kept
    batch = EinCode([['a', 'b'], ['a', 'c'], ['c', 'd']], ['a', 'b', 'd'])
    tb = [rand([sd[j] for j in ixs]...) for ixs in getixsv(batch)]
    ocb = optimize_code(batch, sd, optimizer)
    @test sort(getiyv(ocb)) == sort(getiyv(batch))
    @test decorate(batch)(tb...) ≈ decorate(ocb)(tb...)
    @test flop(ocb, sd) == brute_force_min_flop(batch, sd)
end

@testset "exhaustive search: dimension-1 indices" begin
    optimizer = ExhaustiveSearch()

    # some contracted indices have dimension 1
    sd = Dict('a' => 4, 'b' => 1, 'c' => 3, 'd' => 1, 'e' => 2)
    mixed = EinCode([['a', 'b'], ['b', 'c'], ['c', 'd'], ['d', 'e']], ['a', 'e'])
    tm = [rand([sd[j] for j in ix]...) for ix in getixsv(mixed)]
    ocm = optimize_code(mixed, sd, optimizer)
    @test decorate(mixed)(tm...) ≈ decorate(ocm)(tm...)
    @test flop(ocm, sd) == brute_force_min_flop(mixed, sd)

    # every index has dimension 1 (the costfac == 1 cost-cap edge case)
    sd1 = Dict(c => 1 for c in ['a', 'b', 'c', 'd'])
    allone = EinCode([['a', 'b'], ['b', 'c'], ['c', 'd'], ['a', 'd']], Char[])
    t1 = [rand([sd1[j] for j in ix]...) for ix in getixsv(allone)]
    oc1 = optimize_code(allone, sd1, optimizer)
    @test decorate(allone)(t1...) ≈ decorate(oc1)(t1...)
    @test flop(oc1, sd1) == brute_force_min_flop(allone, sd1)

    # dimension-1 output indices
    sd2 = Dict('a' => 1, 'b' => 2, 'c' => 2, 'd' => 1)
    out1 = EinCode([['a', 'b'], ['b', 'c'], ['c', 'd']], ['a', 'd'])
    to = [rand([sd2[j] for j in ix]...) for ix in getixsv(out1)]
    oco = optimize_code(out1, sd2, optimizer)
    @test sort(getiyv(oco)) == sort(getiyv(out1))
    @test decorate(out1)(to...) ≈ decorate(oco)(to...)
end

@testset "exhaustive search: scope errors" begin
    optimizer = ExhaustiveSearch()
    sd = Dict(c => 2 for c in ['a', 'b', 'c', 'd'])

    # partial trace within a single tensor
    selftrace = EinCode([['a', 'a', 'b'], ['b', 'c'], ['c', 'd']], ['d'])
    @test_throws ArgumentError optimize_code(selftrace, sd, optimizer)

    # dangling index: 'x' appears once but is not an output
    dangling = EinCode([['a', 'b', 'x'], ['b', 'c'], ['c', 'd']], ['a', 'd'])
    @test_throws ArgumentError optimize_code(dangling, sd, optimizer)
end
