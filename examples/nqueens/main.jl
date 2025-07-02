using OMEinsumContractionOrders.CliqueTrees, OMEinsumContractionOrders.Graphs, OMEinsumContractionOrders.AbstractTrees
using OMEinsumContractionOrders
using OMEinsumContractionOrders: EinCode

struct TensorNQ
    labels::Vector{Int}
end
struct TensorNQLattice
    lattice::Matrix{TensorNQ}
    pos10::Vector{Int}
    pos01::Vector{Int}
    pos11::Vector{Int}
end

function generate_TensorNQ_lattice(n::Int)
    # res = Matrix{TensorNQ}(undef,n,n)
    res = fill(TensorNQ(fill(-1,9)),n,n)
    label_count = 0

    pos_vec = [(0,-1), (1,-1), (-1,-1), (-1,0),(0,0), (1,0), (1,1),(-1,1),(0,1)]

    pos11 = Int[]
    pos01 = Int[]
    pos10 = Int[]
    for i in 1:n
        for j in 1:n
            index_vec = fill(-1,9)
            for (ind,pos) in enumerate(pos_vec)
                if checkbounds(Bool,res,i+pos[1],j+pos[2])
                    if res[i+pos[1],j+pos[2]].labels[10-ind] != -1
                        index_vec[ind] = res[i+pos[1],j+pos[2]].labels[10-ind]
                    else
                        label_count += 1
                        index_vec[ind] = label_count
                    end
                else
                    label_count += 1
                    index_vec[ind] = label_count
                    if ind < 5
                        push!(pos10,label_count)
                    elseif ind == 9 || ind == 6
                        push!(pos01,label_count)
                    else
                        push!(pos11,label_count)
                    end
                end
            end
            res[i,j] = TensorNQ(index_vec)
        end
    end
    return TensorNQLattice(res, pos10,pos01,pos11 ∪ [res[i,j].labels[5] for i in 1:n, j in 1:n])
end

function generate_3_tensor_network(n::Int)
    t9_lattice = generate_TensorNQ_lattice(n)
    return generate_3_tensor_network(t9_lattice)
end
function generate_3_tensor_network(t9_lattice::TensorNQLattice)
    lattice,pos10,pos01,pos11 = t9_lattice.lattice,t9_lattice.pos10,t9_lattice.pos01,t9_lattice.pos11
    t3_ixs = Vector{Vector{Int}}()
    t9_ixs = getfield.(vec(lattice),:labels)

    for vec in t9_ixs
        for i in 1:4
            push!(t3_ixs,[vec[i],vec[10-i],vec[5]])
        end
    end
    return EinCode(t3_ixs ∪ [[p] for p in pos10] ∪ [[p] for p in pos01] ∪  [[p] for p in pos11],Int[])
end

function main(optimizer)
    @info "Running N-Queens with optimizer: $(optimizer)"
    n = 28
    code = generate_3_tensor_network(n)
    time_elapsed = @elapsed optcode = optimize_code(code, uniformsize(code, 2), optimizer)
    @info "Contraction complexity: $(contraction_complexity(optcode, uniformsize(optcode, 2))), time cost: $(time_elapsed)s"
    return contraction_complexity(optcode, uniformsize(optcode, 2))
end