using OMEinsum: StaticEinCode, DynamicEinCode
export merge_vectors, embed_simplifier

struct VectorRemover
    operations::Vector{NestedEinsum}
end

function merge_vectors(code::EinCode)
    ixs = OMEinsum.getixs(code)
    mask = trues(length(ixs))
    ops = [NestedEinsum((i,), _similar(code, (ix,), ix)) for (i,ix) in enumerate(ixs)]
    for i in 1:length(ixs)
        if length(ixs[i]) == 1
            for j in 1:length(ixs)
                if i!=j && mask[j] && ixs[i][1] ∈ ixs[j]  # merge i to j
                    mask[i] = false
                    ops[j] = NestedEinsum((ops[i], ops[j]),
                        _similar(code, (ixs[i], ixs[j]), ixs[j]))
                    break
                end
            end
        end
    end
    newcode = _similar(code, ixs[mask], OMEinsum.getiy(code))
    return VectorRemover(ops[mask]), newcode
end
_similar(::DynamicEinCode, ixs, iy) = DynamicEinCode(collect(ixs), iy)
_similar(::StaticEinCode, ixs, iy) = StaticEinCode{ixs, iy}()

function apply_simplifier(s::VectorRemover, xs)
    map(s.operations) do op
        return op(xs...)
    end
end
(s::VectorRemover)(xs...) = apply_simplifier(s, xs)

function embed_simplifier(code::NestedEinsum, simplifier)
    NestedEinsum(map(code.args) do arg
        embed_simplifier(arg, simplifier)
    end, code.eins)
end

function embed_simplifier(code::Integer, simplifier::VectorRemover)
    op = simplifier.operations[code]
    return unwrap_identity(op)
end

function unwrap_identity(op::NestedEinsum)
    ixs = OMEinsum.getixs(op.eins)
    if length(ixs) == 1 && ixs[1] == OMEinsum.getiy(op.eins)  # identity
        if op.args[1] isa Integer
            return op.args[1]
        else
            return unwrap_identity(op.args[1])
        end
    else
        return NestedEinsum(unwrap_identity.(op.args), op.eins)
    end
end
