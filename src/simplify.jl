using OMEinsum: StaticEinCode, DynamicEinCode
export merge_vectors

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
    newcode = _similar(code, ixs[mask], getiy(code))
    return VectorRemover(ops[mask]), newcode
end
_similar(::DynamicEinCode, ixs, iy) = DynamicEinCode(ixs, iy)
_similar(::StaticEinCode, ixs, iy) = StaticEinCode{ixs, iy}()

function apply_simplifier(s::VectorRemover, xs)
    map(s.operations) do op
        return op(xs...)
    end
end
(s::VectorRemover)(xs...) = apply_simplifier(s, xs)
