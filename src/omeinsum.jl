export peak_memory

function decorate(code::ContractionOrderAlgorithms.EinCode)
    DynamicEinCode(code.ixs, code.iy)
end
function decorate(code::ContractionOrderAlgorithms.NestedEinsum{LT}) where LT
    if ContractionOrderAlgorithms.isleaf(code)
        NestedEinsum{DynamicEinCode{LT}}(code.tensorindex)
    else
        NestedEinsum(decorate.(code.args), decorate(code.eins))
    end
end
function decorate(code::ContractionOrderAlgorithms.SlicedEinsum)
    SlicedEinsum(code.slicing, decorate(code.eins))
end

function rawcode(code::EinCode)
    ContractionOrderAlgorithms.DynamicEinCode(getixsv(code), getiyv(code))
end
function rawcode(code::NestedEinsum)
    if isleaf(code)
        ContractionOrderAlgorithms.NestedEinsum{labeltype(code)}(code.tensorindex)
    else
        ContractionOrderAlgorithms.NestedEinsum(decorate.(code.args), decorate(code.eins))
    end
end
function rawcode(code::SlicedEinsum)
    ContractionOrderAlgorithms.SlicedEinsum(code.slicing, rawcode(code.eins))
end