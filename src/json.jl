
function writejson(filename::AbstractString, ne::Union{NestedEinsum, SlicedEinsum})
    dict = _todict(ne)
    open(filename, "w") do f
        JSON.print(f, dict, 0)
    end
end
function _todict(ne::SlicedEinsum)
    dict = _todict(ne.eins)
    dict["slices"] = ne.slicing
    return dict
end
function _todict(ne::NestedEinsum)
    LT = labeltype(ne)
    dict = Dict{String,Any}("label-type"=>LT, "inputs"=>getixsv(ne), "output"=>getiyv(ne))
    dict["tree"] = todict(ne)
    return dict
end
function readjson(filename::AbstractString)
    dict = JSON.parsefile(filename)
    return _fromdict(dict)
end
function _fromdict(dict)
    lt = dict["label-type"]
    LT = if lt == "Char"
        Char
    elseif lt ∈ ("Int64", "Int", "Int32")
        Int
    else
        error("label type `$lt` not known.")
    end
    ne = fromdict(LT, dict["tree"])
    if haskey(dict, "slices")
        return SlicedEinsum(_convert.(LT, dict["slices"]), ne)
    else
        return ne
    end
end

function todict(ne::NestedEinsum)
    dict = Dict{String,Any}()
    if isleaf(ne)
        dict["isleaf"] = true
        dict["tensorindex"] = ne.tensorindex
        return dict
    end
    dict["args"] = collect(todict.(ne.args))
    dict["eins"] = einstodict(ne.eins)
    dict["isleaf"] = false
    return dict
end
function einstodict(eins::EinCode)
    ixs = getixsv(eins)
    iy = getiyv(eins)
    return Dict("ixs"=>ixs, "iy"=>iy)
end

function fromdict(::Type{LT}, dict::Dict) where LT
    if dict["isleaf"]
        return NestedEinsum{LT}(dict["tensorindex"])
    end
    eins = einsfromdict(LT, dict["eins"])
    return NestedEinsum(fromdict.(LT, dict["args"]), eins)
end

function einsfromdict(::Type{LT}, dict::Dict) where LT
    return EinCode([collect(LT, _convert.(LT, ix)) for ix in dict["ixs"]], collect(LT, _convert.(LT, dict["iy"])))
end

_convert(::Type{LT}, x) where LT = convert(LT, x)
_convert(::Type{Char}, x::String) where LT = (@assert length(x)==1; x[1])
