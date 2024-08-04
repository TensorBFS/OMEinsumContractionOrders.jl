function OMEinsumContractionOrders.ein2elimination(code::NestedEinsum{T}) where{T}
    elimination_order = Vector{T}()
    _ein2elimination!(code, elimination_order)
    return elimination_order
end

function OMEinsumContractionOrders.ein2elimination(code::SlicedEinsum{T, NestedEinsum{T}}) where{T}
    elimination_order = Vector{T}()
    _ein2elimination!(code.eins, elimination_order)
    # the slicing indices are eliminated at the end
    return vcat(elimination_order, code.slicing)
end

function _ein2elimination!(code::NestedEinsum{T}, elimination_order::Vector{T}) where{T}
    if code.tensorindex == -1
        for arg in code.args
            _ein2elimination!(arg, elimination_order)
        end
        iy = unique(vcat(getiyv(code.eins)...))
        for ix in unique(vcat(getixsv(code.eins)...))
            if !(ix in iy) && !(ix in elimination_order)
                push!(elimination_order, ix)
            end
        end
    end
    return elimination_order
end

function elimination_frame(gviz::GraphViz, tng::TensorNetworkGraph{TG, TL}, elimination_order::Vector{TL}, i::Int; filename = nothing) where{TG, TL}
    gviz2 = deepcopy(gviz)
    for j in 1:i
        id = _get_key(tng.indices_labels, elimination_order[j])
        gviz2.vertex_colors[id] = (0.5, 0.5, 0.5, 0.5)
    end
    return show_graph(gviz2, filename = filename)
end

function OMEinsumContractionOrders.viz_contraction(code::T, args...; kwargs...) where{T <: AbstractEinsum}
    throw(ArgumentError("Only NestedEinsum and SlicedEinsum{T, NestedEinsum{T}} have contraction order"))
end

"""
    viz_contraction(code::Union{NestedEinsum, SlicedEinsum}; locs=StressLayout(), framerate=10, filename=tempname() * ".mp4", show_progress=true)

Visualize the contraction process of a tensor network.

### Arguments
- `code`: The tensor network to visualize.

### Keyword Arguments
- `locs`: The coordinates or layout algorithm to use for positioning the nodes in the graph. Default is `StressLayout()`.
- `framerate`: The frame rate of the animation. Default is `10`.
- `filename`: The name of the output file, with `.gif` or `.mp4` extension. Default is a temporary file with `.mp4` extension.
- `show_progress`: Whether to show progress information. Default is `true`.

# Returns
- the path of the generated file.
"""
function OMEinsumContractionOrders.viz_contraction(
        code::Union{NestedEinsum, SlicedEinsum}; 
        locs=StressLayout(),
        framerate = 10,
        filename::String = tempname() * ".mp4",
        show_progress::Bool = true)

    # analyze the output format
    @assert endswith(filename, ".gif") || endswith(filename, ".mp4") "Unsupported file format: $filename, only :gif and :mp4 are supported"
    tempdirectory = mktempdir()

    # generate the frames
    elimination_order = ein2elimination(code)
    tng = TensorNetworkGraph(ein2hypergraph(code))
    gviz = GraphViz(tng, locs)

    le = length(elimination_order)
    for i in 0:le
        show_progress && @info "Frame $(i + 1) of $(le + 1)"
        fig_name = joinpath(tempdirectory, "$(lpad(i+1, 10, "0")).png")
        elimination_frame(gviz, tng, elimination_order, i; filename = fig_name)
    end

    if endswith(filename, ".gif")
        Luxor.FFMPEG.exe(`-loglevel panic -r $(framerate) -f image2 -i $(tempdirectory)/%10d.png -filter_complex "[0:v] split [a][b]; [a] palettegen=stats_mode=full:reserve_transparent=on:transparency_color=FFFFFF [p]; [b][p] paletteuse=new=1:alpha_threshold=128" -y $filename`)
    else
        Luxor.FFMPEG.ffmpeg_exe(`
            -loglevel panic
            -r $(framerate) 
            -f image2 
            -i $(tempdirectory)/%10d.png 
            -c:v libx264 
            -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"
            -r $(framerate) 
            -pix_fmt yuv420p 
            -y $filename`)
    end
    show_progress && @info "Generated output at: $filename"
    return filename
end