function OMEinsumContractionOrders.ein2elimination(ein::NestedEinsum{T}) where{T}
    elimination_order = Vector{T}()
    _ein2elimination!(ein, elimination_order)
    return elimination_order
end

function OMEinsumContractionOrders.ein2elimination(ein::SlicedEinsum{T, NestedEinsum{T}}) where{T}
    elimination_order = Vector{T}()
    _ein2elimination!(ein.eins, elimination_order)
    # the slicing indices are eliminated at the end
    return vcat(elimination_order, ein.slicing)
end

function _ein2elimination!(ein::NestedEinsum{T}, elimination_order::Vector{T}) where{T}
    if ein.tensorindex == -1
        for arg in ein.args
            _ein2elimination!(arg, elimination_order)
        end
        iy = unique(vcat(getiyv(ein.eins)...))
        for ix in unique(vcat(getixsv(ein.eins)...))
            if !(ix in iy) && !(ix in elimination_order)
                push!(elimination_order, ix)
            end
        end
    end
    return elimination_order
end

function elimination_frame(GViz, tng::TensorNetworkGraph{TG, TL}, elimination_order::Vector{TL}, i::Int; filename = nothing, color = (0.5, 0.5, 0.5, 0.5)) where{TG, TL}
    GViz2 = deepcopy(GViz)
    for j in 1:i
        id = _get_key(tng.indices_labels, elimination_order[j])
        GViz2.vertex_colors[id] = color
    end
    return show_graph(GViz2, filename = filename)
end

function OMEinsumContractionOrders.viz_contraction(ein::T, args...; kwargs...) where{T <: AbstractEinsum}
    throw(ArgumentError("Only NestedEinsum and SlicedEinsum{T, NestedEinsum{T}} have contraction order"))
end

"""
    viz_contraction(ein::ET; locs=StressLayout(), framerate=30, filename="contraction", pathname=".", create_gif=false, create_video=true, color=(0.5, 0.5, 0.5, 0.5), show_progress=false) where {ET <: Union{NestedEinsum, SlicedEinsum}}

Visualize the contraction process of a tensor network.

# Arguments
- `ein::ET`: The tensor network to visualize.
- `locs`: The layout algorithm to use for positioning the nodes in the graph. Default is `StressLayout()`.
- `framerate`: The frame rate of the animation. Default is 30.
- `filename`: The name of the output files. If the filename ends with ".gif", a GIF file will be generated. If the filename ends with ".mp4", a video file will be generated. If the filename contains path, the file will be saved to the path, or the file will be saved to the tempdirectory. Default is "contraction.mp4".
- `color`: The color of the contraction lines. Default is `(0.5, 0.5, 0.5, 0.5)`.
- `show_progress`: Whether to show progress information. Default is `false`.

# Returns
- the path of the GIF or video file.
"""
function OMEinsumContractionOrders.viz_contraction(
    ein::ET; 
    locs=StressLayout(),
    framerate = 10,
    filename = "contraction.mp4",
    color = (0.5, 0.5, 0.5, 0.5),
    show_progress::Bool = false
    ) where{ET <: Union{NestedEinsum, SlicedEinsum}}

    # analyze the output format
    paths = splitpath(filename)
    file_name = paths[end]
    path_name = length(paths) > 1 ? joinpath(paths[1:end-1]...) : ""
    format = splitext(file_name)[end]
    tempdirectory = mktempdir()

    if format != ".gif" && format != ".mp4"
        throw(ArgumentError("Unsupported format $format, only .gif and .mp4 are supported"))
    end

    # generate the frames
    elimination_order = ein2elimination(ein)
    tng = TensorNetworkGraph(ein2hypergraph(ein))
    GViz = GraphViz(tng, locs)

    filecounter = 1
    le = length(elimination_order)
    @info "Generating frames, $(le + 1) frames in total"
    for i in 0:le
        if show_progress
            @info "Frame $(i + 1) of $(le + 1)"
        end
        fig_name = "$(tempdirectory)/$(lpad(filecounter, 10, "0")).png"
        elimination_frame(GViz, tng, elimination_order, i; filename = fig_name, color = color)
        filecounter += 1
    end

    if format == ".gif"
        Luxor.FFMPEG.exe(`-loglevel panic -r $(framerate) -f image2 -i $(tempdirectory)/%10d.png -filter_complex "[0:v] split [a][b]; [a] palettegen=stats_mode=full:reserve_transparent=on:transparency_color=FFFFFF [p]; [b][p] paletteuse=new=1:alpha_threshold=128" -y $(tempdirectory)/$(file_name)`)

        if !isempty(path_name)
            if !isdir(path_name)
                @error "$path_name is not a directory."
            end
            mv("$(tempdirectory)/$(file_name)", filename, force = true)
            @info("GIF is: $filename")
            giffn = filename
        else
            @info("GIF is: $(tempdirectory)/$(file_name)")
            giffn = tempdirectory * "/" * file_name
        end

        return giffn
    elseif format == ".mp4"
        if !isempty(path_name)
            if !isdir(path_name)
                @error "$path_name is not a directory."
            end
            video_path = filename
        else
            video_path = joinpath("$(tempdirectory)", "$(file_name)")
        end

        @info "Creating video at: $video_path"
        FFMPEG.ffmpeg_exe(`
            -loglevel panic
            -r $(framerate) 
            -f image2 
            -i $(tempdirectory)/%10d.png 
            -c:v libx264 
            -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"
            -r $(framerate) 
            -pix_fmt yuv420p 
            -y $(video_path)`)

        return video_path
    end
end