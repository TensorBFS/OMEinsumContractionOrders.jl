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
- `filename`: The base name of the output files. Default is "contraction".
- `pathname`: The directory path to save the output files. Default is the current directory.
- `create_gif`: Whether to create a GIF animation. Default is `false`.
- `create_video`: Whether to create a video. Default is `true`.
- `color`: The color of the contraction lines. Default is `(0.5, 0.5, 0.5, 0.5)`.
- `show_progress`: Whether to show progress information. Default is `false`.

# Returns
- If `create_gif` is `true`, returns the path to the generated GIF animation.
- If `create_video` is `true`, returns the path to the generated video.
"""
function OMEinsumContractionOrders.viz_contraction(
    ein::ET; 
    locs=StressLayout(),
    framerate = 30,
    filename = "contraction",
    pathname = ".",
    create_gif = false,
    create_video = true,
    color = (0.5, 0.5, 0.5, 0.5),
    show_progress::Bool = false
    ) where{ET <: Union{NestedEinsum, SlicedEinsum}}

    elimination_order = ein2elimination(ein)
    tng = TensorNetworkGraph(ein2hypergraph(ein))
    GViz = GraphViz(tng, locs)

    tempdirectory = mktempdir()
    # @info("Frames for animation \"$(filename)\" are being stored in directory: \n\t $(tempdirectory)")

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

    if create_gif
        Luxor.FFMPEG.exe(`-loglevel panic -r $(framerate) -f image2 -i $(tempdirectory)/%10d.png -filter_complex "[0:v] split [a][b]; [a] palettegen=stats_mode=full:reserve_transparent=on:transparency_color=FFFFFF [p]; [b][p] paletteuse=new=1:alpha_threshold=128" -y $(tempdirectory)/$(filename).gif`)

        if !isempty(pathname)
            if !isdir(pathname)
                @error "$pathname is not a directory."
            end
            fig_path = joinpath(pathname, "$filename.gif")
            mv("$(tempdirectory)/$(filename).gif", fig_path, force = true)
            @info("GIF is: $fig_path")
            giffn = fig_path
        else
            @info("GIF is: $(tempdirectory)/$(filename).gif")
            giffn = tempdirectory * "/" * filename * ".gif"
        end

        return giffn
    elseif create_video
        movieformat = ".mp4"

        if !isempty(pathname)
            if !isdir(pathname)
                @error "$pathname is not a directory."
            end
            pathname = joinpath(pathname, "$(filename)$(movieformat)")
        else
            pathname = joinpath("$(tempdirectory)", "$(filename)$(movieformat)")
        end

        @info "Creating video at: $pathname"
        FFMPEG.ffmpeg_exe(`
            -loglevel panic
            -r $(framerate) 
            -f image2 
            -i $(tempdirectory)/%10d.png 
            -c:v libx264 
            -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"
            -r $(framerate) 
            -pix_fmt yuv420p 
            -y $(pathname)`)

        return pathname
    else
        return tempdirectory
    end
end