using TensorQEC, OMEinsumContractionOrders

function main(optimizer)
    @info "Running QEC with optimizer: $(optimizer)"
    for d in [9, 13, 17, 21]
        @info "Running surface code for d = $d"
        tanner = CSSTannerGraph(SurfaceCode(d, d))
        time_elapsed = @elapsed ct = compile(TNMAP(; optimizer),tanner)

        # NOTE: TreeSA gives sc = 16, tc=21.74 for d = 9
        @info "Contraction complexity: $(contraction_complexity(ct.cd.net)), time cost: $(time_elapsed)s"
        return contraction_complexity(ct.cd.net)
    end
end