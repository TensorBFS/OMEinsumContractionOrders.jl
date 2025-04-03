using OMEinsumContractionOrders
using Documenter
using DocThemeIndigo

DocMeta.setdocmeta!(OMEinsumContractionOrders, :DocTestSetup, :(using OMEinsumContractionOrders); recursive=true)
indigo = DocThemeIndigo.install(OMEinsumContractionOrders)

makedocs(;
    modules=[OMEinsumContractionOrders],
    authors="GiggleLiu <cacate0129@gmail.com> and contributors",
    sitename="OMEinsumContractionOrders.jl",
    format=Documenter.HTML(;
        canonical="https://GiggleLiu.github.io/OMEinsumContractionOrders.jl",
        edit_link="main",
        assets=String[indigo],
    ),
    doctest = ("doctest=true" in ARGS),
    pages=[
        "Home" => "index.md",
        "Reference" => "ref.md",
    ],
)

deploydocs(;
    repo="github.com/GiggleLiu/OMEinsumContractionOrders.jl",
    devbranch="main",
)
