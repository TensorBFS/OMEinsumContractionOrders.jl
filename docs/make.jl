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
        edit_link="master",
        assets=String[indigo],
    ),
    doctest = ("doctest=true" in ARGS),
    pages=[
        "Home" => "index.md",
        "Background Knowledge" => "background.md",
        "Tutorial" => "tutorial.md",
        "Choosing Optimizers" => "optimizers.md",
        "Reference" => "ref.md",
    ],
)

deploydocs(;
    repo="github.com/TensorBFS/OMEinsumContractionOrders.jl",
    devbranch="master",
)
