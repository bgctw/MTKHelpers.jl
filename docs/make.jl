using MTKHelpers
using CairoMakie # to load cairomakie functions
using Documenter

DocMeta.setdocmeta!(MTKHelpers, :DocTestSetup, :(using MTKHelpers); recursive=true)

makedocs(;
    #modules=[MTKHelpers], # do not check for not included docstrings
    authors="Thomas Wutzler <twutz@bgc-jena.mpg.de> and contributors",
    repo="https://github.com/bgctw/MTKHelpers.jl/blob/{commit}{path}#{line}",
    sitename="MTKHelpers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://bgctw.github.io/MTKHelpers.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Update parameters" => "problemparsetter.md",
        "Solution handling" => "solution.md",
        "Embedding a system" => "embed_system.md",
        "Smooth steps" => "smoothstep.md",
        "CairoMakie Helpers" => "cairomakie.md",
        "index" => "zindex.md",
    ],
)

deploydocs(;
    repo="github.com/bgctw/MTKHelpers.jl",
    devbranch="main",
)
