using MTKHelpers
using CairoMakie # to load cairomakie functions
using Documenter

# allow plot to work without display
# https://discourse.julialang.org/t/generation-of-documentation-fails-qt-qpa-xcb-could-not-connect-to-display/60988/2
ENV["GKSwstype"] = "100"

DocMeta.setdocmeta!(MTKHelpers, :DocTestSetup, :(using MTKHelpers); recursive = true)

makedocs(;
    #modules=[MTKHelpers], # do not check for not included docstrings
    authors = "Thomas Wutzler <twutz@bgc-jena.mpg.de> and contributors",
    #repo="https://github.com/bgctw/MTKHelpers.jl/blob/{commit}{path}#{line}",
    repo = Remotes.GitHub("bgctw", "MTKHelpers.jl"),
    sitename = "MTKHelpers.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://bgctw.github.io/MTKHelpers.jl",
        assets = String[]),
    pages = [
        "Home" => "index.md",
        "Update parameters" => [
            "Updating problems" => "updating_problems.md",
            "ProblemUpdater" => "problemupdater.md",
            "ODEProblem" => "odeproblemparsetter.md",
            "Type inference" => "concrete_parupdater.md",
            "PDE support" => "pde.md",
        ],
        "Translating symbols and Nums" => "system_num_dict.md",
        "Embedding a system" => "embed_system.md",
        "Solution handling" => "solution.md",
        "Smooth steps" => "smoothstep.md",
        "CairoMakie Helpers" => "cairomakie.md",
        "Developer notes" => [
            "PDE support" => "pde_dev.md",
        ],
        "index" => "zindex.md",
    ])

deploydocs(; repo = "github.com/bgctw/MTKHelpers.jl", devbranch = "main")
