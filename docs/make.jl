using Pkg
Pkg.activate("../MT.jl/docs/.")

using MT
using Documenter

include("pages.jl")
DocMeta.setdocmeta!(MT, :DocTestSetup, :(using MT); recursive=true)

makedocs(;
    modules=[MT],
    authors="Abhinav Pratap Singh",
    repo="https://github.com/ayushinav/MT.jl", #.jl/blob/{commit}{path}#{line}",
    sitename="MT.jl",
    doctest=true,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ayushinav.github.io/MT.jl",
        edit_link="main",
        assets=String[]
    ),
    pages=pages
)

deploydocs(; repo="github.com/ayushinav/MT.git", devbranch="main")
