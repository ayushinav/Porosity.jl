# using Pkg
# Pkg.activate("docs/.")

using Porosity
using Documenter, DocumenterVitepress, DocumenterMermaid

include("pages.jl")
DocMeta.setdocmeta!(Porosity, :DocTestSetup, :(using Porosity); recursive = true)

makedocs(;
    modules = [Porosity],
    authors = "Abhinav Pratap Singh",
    repo = "https://github.com/ayushinav/Porosity.jl", #.jl/blob/{commit}{path}#{line}",
    sitename = "Porosity.jl",
    doctest = true,
    # format=Documenter.HTML(;
    #     prettyurls=get(ENV, "CI", "false") == "true",
    #     canonical="https://ayushinav.github.io/Porosity.jl",
    #     edit_link="main",
    #     assets=String[]
    # ),

    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/ayushinav/Porosity.jl",
    ),
    pages = pages,
)

# deploydocs(; repo="github.com/ayushinav/Porosity.git", devbranch="main")
