using Porosity
using Documenter, DocumenterVitepress

include("pages.jl")
DocMeta.setdocmeta!(Porosity, :DocTestSetup, :(using Porosity); recursive = true)

makedocs(;
    authors = "Abhinav Pratap Singh et al.",
    sitename = "Porosity.jl",
    doctest = true,
    clean = true,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/ayushinav/Porosity.jl",
        devurl = "dev",
        devbranch = "main",
    ),
    pages = pages,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/ayushinav/Porosity.jl",
    devbranch = "main",
    push_preview = true,
)
