using Porosity
using Documenter, DocumenterVitepress

include("pages.jl")
DocMeta.setdocmeta!(Porosity, :DocTestSetup, :(using Porosity); recursive = true)

makedocs(;
    # modules = [Porosity],
    authors = "Abhinav Pratap Singh et al.",
    repo = "https://github.com/ayushinav/Porosity.jl",
    sitename = "Porosity.jl",
    doctest = true,
    clean = true,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/ayushinav/Porosity.jl",
        devurl = "dev",
        devbranch = "main",
        deploy_url = "/",
    ),
    pages = pages,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/ayushinav/Porosity.git",
    devbranch = "main",
    push_preview = true,
)
