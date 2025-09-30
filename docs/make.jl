using Porosity
using Documenter, DocumenterVitepress

include("pages.jl")
DocMeta.setdocmeta!(Porosity, :DocTestSetup, :(using Porosity); recursive = true)

makedocs(;
    # modules = [Porosity],
    # authors = "Abhinav Pratap Singh",
    # repo = "https://github.com/ayushinav/Porosity.jl", #.jl/blob/{commit}{path}#{line}",
    # sitename = "Porosity.jl",
    # doctest = true,
    clean = true,
    # format=Documenter.HTML(;
    #     prettyurls=get(ENV, "CI", "false") == "true",
    #     canonical="https://ayushinav.github.io/Porosity.jl",
    #     edit_link="main",
    #     assets=String[]
    # ),

    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/ayushinav/Porosity.jl",
        devurl = "dev", #, md_output_path = ".", 
        deploy_url="/"
        # md_output_path = "docs2"
        # build_vitepress = false
    ),
    pages = pages,
)

# DocumenterVitepress.deploydocs(; target = "../docs2", repo="github.com/ayushinav/Porosity.git", devbranch="main")
