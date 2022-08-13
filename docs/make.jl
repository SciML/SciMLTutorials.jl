using Documenter, SciMLTutorialsOutput

dir = @__DIR__() * "/.."

@show dir
@show readdir(dir)

include("pages.jl")

makedocs(
    sitename="The SciML Tutorials",
    authors="Chris Rackauckas",
    modules=[SciMLTutorialsOutput],
    clean=true, doctest=false,
    format=Documenter.HTML(#analytics = "UA-90474609-3",
        assets=["assets/favicon.ico"],
        canonical="https://tutorials.sciml.ai/stable/"),
    pages=pages
)

deploydocs(;
    repo="github.com/SciML/SciMLTutorialsOutput",
    devbranch="main",
    branch="main"
)
