using Documenter, SciMLTutorialsOutput

dir = @__DIR__() * "/.."

@show dir
@show readdir(dir)

include("pages.jl")

mathengine = MathJax3(Dict(:loader => Dict("load" => ["[tex]/require", "[tex]/mathtools"]),
                           :tex => Dict("inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
                                        "packages" => [
                                            "base",
                                            "ams",
                                            "autoload",
                                            "mathtools",
                                            "require",
                                        ])))

makedocs(
    sitename="The SciML Tutorials",
    authors="Chris Rackauckas",
    modules=[SciMLTutorialsOutput],
    clean=true, doctest=false,
    format=Documenter.HTML(#analytics = "UA-90474609-3",
        assets=["assets/favicon.ico"],
        canonical="https://tutorials.sciml.ai/stable/",
        mathengine = mathengine),
    pages=pages
)

deploydocs(;
    repo="github.com/SciML/SciMLTutorialsOutput",
    devbranch="main",
    branch="main"
)
