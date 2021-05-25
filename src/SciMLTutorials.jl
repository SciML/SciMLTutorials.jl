module SciMLTutorials

using Weave, Pkg, InteractiveUtils, IJulia

repo_directory = joinpath(@__DIR__,"..")
cssfile = joinpath(@__DIR__, "..", "templates", "skeleton_css.css")
latexfile = joinpath(@__DIR__, "..", "templates", "julia_tex.tpl")

function weave_file(folder,file,build_list=(:script,:html,:pdf,:github,:notebook); kwargs...)
  tmp = joinpath(repo_directory,"tutorials",folder,file)
  Pkg.activate(dirname(tmp))
  Pkg.instantiate()
  args = Dict{Symbol,String}(:folder=>folder,:file=>file)
  if :script ∈ build_list
    println("Building Script")
    dir = joinpath(repo_directory,"script",folder)
    mkpath(dir)
    args[:doctype] = "script"
    tangle(tmp;out_path=dir)
  end
  if :html ∈ build_list
    println("Building HTML")
    dir = joinpath(repo_directory,"html",folder)
    mkpath(dir)
    args[:doctype] = "html"
    weave(tmp,doctype = "md2html",out_path=dir,args=args; fig_ext=".svg", css=cssfile, kwargs...)
  end
  if :pdf ∈ build_list
    println("Building PDF")
    dir = joinpath(repo_directory,"pdf",folder)
    mkpath(dir)
    args[:doctype] = "pdf"
    try
      weave(tmp,doctype="md2pdf",out_path=dir,args=args; template=latexfile, kwargs...)
    catch ex
      @warn "PDF generation failed" exception=(ex, catch_backtrace())
    end
  end
  if :github ∈ build_list
    println("Building Github Markdown")
    dir = joinpath(repo_directory,"markdown",folder)
    mkpath(dir)
    args[:doctype] = "github"
    weave(tmp,doctype = "github",out_path=dir,args=args; kwargs...)
  end
  if :notebook ∈ build_list
    println("Building Notebook")
    dir = joinpath(repo_directory,"notebook",folder)
    mkpath(dir)
    args[:doctype] = "notebook"
    Weave.convert_doc(tmp,joinpath(dir,file[1:end-4]*".ipynb"))
  end
end

function weave_all()
  for folder in readdir(joinpath(repo_directory,"tutorials"))
    folder == "test.jmd" && continue
    weave_folder(folder)
  end
end

function weave_folder(folder)
  for file in readdir(joinpath(repo_directory,"tutorials",folder))
    println("Building $(joinpath(folder,file)))")
    try
      weave_file(folder,file)
    catch
    end
  end
end

function tutorial_footer(folder=nothing, file=nothing; remove_homedir=true)
    display("text/markdown", """
    ## Appendix

     This tutorial is part of the SciMLTutorials.jl repository, found at: <https://github.com/SciML/SciMLTutorials.jl>.
     For more information on doing scientific machine learning (SciML) with open source software, check out <https://sciml.ai/>.
    """)
    if folder !== nothing && file !== nothing
        display("text/markdown", """
        To locally run this tutorial, do the following commands:
        ```
        using SciMLTutorials
        SciMLTutorials.weave_file("$folder","$file")
        ```
        """)
    end
    display("text/markdown", "Computer Information:")
    vinfo = sprint(InteractiveUtils.versioninfo)
    display("text/markdown",  """
    ```
    $(vinfo)
    ```
    """)
    proj = sprint(io -> Pkg.status(io=io))
    mani = sprint(io -> Pkg.status(io=io, mode = Pkg.PKGMODE_MANIFEST))

    md = """
    ```
    $(chomp(proj))
    ```

    And the full manifest:

    ```
    $(chomp(mani))
    ```
    """
    display("text/markdown", md)
end

function open_notebooks()
  Base.eval(Main, Meta.parse("import IJulia"))
  path = joinpath(repo_directory,"notebook")
  IJulia.notebook(;dir=path)
end

end
