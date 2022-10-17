module SciMLTutorials

using Weave, Pkg, IJulia, InteractiveUtils, Markdown

repo_directory = joinpath(@__DIR__,"..")
cssfile = joinpath(@__DIR__, "..", "templates", "skeleton_css.css")
latexfile = joinpath(@__DIR__, "..", "templates", "julia_tex.tpl")
default_builds = (:script,:github)

function weave_file(folder,file,build_list=default_builds)
  target = joinpath(folder, file)
  @info("Weaving $(target)")

  if isfile(joinpath(folder, "Project.toml")) && build_list != (:notebook,)
    @info("Instantiating", folder)
    Pkg.activate(joinpath(folder))
    Pkg.instantiate()
    Pkg.build()

    @info("Printing out `Pkg.status()`")
    Pkg.status()
  end

  args = Dict{Symbol,String}(:folder=>folder,:file=>file)
  if :script ∈ build_list
    println("Building Script")
    dir = joinpath(repo_directory,"script",basename(folder))
    mkpath(dir)
    tangle(target; out_path=dir)
  end
  if :html ∈ build_list
    println("Building HTML")
    dir = joinpath(repo_directory,"html",basename(folder))
    mkpath(dir)
    weave(target,doctype = "md2html",out_path=dir,args=args,css=cssfile,fig_ext=".svg")
  end
  if :pdf ∈ build_list
    println("Building PDF")
    dir = joinpath(repo_directory,"pdf",basename(folder))
    mkpath(dir)
    try
      weave(target,doctype="md2pdf",out_path=dir,template=latexfile,args=args)
    catch ex
      @warn "PDF generation failed" exception=(ex, catch_backtrace())
    end
  end
  if :github ∈ build_list
    println("Building Github Markdown")
    dir = joinpath(repo_directory,"markdown",basename(folder))
    mkpath(dir)
    weave(target,doctype = "github",out_path=dir,args=args)
  end
  if :notebook ∈ build_list
    println("Building Notebook")
    dir = joinpath(repo_directory,"notebook",basename(folder))
    mkpath(dir)
    Weave.convert_doc(target,joinpath(dir,file[1:end-4]*".ipynb"))
  end
end

function weave_all(build_list=default_builds)
  for folder in readdir(joinpath(repo_directory,"tutorials"))
    folder == "test.jmd" && continue
    weave_folder(joinpath(repo_directory,"tutorials",folder),build_list)
  end
end

function weave_folder(folder,build_list=default_builds)
  for file in readdir(joinpath(folder))
    # Skip non-`.jmd` files
    if !endswith(file, ".jmd")
      continue
    end

    try
      weave_file(folder,file,build_list)
    catch e
      @error(e)
    end
  end
end

function tutorial_footer(folder=nothing, file=nothing)
    display(md"""
    ## Appendix

    These tutorials are a part of the SciMLTutorials.jl repository, found at: <https://github.com/SciML/SciMLTutorials.jl>.
    For more information on high-performance scientific machine learning, check out the SciML Open Source Software Organization <https://sciml.ai>.

    """)
    if folder !== nothing && file !== nothing
        display(Markdown.parse("""
        To locally run this tutorial, do the following commands:
        ```
        using SciMLTutorials
        SciMLTutorials.weave_file("$folder","$file")
        ```
        """))
    end
    display(md"Computer Information:")
    vinfo = sprint(InteractiveUtils.versioninfo)
    display(Markdown.parse("""
    ```
    $(vinfo)
    ```
    """))

    display(md"""
    Package Information:
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
    display(Markdown.parse(md))
end

function open_notebooks()
  Base.eval(Main, Meta.parse("import IJulia"))
  weave_all((:notebook,))
  path = joinpath(repo_directory,"notebook")
  newpath = joinpath(pwd(),"generated_notebooks")
  mv(path, newpath)
  IJulia.notebook(;dir=newpath)
end

end
