module DiffEqTutorials

using Weave, Pkg, InteractiveUtils, IJulia

repo_directory = joinpath(@__DIR__,"..")

function weave_file(folder,file,build_list=(:script,:html,:pdf,:notebook))
  println("File: $file")
  tmp = joinpath(repo_directory,"tutorials",folder,file)
  args = Dict{Symbol,String}(:folder=>folder,:file=>file)
  if :script ∈ build_list
    println("Building Script")
    dir = joinpath(repo_directory,"script",folder)
    isdir(dir) || mkdir(dir)
    tangle(tmp;out_path=dir)
  end
  if :html ∈ build_list
    println("Building HTML")
    dir = joinpath(repo_directory,"html",folder)
    isdir(dir) || mkdir(dir)
    weave(tmp,doctype = "md2html",out_path=dir,args=args)
  end
  if :pdf ∈ build_list
    println("Building PDF")
    dir = joinpath(repo_directory,"pdf",folder)
    isdir(dir) || mkdir(dir)
    weave(tmp,doctype="md2pdf",out_path=dir,args=args)
  end
  if :github ∈ build_list
    println("Building Github Markdown")
    dir = joinpath(repo_directory,"markdown",folder)
    isdir(dir) || mkdir(dir)
    weave(tmp,doctype = "github",out_path=dir,args=args)
  end
  if :notebook ∈ build_list
    println("Building Notebook")
    dir = joinpath(repo_directory,"notebook",folder)
    isdir(dir) || mkdir(dir)
    Weave.notebook(tmp,dir)
  end
end

#=
# Needs two arg form
function weave_all()
  foreach(weave_file,
          file for file in readdir("tutorials") if endswith(file, ".jmd"))
end
=#

function tutorial_footer(folder=nothing,file=nothing)
  println("""
  These benchmarks are part of the DiffEqTutorials.jl repository, found at:

  https://github.com/JuliaDiffEq/DiffEqTutorials.jl
  """)
  if folder !== nothing && file !== nothing
    println("""
    To locally run this tutorial, do the following commands:

    using DiffEqTutorials
    DiffEqTutorials.weave_file("$folder","$file")
    """)
  end
  println("Computer Information:\n")
  InteractiveUtils.versioninfo()
  println()
  println("Package Information:\n")
  DiffEqTutorials.Pkg.status()
end

function open_notebooks()
  Base.eval(Main, Meta.parse("import IJulia"))
  path = joinpath(repo_directory,"notebook")
  IJulia.notebook(;dir=path)
end

end
