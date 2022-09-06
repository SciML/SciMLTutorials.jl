# SciMLTutorials.jl: Tutorials for Scientific Machine Learning and Differential Equations

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](http://tutorials.sciml.ai/stable/)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/dev/highlevels/learning_resources/#SciMLTutorials)

[![Build status](https://badge.buildkite.com/8a39c2e1b44511eb84bdcd9019663cad757ae2479abd340508.svg)](https://buildkite.com/julialang/scimltutorials-dot-jl)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

SciMLTutorials.jl holds PDFs, webpages, and interactive Jupyter notebooks
showing how to utilize the software in the [SciML Scientific Machine Learning ecosystem](https://sciml.ai/).
This set of tutorials was made to complement the [documentation](https://sciml.ai/documentation/)
and the [devdocs](http://devdocs.sciml.ai/latest/)
by providing practical examples of the concepts. For more details, please
consult the docs.

## Results

To view the SciML Tutorials, go to [tutorials.sciml.ai](https://tutorials.sciml.ai/stable/). By default, this
will lead to the latest tagged version of the tutorials. To see the in-development version of the tutorials, go to
[https://tutorials.sciml.ai/dev/](https://tutorials.sciml.ai/dev/).

Static outputs in pdf, markdown, and html reside in [SciMLTutorialsOutput](https://github.com/SciML/SciMLTutorialsOutput).

## Video Tutorial

[![Video Tutorial](https://user-images.githubusercontent.com/1814174/36342812-bdfd0606-13b8-11e8-9eff-ff219de909e5.PNG)](https://youtu.be/KPEqYtEd-zY)

## Interactive Notebooks

To run the tutorials interactively via Jupyter notebooks, install the package
and open the tutorials like:

```julia
using Pkg
pkg"add https://github.com/SciML/SciMLTutorials.jl"
using SciMLTutorials
SciMLTutorials.open_notebooks()
```

## Contributing

First of all, make sure that your current directory is `SciMLTutorials`. All
of the files are generated from the Weave.jl files in the `tutorials` folder.
To run the generation process, do for example:

```julia
using Pkg, SciMLTutorials
cd(joinpath(dirname(pathof(SciMLTutorials)), ".."))
Pkg.pkg"activate ."
Pkg.pkg"instantiate"
SciMLTutorials.weave_file("introduction","01-ode_introduction.jmd")
```

To generate all of the notebooks, do:

```julia
SciMLTutorials.weave_all()
```

If you add new tutorials which require new packages, simply updating your local
environment will change the project and manifest files. When this occurs, the
updated environment files should be included in the PR.
