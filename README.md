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

To generate the interactive notebooks, first install the SciMLTutorials, instantiate the
environment, and then run `SciMLTutorials.open_notebooks()`. This looks as follows:

```julia
]add SciMLTutorials#master
]activate SciMLTutorials
]instantiate
using SciMLTutorials
SciMLTutorials.open_notebooks()
```

The tutorials will be generated at your `pwd()` in a folder called `generated_notebooks`.

Note that when running the tutorials, the packages are not automatically added. Thus you
will need to add the packages manually or use the internal Project/Manifest tomls to
instantiate the correct packages. This can be done by activating the folder of the tutorials.
For example,

```julia
using Pkg
Pkg.activate(joinpath(pkgdir(SciMLTutorials),"tutorials","models"))
Pkg.instantiate()
```

will add all of the packages required to run any tutorial in the `models` folder.

## Contributing

All of the files are generated from the Weave.jl files in the `tutorials` folder. The generation process runs automatically,
and thus one does not necessarily need to test the Weave process locally. Instead, simply open a PR that adds/updates a
file in the "tutorials" folder and the PR will generate the tutorial on demand. Its artifacts can then be inspected in the
Buildkite as described below before merging. Note that it will use the Project.toml and Manifest.toml of the subfolder, so
any changes to dependencies requires that those are updated.

### Reporting Bugs and Issues

Report any bugs or issues at [the SciMLTutorials repository](https://github.com/SciML/SciMLTutorials.jl/issues).

### Inspecting Tutorial Results

To see tutorial results before merging, click into the BuildKite, click onto
Artifacts, and then investigate the trained results.

![](https://user-images.githubusercontent.com/1814174/118359358-02ddc980-b551-11eb-8a9b-24de947cefee.PNG)

### Manually Generating Files

To run the generation process, do for example:

```julia
]activate SciMLTutorials # Get all of the packages
using SciMLTutorials
SciMLTutorials.weave_file("models","01-classical_physics.jmd")
```

To generate all of the files in a folder, for example, run:

```julia
SciMLTutorials.weave_folder("models")
```

To generate all of the notebooks, do:

```julia
SciMLTutorials.weave_all()
```

Each of the tuturials displays the computer characteristics at the bottom of
the benchmark.
