# DiffEqTutorials.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

DiffEqTutorials.jl holds Jupyter notebooks showing how to utilize the software in
the JuliaDiffEq ecosystem.

## Viewing the Notebooks Locally

To view the notebooks locally and interact with the contents, use the following
commands (requires [IJulia](https://github.com/JuliaLang/IJulia.jl)):

```julia
Pkg.clone("https://github.com/JuliaDiffEq/DiffEqTutorials.jl")
using IJulia
notebook(dir=Pkg.dir("DiffEqTutorials"))
```

## Table of Contents

The notebooks can be viewed remotely on Github or via [nbviewer](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqTutorials.jl/tree/master/)

- Plotting
  - [Formatting the Plots](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqTutorials.jl/blob/master/Plotting/Formatting%20the%20Plots.ipynb)
- Extra ODE Features
  - [Faster Calculations via In-Place Operations](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqTutorials.jl/blob/master/ExtraODEFeatures/Faster%20Calculations%20via%20In-Place%20Operations.ipynb)
  - [Solving Equations with Julia-Defined Types](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqTutorials.jl/blob/master/ExtraODEFeatures/Solving%20Equations%20in%20With%20Julia-Defined%20Types.ipynb)
  - [Unit Check Arithmetic via Unitful.jl](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqTutorials.jl/blob/master/ExtraODEFeatures/Unit%20Checked%20Arithmetic%20via%20Unitful.ipynb)
  - [Using General Tableau Methods](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqTutorials.jl/blob/master/ExtraODEFeatures/Using%20General%20Tableau%20Methods.ipynb)
  - [Feagin's Order 10, 12, and 14 Methods](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqTutorials.jl/blob/master/ExtraODEFeatures/Feagin%27s%20Order%2010%2C%2012%2C%20and%2014%20methods.ipynb)
- Extra FEM Features
  - [Heat Equation System Differing Diffusion Constants](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqTutorials.jl/blob/master/ExtraFEMFeatures/Heat%20Equation%20System%20Differing%20Diffusion%20Constants.ipynb)
  - [Solving the Gierer-Meinhardt Equations](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqTutorials.jl/blob/master/ExtraFEMFeatures/Solving%20the%20Gierer-Meinhardt%20Equations.ipynb)
