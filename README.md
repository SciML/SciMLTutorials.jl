# DiffEqTutorials.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

DiffEqTutorials.jl holds PDFs, webpages, and interactive Jupyter notebooks
showing how to utilize the software in the JuliaDiffEq ecosystem. This set of
tutorials was made to complement the
[documentation](http://docs.juliadiffeq.org/latest/) and the
[devdocs](http://devdocs.juliadiffeq.org/latest/)
by providing practical examples of the concepts. For more details, please
consult the docs.

## Interactive Notebooks

To run the tutorials interactively via Jupyter notebooks, install the package
and open the tutorials like:

```julia
using Pkg
pkg"add https://github.com/JuliaDiffEq/DiffEqTutorials.jl"
using DiffEqTutorials
DiffEqTutorials.open_notebooks()
```

## Video Tutorial

[![Video Tutorial](https://user-images.githubusercontent.com/1814174/36342812-bdfd0606-13b8-11e8-9eff-ff219de909e5.PNG)](https://youtu.be/KPEqYtEd-zY)

## Table of Contents

- Introduction
  - [Introduction to DifferentialEquations.jl through ODEs](http://juliadiffeq.org/DiffEqTutorials.jl/html/introduction/ode_introduction.html)
  - [Detecting Stiffness and Choosing an ODE Algorithm](http://juliadiffeq.org/DiffEqTutorials.jl/html/introduction/choosing_algs.html)
  - [Optimizing your DiffEq Code](http://juliadiffeq.org/DiffEqTutorials.jl/html/introduction/optimizing_diffeq_code.html)
  - [Callbacks and Event Handling](http://juliadiffeq.org/DiffEqTutorials.jl/html/introduction/callbacks_and_events.html)
  - [Formatting Plots](http://juliadiffeq.org/DiffEqTutorials.jl/html/introduction/formatting_plots.html)
- Modeling Examples
  - [Classical Physics Models](http://juliadiffeq.org/DiffEqTutorials.jl/html/models/classical_physics.html)
  - [Conditional Dosing Example](http://juliadiffeq.org/DiffEqTutorials.jl/html/models/conditional_dosing.html)
  - [Kepler Problem Orbit](http://juliadiffeq.org/DiffEqTutorials.jl/html/models/kepler_problem.html)
- Advanced ODE Features
  - [Feagin's Order 10, 12, and 14 Methods](http://juliadiffeq.org/DiffEqTutorials.jl/html/ode_extras/feagin.html)
  - [Finding Maxima and Minima of DiffEq Solutions](http://juliadiffeq.org/DiffEqTutorials.jl/html/ode_extras/ode_minmax.html)
  - [Monte Carlo Parameter Estimation from Data](http://juliadiffeq.org/DiffEqTutorials.jl/html/ode_extras/monte_carlo_parameter_estim.html)
- Type Handling
  - [Solving Equations with Julia-Defined Types](http://juliadiffeq.org/DiffEqTutorials.jl/html/type_handling/number_types.html)
  - [Numbers with Uncertainties](http://juliadiffeq.org/DiffEqTutorials.jl/html/type_handling/uncertainties.html)
  - [Unit Check Arithmetic via Unitful.jl](http://juliadiffeq.org/DiffEqTutorials.jl/html/type_handling/unitful.html)
- Advanced
  - [A 2D Cardiac Electrophysiology Model (CUDA-accelerated PDE solver)](http://juliadiffeq.org/DiffEqTutorials.jl/html/advanced/beeler_reuter.html)

## Contributing

First of all, make sure that your current directory is `DiffEqTutorials`. All
of the files are generated from the Weave.jl files in the `tutorials` folder.
To run the generation process, do for example:

```julia
using Pkg, DiffEqTutorials
cd(joinpath(dirname(pathof(DiffEqTutorials)), ".."))
Pkg.pkg"activate ."
Pkg.pkg"instantiate"
DiffEqTutorials.weave_file("introduction","ode_introduction.jmd")
```

To generate all of the notebooks, do:

```julia
DiffeqTutorials.weave_all()
```

If you add new tutorials which require new packages, simply updating your local
environment will change the project and manifest files. When this occurs, the
updated environment files should be included in the PR.
