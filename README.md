# DiffEqTutorials.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

DiffEqTutorials.jl holds PDFs, webpages, and interactive Jupyter notebooks
showing how to utilize the software in the JuliaDiffEq ecosystem. This set of
tutorials was made to complement the
[documentation](http://docs.juliadiffeq.org/dev/) and the
[devdocs](http://devdocs.juliadiffeq.org/dev/)
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
  - [Introduction to DifferentialEquations.jl through ODEs](http://tutorials.juliadiffeq.org/html/introduction/01-ode_introduction.html)
  - [Detecting Stiffness and Choosing an ODE Algorithm](http://tutorials.juliadiffeq.org/html/introduction/02-choosing_algs.html)
  - [Optimizing your DiffEq Code](http://tutorials.juliadiffeq.org/html/introduction/03-optimizing_diffeq_code.html)
  - [Callbacks and Event Handling](http://tutorials.juliadiffeq.org/html/introduction/04-callbacks_and_events.html)
  - [Formatting Plots](http://tutorials.juliadiffeq.org/html/introduction/05-formatting_plots.html)
- Exercise Sheets
  - [DifferentialEquations.jl Workshop Exercises](http://tutorials.juliadiffeq.org/html/exercises/01-workshop_exercises.html)
  - [DifferentialEquations.jl Workshop Exercise Solutions](http://tutorials.juliadiffeq.org/html/exercises/02-workshop_solutions.html)
- Modeling Examples
  - [Classical Physics Models](http://tutorials.juliadiffeq.org/html/models/01-classical_physics.html)
  - [Conditional Dosing Example](http://tutorials.juliadiffeq.org/html/models/02-conditional_dosing.html)
  - [DiffEqBiological Tutorial I: Introduction](http://tutorials.juliadiffeq.org/html/models/03-diffeqbio_I_introduction.html)
  - [DiffEqBiological Tutorial II: Network Properties API](http://tutorials.juliadiffeq.org/html/models/04-diffeqbio_II_networkproperties.html)
  - [DiffEqBiological Tutorial III: Steady-States and Bifurcations](http://tutorials.juliadiffeq.org/html/models/04b-diffeqbio_III_steadystates.html)
  - [Kepler Problem Orbit](http://tutorials.juliadiffeq.org/html/models/05-kepler_problem.html)
  - [Kolmogorov Backward Equations](http://tutorials.juliadiffeq.org/html/models/08-kolmogorov_equations.html)
- Advanced ODE Features
  - [ModelingToolkit.jl, An IR and Compiler for Scientific Models](http://tutorials.juliadiffeq.org/html/ode_extras/01-ModelingToolkit.html)
  - [Feagin's Order 10, 12, and 14 Methods](http://tutorials.juliadiffeq.org/html/ode_extras/02-feagin.html)
  - [Finding Maxima and Minima of DiffEq Solutions](http://tutorials.juliadiffeq.org/html/ode_extras/03-ode_minmax.html)

- Model Inference
  - [Bayesian Inference of Pendulum Parameters](http://tutorials.juliadiffeq.org/html/model_inference/01-pendulum_bayesian_inference.html)
  - [Monte Carlo Parameter Estimation from Data](http://tutorials.juliadiffeq.org/html/model_inference/02-monte_carlo_parameter_estim.html)
- Type Handling
  - [Solving Equations with Julia-Defined Types](http://tutorials.juliadiffeq.org/html/type_handling/01-number_types.html)
  - [Numbers with Uncertainties](http://tutorials.juliadiffeq.org/html/type_handling/02-uncertainties.html)
  - [Unit Check Arithmetic via Unitful.jl](http://tutorials.juliadiffeq.org/html/type_handling/03-unitful.html)
- Advanced
  - [A 2D Cardiac Electrophysiology Model (CUDA-accelerated PDE solver)](http://tutorials.juliadiffeq.org/html/advanced/01-beeler_reuter.html)
  - [Solving Stiff Equations](http://tutorials.juliadiffeq.org/html/advanced/02-advanced_ODE_solving.html)

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
DiffEqTutorials.weave_all()
```

If you add new tutorials which require new packages, simply updating your local
environment will change the project and manifest files. When this occurs, the
updated environment files should be included in the PR.
