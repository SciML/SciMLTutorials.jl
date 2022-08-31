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

## Interactive Notebooks

To run the tutorials interactively via Jupyter notebooks, install the package
and open the tutorials like:

```julia
using Pkg
pkg"add https://github.com/SciML/SciMLTutorials.jl"
using SciMLTutorials
SciMLTutorials.open_notebooks()
```

## Video Tutorial

[![Video Tutorial](https://user-images.githubusercontent.com/1814174/36342812-bdfd0606-13b8-11e8-9eff-ff219de909e5.PNG)](https://youtu.be/KPEqYtEd-zY)

## Table of Contents

- Introduction
  - [Introduction to DifferentialEquations.jl through ODEs](http://tutorials.sciml.ai/html/introduction/01-ode_introduction.html)
  - [Detecting Stiffness and Choosing an ODE Algorithm](http://tutorials.sciml.ai/html/introduction/02-choosing_algs.html)
  - [Optimizing your DiffEq Code](http://tutorials.sciml.ai/html/introduction/03-optimizing_diffeq_code.html)
  - [Callbacks and Event Handling](http://tutorials.sciml.ai/html/introduction/04-callbacks_and_events.html)
  - [Formatting Plots](http://tutorials.sciml.ai/html/introduction/05-formatting_plots.html)
- Exercise Sheets
  - [DifferentialEquations.jl Workshop Exercises](http://tutorials.sciml.ai/html/exercises/01-workshop_exercises.html)
  - [DifferentialEquations.jl Workshop Exercise Solutions](http://tutorials.sciml.ai/html/exercises/02-workshop_solutions.html)
- Modeling Examples
  - [Classical Physics Models](http://tutorials.sciml.ai/html/models/01-classical_physics.html)
  - [Conditional Dosing Example](http://tutorials.sciml.ai/html/models/02-conditional_dosing.html)
  - [DiffEqBiological Tutorial I: Introduction](http://tutorials.sciml.ai/html/models/03-diffeqbio_I_introduction.html)
  - [DiffEqBiological Tutorial II: Network Properties API](http://tutorials.sciml.ai/html/models/04-diffeqbio_II_networkproperties.html)
  - [DiffEqBiological Tutorial III: Steady-States and Bifurcations](http://tutorials.sciml.ai/html/models/04b-diffeqbio_III_steadystates.html)
  - [Tutorial on using spatial SSAs in DiffEqJump](http://tutorials.sciml.ai/html/jumps/spatial.html)
  - [Kepler Problem Orbit](http://tutorials.sciml.ai/html/models/05-kepler_problem.html)
  - [Spiking Neural Systems](http://tutorials.sciml.ai/html/models/08-spiking_neural_systems.html)
- Advanced ODE Features
  - [Feagin's Order 10, 12, and 14 Methods](http://tutorials.sciml.ai/html/ode_extras/01-feagin.html)
  - [Finding Maxima and Minima of DiffEq Solutions](http://tutorials.sciml.ai/html/ode_extras/02-ode_minmax.html)
- Model Inference
  - [Bayesian Inference of Pendulum Parameters](http://tutorials.sciml.ai/html/model_inference/01-pendulum_bayesian_inference.html)
  - [Monte Carlo Parameter Estimation from Data](http://tutorials.sciml.ai/html/model_inference/02-monte_carlo_parameter_estim.html)
- Type Handling
  - [Solving Equations with Julia-Defined Types](http://tutorials.sciml.ai/html/type_handling/01-number_types.html)
  - [Numbers with Uncertainties](http://tutorials.sciml.ai/html/type_handling/02-uncertainties.html)
  - [Unit Check Arithmetic via Unitful.jl](http://tutorials.sciml.ai/html/type_handling/03-unitful.html)
- DiffEqUncertainty
  - [An Intro to Expectations via DiffEqUncertainty.jl](http://tutorials.sciml.ai/html/DiffEqUncertainty/01-expectation_introduction.html)
  - [Optimization Under Uncertainty with DiffEqUncertainty.jl](http://tutorials.sciml.ai/html/DiffEqUncertainty/02-AD_and_optimization.html)
  - [GPU-Accelerated Data-Driven Bayesian Uncertainty Quantification with Koopman Operators](http://tutorials.sciml.ai/html/DiffEqUncertainty/03-GPU_Bayesian_Koopman.html)
- Advanced
  - [A 2D Cardiac Electrophysiology Model (CUDA-accelerated PDE solver)](http://tutorials.sciml.ai/html/advanced/01-beeler_reuter.html)
  - [Solving Stiff Equations](http://tutorials.sciml.ai/html/advanced/02-advanced_ODE_solving.html)
  - [Solving the heat equation with diffusion-implicit time-stepping](http://tutorials.sciml.ai/html/advanced/04-diffusion_implicit_heat_equation.html)
  - [Kolmogorov Backward Equations](http://tutorials.sciml.ai/html/advanced/03-kolmogorov_equations.html)
- Perturbation Theory
  - [Mixed Symbolic/Numerical Methods for Perturbation Theory - Algebraic Equations](http://tutorials.sciml.ai/html/perturbation/01-perturbation_algebraic.html)
  - [Mixed Symbolic/Numerical Methods for Perturbation Theory - Differential Equations](http://tutorials.sciml.ai/html/perturbation/02-perturbation_differential.html)


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
