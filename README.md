# DiffEqTutorials.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

DiffEqTutorials.jl holds PDFs, webpages, and interactive Jupyter notebooks
showing how to utilize the software in the JuliaDiffEq ecosystem. This set of
tutorials was made to complement the
[documentation](http://docs.juliadiffeq.org/latest/) and the
[devdocs](http://devdocs.juliadiffeq.org/latest/)
by providing practical examples of the concepts. For more details, please
consult the docs.

## Video Tutorial

[![Video Tutorial](https://user-images.githubusercontent.com/1814174/36342812-bdfd0606-13b8-11e8-9eff-ff219de909e5.PNG)](https://youtu.be/KPEqYtEd-zY)

## Table of Contents

- Introduction
  - [Introduction to DifferentialEquations.jl through ODEs](https://htmlpreview.github.io/?https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/html/ode_introduction.html)
  - [Detecting Stiffness and Choosing an ODE Algorithm](https://htmlpreview.github.io/?https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/html/choosing_algs.html)
  - [Optimizing your DiffEq Code](https://htmlpreview.github.io/?https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/html/optimizing_diffeq_code.html)
  - [Callbacks and Event Handling](https://htmlpreview.github.io/?https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/html/callbacks_and_events.html)
  - [Formatting Plots](https://htmlpreview.github.io/?https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/html/formatting_plots.html)
- Modeling Examples
  - [Classical Physics Models](https://htmlpreview.github.io/?https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/html/classical_physics.html)
  - [Conditional Dosing Example](https://htmlpreview.github.io/?https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/html/conditional_dosing.html)
  - [Kepler Problem Orbit](https://htmlpreview.github.io/?https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/html/kepler_problem.html)
- Advanced ODE Features
  - [Feagin's Order 10, 12, and 14 Methods](https://htmlpreview.github.io/?https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/html/feagin.html)
  - [Finding Maxima and Minima of DiffEq Solutions](https://htmlpreview.github.io/?https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/html/ode_minmax.html)
  - [Monte Carlo Parameter Estimation from Data](https://htmlpreview.github.io/?https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/html/monte_carlo_parameter_estim.html)
- Type Handling
  - [Solving Equations with Julia-Defined Types](https://htmlpreview.github.io/?https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/html/number_types.html)
  - [Numbers with Uncertainties](https://htmlpreview.github.io/?https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/html/uncertainties.html)
  - [Unit Check Arithmetic via Unitful.jl](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqTutorials.jl/blob/master/ExtraODEFeatures/Unit%20Checked%20Arithmetic%20via%20Unitful.ipynb)
- Advanced
  - [A 2D Cardiac Electrophysiology Model (CUDA-accelerated PDE solver)](https://htmlpreview.github.io/?https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/html/beeler_reuter.html)

## Contributing

All of the files are generated from the Weave.jl files in the `tutorials` folder. To run the generation process, do for example:

```julia
using DiffEqTutorials
DiffEqTutorials.weave_file(".","introduction/ode_introduction.jmd")
```
