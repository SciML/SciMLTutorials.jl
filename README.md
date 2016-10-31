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

The notebooks can be viewed remotely on Github or via [nbviewer]()

- Introductions
  - Introduction to Solving ODEs
  - Introduction to Solving SDEs
  - Finite Element Solutions to the Poisson Equation
- Plotting
  - Formatting the Plots
  - Animating the Solution to the Stochastic Heat Equation
- Extra ODE Features
  - Faster Calculations via In-Place Operations
  - Defining ODEs using ParameterizedFunctions.jl
  - Solving Equations with Julia-Defined Types
  - Unit Check Arithmetic via Unitful.jl
  - Using General Tableau Methods
  - Feagin's Order 10, 12, and 14 Methods
  - Calling External Solvers - ODE.jl and ODEInterface
- Extra FEM Features
  - Heat Equation System Differing Diffusion Constants
  - Solving the Gierer-Meinhardt Equations
