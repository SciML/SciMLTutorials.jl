---
author: "Mosè Giordano, Chris Rackauckas"
title: "Numbers with Uncertainties"
---


The result of a measurement should be given as a number with an attached uncertainties, besides the physical unit, and all operations performed involving the result of the measurement should propagate the uncertainty, taking care of correlation between quantities.

There is a Julia package for dealing with numbers with uncertainties: [`Measurements.jl`](https://github.com/JuliaPhysics/Measurements.jl).  Thanks to Julia's features, `DifferentialEquations.jl` easily works together with `Measurements.jl` out-of-the-box.

This notebook will cover some of the examples from the tutorial about classical Physics.

## Caveat about `Measurement` type

Before going on with the tutorial, we must point up a subtlety of `Measurements.jl` that you should be aware of:

````julia

using Measurements

5.23 ± 0.14 === 5.23 ± 0.14
````


````
false
````



````julia

(5.23± 0.14) - (5.23 ± 0.14)
````


````
0.0 ± 0.2
````



````julia

(5.23 ± 0.14) / (5.23 ± 0.14)
````


````
1.0 ± 0.038
````





The two numbers above, even though have the same nominal value and the same uncertainties, are actually two different measurements that only by chance share the same figures and their difference and their ratio have a non-zero uncertainty.  It is common in physics to get very similar, or even equal, results for a repeated measurement, but the two measurements are not the same thing.

Instead, if you have *one measurement* and want to perform some operations involving it, you have to assign it to a variable:

````julia

x = 5.23 ± 0.14
x === x
````


````
true
````



````julia

x - x
````


````
0.0 ± 0.0
````



````julia

x / x
````


````
1.0 ± 0.0
````





## Radioactive Decay of Carbon-14

The rate of decay of carbon-14 is governed by a first order linear ordinary differential equation

$$\frac{\mathrm{d}u(t)}{\mathrm{d}t} = -\frac{u(t)}{\tau}$$

where $\tau$ is the mean lifetime of carbon-14, which is related to the half-life $t_{1/2} = (5730 \pm 40)$ years by the relation $\tau = t_{1/2}/\ln(2)$.

````julia

using DifferentialEquations, Measurements, Plots

# Half-life and mean lifetime of radiocarbon, in years
t_12 = 5730 ± 40
τ = t_12 / log(2)

#Setup
u₀ = 1 ± 0
tspan = (0.0, 10000.0)

#Define the problem
radioactivedecay(u,p,t) = - u / τ

#Pass to solver
prob = ODEProblem(radioactivedecay, u₀, tspan)
sol = solve(prob, Tsit5(), reltol = 1e-8)

# Analytic solution
u = exp.(- sol.t / τ)

plot(sol.t, sol.u, label = "Numerical", xlabel = "Years", ylabel = "Fraction of Carbon-14")
plot!(sol.t, u, label = "Analytic")
````


![](figures/02-uncertainties_7_1.png)



The two curves are perfectly superimposed, indicating that the numerical solution matches the analytic one.  We can check that also the uncertainties are correctly propagated in the numerical solution:

````julia

println("Quantity of carbon-14 after ",  sol.t[11], " years:")
println("Numerical: ", sol[11])
println("Analytic:  ", u[11])
````


````
Quantity of carbon-14 after 5207.541347908455 years:
Numerical: 0.5326 ± 0.0023
Analytic:  0.5326 ± 0.0023
````





Both the value of the numerical solution and its uncertainty match the analytic solution within the requested tolerance.  We can also note that close to 5730 years after the beginning of the decay (half-life of the radioisotope), the fraction of carbon-14 that survived is about 0.5.

## Simple pendulum

### Small angles approximation

The next problem we are going to study is the simple pendulum in the approximation of small angles.  We address this simplified case because there exists an easy analytic solution to compare.

The differential equation we want to solve is

$$\ddot{\theta} + \frac{g}{L} \theta = 0$$

where $g = (9.79 \pm 0.02)~\mathrm{m}/\mathrm{s}^2$ is the gravitational acceleration measured where the experiment is carried out, and $L = (1.00 \pm 0.01)~\mathrm{m}$ is the length of the pendulum.

When you set up the problem for `DifferentialEquations.jl` remember to define the measurements as variables, as seen above.

````julia

using DifferentialEquations, Measurements, Plots

g = 9.79 ± 0.02; # Gravitational constants
L = 1.00 ± 0.01; # Length of the pendulum

#Initial Conditions
u₀ = [0 ± 0, π / 60 ± 0.01] # Initial speed and initial angle
tspan = (0.0, 6.3)

#Define the problem
function simplependulum(du,u,p,t)
    θ  = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g/L)*θ
end

#Pass to solvers
prob = ODEProblem(simplependulum, u₀, tspan)
sol = solve(prob, Tsit5(), reltol = 1e-6)

# Analytic solution
u = u₀[2] .* cos.(sqrt(g / L) .* sol.t)

plot(sol.t, getindex.(sol.u, 2), label = "Numerical")
plot!(sol.t, u, label = "Analytic")
````


![](figures/02-uncertainties_9_1.png)



Also in this case there is a perfect superimposition between the two curves, including their uncertainties.

We can also have a look at the difference between the two solutions:

````julia

plot(sol.t, getindex.(sol.u, 2) .- u, label = "")
````


![](figures/02-uncertainties_10_1.png)



## Arbitrary amplitude

Now that we know how to solve differential equations involving numbers with uncertainties we can solve the simple pendulum problem without any approximation.  This time the differential equation to solve is the following:

$$\ddot{\theta} + \frac{g}{L} \sin(\theta) = 0$$

````julia

g = 9.79 ± 0.02; # Gravitational constants
L = 1.00 ± 0.01; # Length of the pendulum

#Initial Conditions
u₀ = [0 ± 0, π / 3 ± 0.02] # Initial speed and initial angle
tspan = (0.0, 6.3)

#Define the problem
function simplependulum(du,u,p,t)
    θ  = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g/L) * sin(θ)
end

#Pass to solvers
prob = ODEProblem(simplependulum, u₀, tspan)
sol = solve(prob, Tsit5(), reltol = 1e-6)

plot(sol.t, getindex.(sol.u, 2), label = "Numerical")
````


![](figures/02-uncertainties_11_1.png)



We note that in this case the period of the oscillations is not constant.


## Appendix

 This tutorial is part of the SciMLTutorials.jl repository, found at: <https://github.com/SciML/SciMLTutorials.jl>.
 For more information on doing scientific machine learning (SciML) with open source software, check out <https://sciml.ai/>.

To locally run this tutorial, do the following commands:
```
using SciMLTutorials
SciMLTutorials.weave_file("type_handling","02-uncertainties.jmd")
```

Computer Information:
```
Julia Version 1.4.2
Commit 44fa15b150* (2020-05-23 18:35 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Core(TM) i7-9700K CPU @ 3.60GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-8.0.1 (ORCJIT, skylake)
Environment:
  JULIA_LOAD_PATH = /builds/JuliaGPU/DiffEqTutorials.jl:
  JULIA_DEPOT_PATH = /builds/JuliaGPU/DiffEqTutorials.jl/.julia
  JULIA_CUDA_MEMORY_LIMIT = 2147483648
  JULIA_NUM_THREADS = 8

```

Package Information:

```
Status `/builds/JuliaGPU/DiffEqTutorials.jl/tutorials/type_handling/Project.toml`
[7e558dbc-694d-5a72-987c-6f4ebed21442] ArbNumerics 1.2.1
[55939f99-70c6-5e9b-8bb0-5071ed7d61fd] DecFP 1.0.0
[abce61dc-4473-55a0-ba07-351d65e31d42] Decimals 0.4.1
[0c46a032-eb83-5123-abaf-570d42b7fbaa] DifferentialEquations 6.15.0
[497a8b3b-efae-58df-a0af-a86822472b78] DoubleFloats 1.1.13
[eff96d63-e80a-5855-80a2-b1b0885c5ab7] Measurements 2.2.1
[1dea7af3-3e70-54e6-95c3-0bf5283fa5ed] OrdinaryDiffEq 5.42.5
[91a5bcdd-55d7-5caf-9e0b-520d859cae80] Plots 1.6.1
[1986cc42-f94f-5a68-af5c-568840ba703d] Unitful 1.4.0
```
