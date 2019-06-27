
using Measurements

5.23 ± 0.14 === 5.23 ± 0.14


(5.23± 0.14) - (5.23 ± 0.14)


(5.23 ± 0.14) / (5.23 ± 0.14)


x = 5.23 ± 0.14
x === x


x - x


x / x


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


println("Quantity of carbon-14 after ",  sol.t[11], " years:")
println("Numerical: ", sol[11])
println("Analytic:  ", u[11])


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


plot(sol.t, getindex.(sol.u, 2) .- u, label = "")


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


using DiffEqTutorials
DiffEqTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

