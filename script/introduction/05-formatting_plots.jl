
using DifferentialEquations, Plots, ParameterizedFunctions
gr()
lorenz = @ode_def Lorenz begin
  dx = σ*(y-x)
  dy = ρ*x-y-x*z
  dz = x*y-β*z
end σ β ρ

p = [10.0,8/3,28]
u0 = [1., 5., 10.]
tspan = (0., 100.)
prob = ODEProblem(lorenz, u0, tspan, p)
sol = solve(prob)


plot(sol)


plot(sol,vars=(1, 2, 3))


plot(sol,vars=[:x])


plot(sol,vars=(1,2,3))
plot(sol,vars=[1])


plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
xaxis="Time (t)",yaxis="u(t) (in mm)",label=["X","Y","Z"])


scatter(sol,vars=[:x])


plot(sol,vars=(1,2,3),denseplot=false)


plot(sol,vars=(1,2,3),plotdensity=100)


plot(sol,vars=(1,2,3),plotdensity=10000)


plot(sol,vars=(1,2,3))
scatter!(sol,vars=(1,2,3),plotdensity=100)


p = plot(sol,vars=(1,2,3))
scatter!(p,sol,vars=(1,2,3),plotdensity=100)
title!("I added a title")


using DiffEqTutorials
DiffEqTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

