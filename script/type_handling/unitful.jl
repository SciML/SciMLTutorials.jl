
using Unitful
t = 1.0u"s"


t2 = 1.02u"s"
t+t2


t*t2


sqrt(t)


t + sqrt(t)


using DifferentialEquations
f = (y,p,t) -> 0.5*y
u0 = 1.5u"N"
prob = ODEProblem(f,u0,(0.0u"s",1.0u"s"))
sol = solve(prob,Tsit5())


f = (y,p,t) -> 0.5*y/3.0u"s"
prob = ODEProblem(f,u0,(0.0u"s",1.0u"s"))
sol = solve(prob,Tsit5())


print(sol[:])


#Pkg.clone("https://github.com/ajkeller34/UnitfulPlots.jl")
using UnitfulPlots, Plots
gr()
plot(sol.t,sol[:],lw=3)


using DiffEqTutorials
DiffEqTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

