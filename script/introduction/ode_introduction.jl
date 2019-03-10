
f(u,p,t) = 0.98u


using DifferentialEquations
f(u,p,t) = 0.98u
u0 = 1.0
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)


sol = solve(prob)


using Plots; gr()
plot(sol)


plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false


plot!(sol.t, t->1.0*exp(0.98t),lw=3,ls=:dash,label="True Solution!")


sol.t


sol.u


[t+u for (u,t) in tuples(sol)]


sol


sol(0.45)


sol = solve(prob,abstol=1e-8,reltol=1e-8)


plot(sol)
plot!(sol.t, t->1.0*exp(0.98t),lw=3,ls=:dash,label="True Solution!")


sol = solve(prob,saveat=0.1)


sol = solve(prob,saveat=[0.2,0.7,0.9])


sol = solve(prob,dense=false)


sol = solve(prob,save_everystep=false)


sol = solve(prob,save_everystep=false,save_start = false)


sol = solve(prob,alg_hints=[:stiff])


sol = solve(prob,Tsit5(),reltol=1e-6)


function lorenz!(du,u,p,t)
    σ,ρ,β = p
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end


u0 = [1.0,0.0,0.0]


p = (10,28,8/3) # we could also make this an array, or any other type!


tspan = (0.0,100.0)
prob = ODEProblem(lorenz!,u0,tspan,p)


sol = solve(prob)


sol.t[10],sol[10]


sol[2,10]


A = Array(sol)


plot(sol)


plot(sol,vars=(1,2,3))


plot(sol,vars=(1,2,3),denseplot=false)


plot(sol,vars=(0,2))


function lotka_volterra!(du,u,p,t)
  du[1] = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] = -p[3]*u[2] + p[4]*u[1]*u[2]
end


using ParameterizedFunctions
lv! = @ode_def LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d


u0 = [1.0,1.0]
p = (1.5,1.0,3.0,1.0)
tspan = (0.0,10.0)
prob = ODEProblem(lv!,u0,tspan,p)
sol = solve(prob)
plot(sol)


lv!.Jex


A  = [1. 0  0 -5
      4 -2  4 -3
     -4  0  0  1
      5 -2  2  3]
u0 = rand(4,2)
tspan = (0.0,1.0)
f(u,p,t) = A*u
prob = ODEProblem(f,u0,tspan)
sol = solve(prob)


sol[3]


big_u0 = big.(u0)


prob = ODEProblem(f,big_u0,tspan)
sol = solve(prob)


sol[1,3]


prob = ODEProblem(f,big_u0,big.(tspan))
sol = solve(prob)


using StaticArrays
A  = @SMatrix [ 1.0  0.0 0.0 -5.0
                4.0 -2.0 4.0 -3.0
               -4.0  0.0 0.0  1.0
                5.0 -2.0 2.0  3.0]
u0 = @SMatrix rand(4,2)
tspan = (0.0,1.0)
f(u,p,t) = A*u
prob = ODEProblem(f,u0,tspan)
sol = solve(prob)


sol[3]


using DiffEqTutorials
DiffEqTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

