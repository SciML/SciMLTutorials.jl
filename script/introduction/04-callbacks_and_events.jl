
using DifferentialEquations, ParameterizedFunctions
ball! = @ode_def BallBounce begin
  dy =  v
  dv = -g
end g


function condition(u,t,integrator)
  u[1]
end


function affect!(integrator)
    integrator.u[2] = -integrator.p[2] * integrator.u[2]
end


bounce_cb = ContinuousCallback(condition,affect!)


u0 = [50.0,0.0]
tspan = (0.0,15.0)
p = (9.8,0.9)
prob = ODEProblem(ball!,u0,tspan,p,callback=bounce_cb)


sol = solve(prob,Tsit5())
using Plots; gr()
plot(sol)


function condition_kick(u,t,integrator)
    t == 2
end


function affect_kick!(integrator)
    integrator.u[2] += 50
end


kick_cb = DiscreteCallback(condition_kick,affect_kick!)
u0 = [50.0,0.0]
tspan = (0.0,10.0)
p = (9.8,0.9)
prob = ODEProblem(ball!,u0,tspan,p,callback=kick_cb)


sol = solve(prob,Tsit5(),tstops=[2.0])
plot(sol)


cb = CallbackSet(bounce_cb,kick_cb)


u0 = [50.0,0.0]
tspan = (0.0,15.0)
p = (9.8,0.9)
prob = ODEProblem(ball!,u0,tspan,p,callback=cb)
sol = solve(prob,Tsit5(),tstops=[2.0])
plot(sol)


u0 = [1.,0.]
harmonic! = @ode_def HarmonicOscillator begin
   dv = -x
   dx = v
end
tspan = (0.0,10.0)
prob = ODEProblem(harmonic!,u0,tspan)
sol = solve(prob)
plot(sol)


function terminate_affect!(integrator)
    terminate!(integrator)
end


function terminate_condition(u,t,integrator)
    u[2]
end
terminate_cb = ContinuousCallback(terminate_condition,terminate_affect!)


sol = solve(prob,callback=terminate_cb)
plot(sol)


sol.t[end]


terminate_upcrossing_cb = ContinuousCallback(terminate_condition,terminate_affect!,nothing)


sol = solve(prob,callback=terminate_upcrossing_cb)
plot(sol)


tspan = (0.0,10000.0)
prob = ODEProblem(harmonic!,u0,tspan)
sol = solve(prob)
gr(fmt=:png) # Make it a PNG instead of an SVG since there's a lot of points!
plot(sol,vars=(1,2))


plot(sol,vars=(0,1),denseplot=false)


plot(sol.t,[u[2]^2 + u[1]^2 for u in sol.u]) # Energy ~ x^2 + v^2


function g(resid,u,p,t)
  resid[1] = u[2]^2 + u[1]^2 - 1
  resid[2] = 0
end


cb = ManifoldProjection(g)
sol = solve(prob,callback=cb)
plot(sol,vars=(1,2))


plot(sol,vars=(0,1),denseplot=false)


u1,u2 = sol[500]
u2^2 + u1^2


prob = ODEProblem((du,u,p,t)->du.=u,rand(1000,1000),(0.0,1.0))


saved_values = SavedValues(Float64, Tuple{Float64,Float64})


using LinearAlgebra
cb = SavingCallback((u,t,integrator)->(tr(u),norm(u)), saved_values)


sol = solve(prob, Tsit5(), callback=cb, save_everystep=false, save_start=false, save_end = false) # Turn off normal saving


saved_values.t


saved_values.saveval


saved_values = SavedValues(Float64, Tuple{Float64,Float64}) # New cache
cb = SavingCallback((u,t,integrator)->(tr(u),norm(u)), saved_values, saveat = 0.0:0.1:1.0)
sol = solve(prob, Tsit5(), callback=cb, save_everystep=false, save_start=false, save_end = false) # Turn off normal saving


saved_values.t


saved_values.saveval


using DiffEqTutorials
DiffEqTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

