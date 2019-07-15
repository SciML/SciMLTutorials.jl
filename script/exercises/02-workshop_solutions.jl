
using DifferentialEquations, Plots
function orego(du,u,p,t)
  s,q,w = p
  y1,y2,y3 = u
  du[1] = s*(y2+y1*(1-q*y1-y2))
  du[2] = (y3-(1+y1)*y2)/s
  du[3] = w*(y1-y3)
end
p = [77.27,8.375e-6,0.161]
prob = ODEProblem(orego,[1.0,2.0,3.0],(0.0,360.0),p)
sol = solve(prob)
plot(sol)


plot(sol,vars=(1,2,3))


using BenchmarkTools
prob = ODEProblem(orego,[1.0,2.0,3.0],(0.0,50.0),p)
@btime sol = solve(prob,Tsit5())


@btime sol = solve(prob,Rodas5())


function orego(du,u,p,t)
  s,q,w = p
  y1,y2,y3 = u
  du[1] = s*(y2+y1*(1-q*y1-y2))
  du[2] = (y3-(1+y1)*y2)/s
  du[3] = w*(y1-y3)
end
function g(du,u,p,t)
  du[1] = 0.1u[1]
  du[2] = 0.1u[2]
  du[3] = 0.1u[3]
end
p = [77.27,8.375e-6,0.161]
prob = SDEProblem(orego,g,[1.0,2.0,3.0],(0.0,30.0),p)
sol = solve(prob,SOSRI())
plot(sol)


sol = solve(prob,ImplicitRKMil()); plot(sol)


sol = solve(prob,ImplicitRKMil()); plot(sol)


function orego(du,u,p,t)
  s,q,w = p
  y1,y2,y3 = u
  du[1] = s*(y2+y1*(1-q*y1-y2))
  du[2] = (y3-(1+y1)*y2)/s
  du[3] = w*(y1-y3)
end
p = [60.0,1e-5,0.2]
prob = ODEProblem(orego,[1.0,2.0,3.0],(0.0,30.0),p)
sol = solve(prob,Rodas5(),abstol=1/10^14,reltol=1/10^14)


function onecompartment(du,u,p,t)
  Ka,Ke = p
  du[1] = -Ka*u[1]
  du[2] =  Ka*u[1] - Ke*u[2]
end
p = (Ka=2.268,Ke=0.07398)
prob = ODEProblem(onecompartment,[100.0,0.0],(0.0,90.0),p)

tstops = [24,48,72]
condition(u,t,integrator) = t ∈ tstops
affect!(integrator) = (integrator.u[1] += 100)
cb = DiscreteCallback(condition,affect!)
sol = solve(prob,Tsit5(),callback=cb,tstops=tstops)
plot(sol)


function onecompartment_delay(du,u,h,p,t)
  Ka,Ke,τ = p
  delayed_depot = h(p,t-τ)[1]
  du[1] = -Ka*u[1]
  du[2] =  Ka*delayed_depot - Ke*u[2]
end
p = (Ka=2.268,Ke=0.07398,τ=6.0)
h(p,t) = [0.0,0.0]
prob = DDEProblem(onecompartment_delay,[100.0,0.0],h,(0.0,90.0),p)

tstops = [24,48,72]
condition(u,t,integrator) = t ∈ tstops
affect!(integrator) = (integrator.u[1] += 100)
cb = DiscreteCallback(condition,affect!)
sol = solve(prob,MethodOfSteps(Rosenbrock23()),callback=cb,tstops=tstops)
plot(sol)


p = (Ka = 0.5, Ke = 0.1, τ = 4.0)

