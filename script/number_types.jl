
using DifferentialEquations
f = (u,p,t) -> (p*u)
prob_ode_linear = ODEProblem(f,1/2,(0.0,1.0),1.01);


prob = prob_ode_linear
sol =solve(prob,Tsit5())
println(sol)


prob = ODEProblem(f,1/2,(0//1,1//1),101//100);
sol = solve(prob,RK4(),dt=1//2^(6),adaptive=false)
println(sol)


prob = ODEProblem(f,BigInt(1)//BigInt(2),(0//1,1//1),101//100);
sol =solve(prob,RK4(),dt=1//2^(6),adaptive=false)
println(sol[end])


prob_ode_biglinear = ODEProblem(f,big(1.0)/big(2.0),(big(0.0),big(1.0)),big(1.01))
sol =solve(prob_ode_biglinear,Tsit5())
println(sol[end])


using DoubleFloats
prob_ode_doublelinear = ODEProblem(f,Double64(1)/Double64(2),(Double64(0),Double64(1)),Double64(1.01))
sol =solve(prob_ode_doublelinear,Tsit5())
println(sol[end])


using ArbNumerics
prob_ode_arbfloatlinear = ODEProblem(f,ArbFloat(1)/ArbFloat(2),(ArbFloat(0.0),ArbFloat(1.0)),ArbFloat(1.01))
sol =solve(prob_ode_arbfloatlinear,Tsit5())
println(sol[end])


using DecFP
prob_ode_decfplinear = ODEProblem(f,Dec128(1)/Dec128(2),(Dec128(0.0),Dec128(1.0)),Dec128(1.01))
sol =solve(prob_ode_decfplinear,Tsit5())
println(sol[end]); println(typeof(sol[end]))


using Decimals
prob_ode_decimallinear = ODEProblem(f,[decimal("1.0")]./[decimal("2.0")],(0//1,1//1),decimal(1.01))
sol =solve(prob_ode_decimallinear,RK4(),dt=1/2^(6)) #Fails
println(sol[end]); println(typeof(sol[end]))


using DiffEqTutorials
DiffEqTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

