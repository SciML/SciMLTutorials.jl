
using DiffEqBayes, OrdinaryDiffEq, RecursiveArrayTools, Distributions, Plots, StatsPlots, BenchmarkTools, TransformVariables, CmdStan, DynamicHMC


function pendulum(du,u,p,t)
    ω,L = p
    x,y = u
    du[1] = y
    du[2] = - ω*y -(9.8/L)*sin(x)
end

u0 = [1.0,0.1]
tspan = (0.0,10.0)
prob1 = ODEProblem(pendulum,u0,tspan,[1.0,2.5])


sol = solve(prob1,Tsit5())
plot(sol)


t = collect(range(1,stop=10,length=10))
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)


scatter!(data')


priors = [Uniform(0.1,3.0), Normal(3.0,1.0)]


bayesian_result = turing_inference(prob1,Tsit5(),t,data,priors;num_samples=10_000,
                                   syms = [:omega,:L])


plot(bayesian_result)


plot(bayesian_result, colordim = :parameter)


@btime bayesian_result = turing_inference(prob1,Tsit5(),t,data,priors;syms = [:omega,:L],num_samples=10_000)


@btime bayesian_result = stan_inference(prob1,t,data,priors;num_samples=10_000,printsummary=false)


@btime bayesian_result = dynamichmc_inference(prob1,Tsit5(),t,data,priors;num_samples = 10_000)

