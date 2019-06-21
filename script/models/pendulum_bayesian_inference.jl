
using DiffEqBayes, OrdinaryDiffEq, ParameterizedFunctions, RecursiveArrayTools, Distributions, Plots


f1 = @ode_def begin
        dx = y
        dy = -(g/L)*sin(x)
    end g L

u0 = [1.0,0.1]
tspan = (0.0,10.0)
prob1 = ODEProblem(f1,u0,tspan,[9.79, 1.0])


sol = solve(prob1,Tsit5())
plot(sol)


t = collect(range(1,stop=10,length=10))
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)


scatter!(data')


priors = [Uniform(9.5,10.0), Uniform(0.9,1.1)]


bayesian_result = turing_inference(prob1,Tsit5(),t,data,priors;num_samples=15000)

