
f(u,p,t) = p.*u


using DifferentialEquations, Plots
u0 = [10.0]
p = [-0.3]
tspan = (0.0,10.0)
prob = ODEProblem(f,u0,tspan,p)
sol = solve(prob)
plot(sol)
ylims!(0.0,10.0)


using Distributions
u0_dist = [Uniform(-10.0,10.0)]


prob_func(prob,i,repeat) = remake(prob, u0 = rand.(u0_dist))
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)

ensemblesol = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=100000)


plot(ensemblesol, vars = (0,1), lw=1,alpha=0.1, label=nothing, idxs = 1:250)


g(sol) = sol(4.0)
mean([g(sol) for sol in ensemblesol])


using DiffEqUncertainty
expectation(g, prob, u0_dist, p, MonteCarlo(), Tsit5(); trajectories=100000)


expectation(g, prob, u0_dist, p, Koopman(), Tsit5())


exp(p[1]*4.0)*mean(u0_dist[1])


@time expectation(g, prob, u0_dist, p, MonteCarlo(), Tsit5(); trajectories=100000)


@time expectation(g, prob, u0_dist, p, Koopman(), Tsit5())


u0_dist = [Uniform(0.0,10.0)]
@time expectation(g, prob, u0_dist, p, MonteCarlo(), Tsit5(); trajectories=100_000)


@time expectation(g, prob, u0_dist, p, Koopman(), Tsit5())[1]


exp(p[1]*4.0)*mean(u0_dist[1])


u0_dist = [Normal(3.0,2.0)]
expectation(g, prob, u0_dist, p, Koopman(), Tsit5())


exp(p[1]*4.0)*mean(u0_dist[1])


u0_dist = [truncated(Normal(3.0,2.0),-1000,1000)]
expectation(g, prob, u0_dist, p, Koopman(), Tsit5())


u0_dist = [truncated(Normal(3.0,2.0),-5,11)]
expectation(g, prob, u0_dist, p, Koopman(), Tsit5())


u0_dist = [truncated(Normal(3.0,2.0),3-1000,3+1000)]
expectation(g, prob, u0_dist, p, Koopman(), Tsit5())


g(sol) = [sol(4.0)[1], sol(6.0)[1]]
expectation(g, prob, u0_dist, p, Koopman(), Tsit5(); nout = 2)


exp.(p.*[4.0,6.0])*mean(u0_dist[1])


saveat = tspan[1]:.5:tspan[2]
g(sol) = Matrix(sol)
mean_koop = expectation(g, prob, u0_dist, p, Koopman(), Tsit5(); nout = length(saveat), saveat=saveat)


plot(t->exp(p[1]*t)*mean(u0_dist[1]),tspan..., xlabel="t", label="analytical")
scatter!(collect(saveat),mean_koop.u[:],marker=:o, label=nothing)


function g(sol, power, counter)
    counter[power] = counter[power] + 1
    sol(4.0)[1]^power
end

counters = [0,0,0]
x_koop = expectation(s->g(s,1,counters), prob, u0_dist, p, Koopman(), Tsit5())
x2_koop = expectation(s->g(s,2,counters), prob, u0_dist, p, Koopman(), Tsit5())
x3_koop = expectation(s->g(s,3,counters), prob, u0_dist, p, Koopman(), Tsit5())
counters


function g(sol, counter) 
    counter[1] = counter[1] + 1
    v = sol(4.0)[1]
    [v, v^2, v^3]
end

counter = [0]
expectation(s->g(s,counter), prob, u0_dist, p, Koopman(), Tsit5(); nout = 3)
counter


function g(sol) 
    x = sol(4.0)[1]
    [x, x^2]
end

koop = expectation(g, prob, u0_dist, p, Koopman(), Tsit5(); nout = 2)
mean_koop = koop[1]
var_koop = koop[2] - mean_koop^2


exp(2*p[1]*4.0)*var(u0_dist[1])


saveat = tspan[1]:.5:tspan[2]
g(sol) = [Matrix(sol)'; (Matrix(sol).^2)']

koop = expectation(g, prob, u0_dist, p, Koopman(), Tsit5(); nout = length(saveat)*2, saveat=saveat)
μ = koop.u[1:length(saveat)]
σ = sqrt.(koop.u[length(saveat)+1:end] - μ.^2)

plot(t->exp(p[1]*t)*mean(u0_dist[1]),tspan..., ribbon = t->-sqrt(exp(2*p[1]*t)*var(u0_dist[1])), label="Analytical Mean, 1 std bounds")
scatter!(collect(saveat),μ,marker=:x, yerror = σ, c=:black, label = "Koopman Mean, 1 std bounds")


function g(sol) 
    v = sol(4.0)[1]
    [v, v^2, v^3]
end

koop = expectation(g, prob, u0_dist, p, Koopman(), Tsit5(); nout = 3)
mean_koop = koop[1]
var_koop = koop[2] - mean_koop^2
(koop[3] - 3.0*mean_koop*var_koop - mean_koop^3) / var_koop^(3/2)


g(sol) = sol(4.0)[1]
centralmoment(5, g, prob, u0_dist, p, Koopman(), Tsit5(),
                ireltol = 1e-9, iabstol = 1e-9)


using Quadrature, Cuba


u0_dist = [truncated(Normal(3.0,2.0),-5,11)]
p_dist = [truncated(Normal(-.7, .1), -1,0)]

g(sol) = sol(6.0)[1]

expectation(g, prob, u0_dist, p_dist, Koopman(), Tsit5(), EnsembleThreads(); 
                quadalg = CubaSUAVE(), batch=10)[1]


using BenchmarkTools

@btime expectation(g, prob, u0_dist, p_dist, Koopman(), Tsit5(); 
                quadalg = CubaSUAVE())[1]


@btime expectation(g, prob, u0_dist, p_dist, Koopman(), Tsit5(), EnsembleThreads(); 
                quadalg = CubaSUAVE(), batch=10)[1]


using DiffEqGPU

function f(du, u,p,t) 
    @inbounds begin
        du[1] = p[1]*u[1];
    end
    nothing
end

u0 = Float32[10.0]
p = Float32[-0.3]
tspan = (0.0f0,10.0f0)
prob = ODEProblem(f,u0,tspan,p)

g(sol) = sol(6.0)[1]

u0_dist = [truncated(Normal(3.0f0,2.0f0),-5f0,11f0)]
p_dist = [truncated(Normal(-.7f0, .1f0), -1f0,0f0)]

@btime expectation(g, prob, u0_dist, p_dist, Koopman(), Tsit5(), EnsembleGPUArray(); 
                   quadalg = CubaSUAVE(), batch=1000)[1]


using SciMLTutorials
SciMLTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])


using DiffEqGPU

function f(du, u,p,t) 
    @inbounds begin
        du[1] = p[1]*u[1];
    end
    nothing
end

u0 = Float32[10.0]
p = Float32[-0.3]
tspan = (0.0f0,10.0f0)
prob = ODEProblem(f,u0,tspan,p)

g(sol) = [sol(4.0)[1], sol(6.0)[1]]

u0_dist = [truncated(Normal(3.0f0,2.0f0),-5f0,11f0)]
p_dist = [truncated(Normal(-.7f0, .1f0), -1f0,0f0)]

expectation(g, prob, u0_dist, p_dist, Koopman(), Tsit5(), EnsembleGPUArray(); 
            quadalg = CubaSUAVE(), batch=1000, nout=2)[1]

