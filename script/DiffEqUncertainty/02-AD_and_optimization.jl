
using OrdinaryDiffEq, Plots

function ball!(du,u,p,t) 
    du[1] = u[2]
    du[2] = 0.0
    du[3] = u[4]
    du[4] = -p[1]
end

ground_condition(u,t,integrator) = u[3]
ground_affect!(integrator) = integrator.u[4] = -integrator.p[2] * integrator.u[4]
ground_cb = ContinuousCallback(ground_condition, ground_affect!)

u0 = [0.0,2.0,50.0,0.0]
tspan = (0.0,50.0)
p = [9.807, 0.9]

prob = ODEProblem(ball!,u0,tspan,p)
sol = solve(prob,Tsit5(),callback=ground_cb)
plot(sol, vars=(1,3), label = nothing, xlabel="x", ylabel="y")


stop_condition(u,t,integrator) = u[1] - 25.0
stop_cb = ContinuousCallback(stop_condition, terminate!)
cbs = CallbackSet(ground_cb, stop_cb)

tspan = (0.0, 1500.0)
prob = ODEProblem(ball!,u0,tspan,p)
sol = solve(prob,Tsit5(),callback=cbs)


rectangle(xc, yc, w, h) = Shape(xc .+ [-w,w,w,-w]./2.0, yc .+ [-h,-h,h,h]./2.0)

begin
    plot(sol, vars=(1,3), label=nothing, lw = 3, c=:black)
    xlabel!("x [m]")
    ylabel!("y [m]")
    plot!(rectangle(27.5, 25, 5, 50), c=:red, label = nothing)
    scatter!([25],[25],marker=:star, ms=10, label = nothing,c=:green)
    ylims!(0.0,50.0)
end


using Distributions

cor_dist = truncated(Normal(0.9, 0.02), 0.9-3*0.02, 1.0)
trajectories = 100000

prob_func(prob,i,repeat) = remake(prob, p = [p[1], rand(cor_dist)])
ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)
ensemblesol = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=trajectories, callback=cbs)

begin # plot
    plot(ensemblesol, vars = (1,3), lw=1,alpha=0.2, label=nothing, idxs = 1:350)
    xlabel!("x [m]")
    ylabel!("y [m]")
    plot!(rectangle(27.5, 25, 5, 50), c=:red, label = nothing)
    scatter!([25],[25],marker=:star, ms=10, label = nothing, c=:green)
    plot!(sol, vars=(1,3), label=nothing, lw = 3, c=:black, ls=:dash)
    xlims!(0.0,27.5)
end


obs(sol) = abs2(sol[3,end]-25)


mean_ensemble = mean([obs(sol) for sol in ensemblesol])


using DiffEqUncertainty

p_uncertain = [9.807, cor_dist]
expectation(obs, prob, u0, p_uncertain, Koopman(), Tsit5();
            ireltol = 1e-5, callback=cbs)


using NLopt, DiffEqSensitivity, ForwardDiff

make_u0(θ) = [θ[1],θ[2],θ[3], 0.0]

function 𝔼_loss(θ)   # \bbE
    u0 = make_u0(θ)
    expectation(obs, prob, u0, p_uncertain, Koopman(), Tsit5();
                 ireltol = 1e-5, callback=cbs)[1]
end


function 𝔼_loss_nlopt(x,∇)
    length(∇) > 0 ? ForwardDiff.gradient!(∇, 𝔼_loss,x) : nothing
    𝔼_loss(x)
end


opt = Opt(:LD_MMA, 3)
opt.lower_bounds = [-100.0,1.0, 10.0]
opt.upper_bounds = [0.0,3.0, 50.0]
opt.xtol_rel = 1e-3
opt.min_objective = 𝔼_loss_nlopt
(minf,minx,ret) = NLopt.optimize(opt, [-1.0, 2.0, 50.0])


ensembleprob = EnsembleProblem(remake(prob,u0 = make_u0(minx)),prob_func=prob_func)
ensemblesol = solve(ensembleprob,Tsit5(),EnsembleThreads(), trajectories=100_000, callback=cbs)

begin
    plot(ensemblesol, vars = (1,3), lw=1,alpha=0.1, label=nothing, idxs = 1:350)
    plot!(solve(remake(prob, u0=make_u0(minx)),Tsit5(), callback=cbs), 
            vars=(1,3),label = nothing, c=:black, lw=3,ls=:dash)
    xlabel!("x [m]")
    ylabel!("y [m]")
    plot!(rectangle(27.5, 25, 5, 50), c=:red, label = nothing)
    scatter!([25],[25],marker=:star, ms=10, label = nothing,c=:green)
    ylims!(0.0,50.0)
    xlims!(minx[1], 27.5)
end


using BenchmarkTools

@btime NLopt.optimize($opt, $[-1.0, 2.0, 50.0])


constraint = [20.0, 25.0]
begin
    plot(rectangle(27.5, 25, 5, 50), c=:red, label = nothing)
    xlabel!("x [m]")
    ylabel!("y [m]")
    plot!([constraint[1], constraint[1]],[0.0,constraint[2]], lw=5, c=:black, label=nothing)
    scatter!([25],[25],marker=:star, ms=10, label = nothing,c=:green)
    ylims!(0.0,50.0)
    xlims!(minx[1], 27.5)
end


constraint_condition(u,t,integrator) = u[1] - constraint[1]
constraint_affect!(integrator) = integrator.u[3] < constraint[2] ? terminate!(integrator) : nothing
constraint_cb = ContinuousCallback(constraint_condition, constraint_affect!, save_positions=(true,false));
constraint_cbs = CallbackSet(ground_cb, stop_cb, constraint_cb)

ensemblesol = solve(ensembleprob,Tsit5(),EnsembleThreads(), trajectories=350, callback=constraint_cbs, maxstep=0.1)

begin
    plot(ensemblesol, vars = (1,3), lw=1,alpha=0.1, label=nothing)
    plot!(solve(remake(prob, u0=make_u0(minx)),Tsit5(), callback=constraint_cbs), 
            vars=(1,3),label = nothing, c=:black, lw=3, ls=:dash)

    xlabel!("x [m]")
    ylabel!("y [m]")
    plot!(rectangle(27.5, 25, 5, 50), c=:red, label = nothing)
    plot!([constraint[1], constraint[1]],[0.0,constraint[2]], lw=5, c=:black)
    scatter!([25],[25],marker=:star, ms=10, label = nothing,c=:green)
    ylims!(0.0,50.0)
    xlims!(minx[1], 27.5)
end


constraint_obs(sol) = sol[1,end] ≈ constraint[1] ? one(sol[1,end]) : zero(sol[1,end])


expectation(constraint_obs, prob, make_u0(minx), p_uncertain, Koopman(), Tsit5();
            ireltol= 1e-9, iabstol = 1e-9, callback=constraint_cbs)[1]


function 𝔼_constraint(θ)
    u0 = [θ[1],θ[2],θ[3], 0.0]
    expectation(constraint_obs, prob, u0, p_uncertain, Koopman(), Tsit5(),
                ireltol= 1e-9, iabstol = 1e-9,callback=constraint_cbs)[1]
end

function 𝔼_constraint_nlopt(x,∇)
    length(∇) > 0 ? ForwardDiff.gradient!(∇, 𝔼_constraint,x) : nothing
    𝔼_constraint(x) - 0.01
end


opt = Opt(:LD_MMA, 3)
opt.lower_bounds = [-100.0, 1.0, 10.0]
opt.upper_bounds = [0.0, 3.0, 50.0]
opt.xtol_rel = 1e-3
opt.min_objective = 𝔼_loss_nlopt
inequality_constraint!(opt,𝔼_constraint_nlopt, 1e-5)
(minf2,minx2,ret2) = NLopt.optimize(opt, [-1.0, 2.0, 50.0])


λ = 𝔼_constraint(minx2)


λ - 0.01 <= 1e-5


ensembleprob = EnsembleProblem(remake(prob,u0 = make_u0(minx2)),prob_func=prob_func)
ensemblesol = solve(ensembleprob,Tsit5(),EnsembleThreads(), 
                    trajectories=350, callback=constraint_cbs)

begin
    plot(ensemblesol, vars = (1,3), lw=1,alpha=0.1, label=nothing)
    plot!(solve(remake(prob, u0=make_u0(minx2)),Tsit5(), callback=constraint_cbs), 
            vars=(1,3),label = nothing, c=:black, lw=3, ls=:dash)
    plot!([constraint[1], constraint[1]],[0.0,constraint[2]], lw=5, c=:black)

    xlabel!("x [m]")
    ylabel!("y [m]")
    plot!(rectangle(27.5, 25, 5, 50), c=:red, label = nothing)
    scatter!([25],[25],marker=:star, ms=10, label = nothing,c=:green)
    ylims!(0.0,50.0)
    xlims!(minx[1], 27.5)
end


using SciMLTutorials
SciMLTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

