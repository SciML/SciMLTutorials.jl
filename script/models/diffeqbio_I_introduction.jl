
# If not already installed, first hit "]" within a Julia REPL. Then type:
# add DifferentialEquations DiffEqBiological PyPlot Plots Latexify 

using DifferentialEquations, DiffEqBiological, Plots, Latexify
pyplot(fmt=:svg);


repressilator = @reaction_network begin
    hillr(P₃,α,K,n), ∅ --> m₁
    hillr(P₁,α,K,n), ∅ --> m₂
    hillr(P₂,α,K,n), ∅ --> m₃
    (δ,γ), m₁ ↔ ∅
    (δ,γ), m₂ ↔ ∅
    (δ,γ), m₃ ↔ ∅
    β, m₁ --> m₁ + P₁
    β, m₂ --> m₂ + P₂
    β, m₃ --> m₃ + P₃
    μ, P₁ --> ∅
    μ, P₂ --> ∅
    μ, P₃ --> ∅
end α K n δ γ β μ;


latexify(repressilator; env=:chemical)


x = latexify(repressilator; env=:chemical, starred=true, mathjax=true);
display("text/latex", "$x");


latexify(repressilator)


x = latexify(repressilator, starred=true);
display("text/latex", "$x");


speciesmap(repressilator)


paramsmap(repressilator)


# parameters [α,K,n,δ,γ,β,μ]
p = (.5, 40, 2, log(2)/120, 5e-3, 20*log(2)/120, log(2)/60)

# initial condition [m₁,m₂,m₃,P₁,P₂,P₃]
u₀ = [0.,0.,0.,20.,0.,0.]

# time interval to solve on
tspan = (0., 10000.)

# create the ODEProblem we want to solve
oprob = ODEProblem(repressilator, u₀, tspan, p)


sol = solve(oprob, saveat=10.)
plot(sol, fmt=:svg)


# first we redefine the initial condition to be integer valued
u₀ = [0,0,0,20,0,0]

# next we create a discrete problem to encode that our species are integer valued:
dprob = DiscreteProblem(repressilator, u₀, tspan, p)

# now we create a JumpProblem, and specify Gillespie's Direct Method as the solver:
jprob = JumpProblem(dprob, Direct(), repressilator, save_positions=(false,false))

# now let's solve and plot the jump process:
sol = solve(jprob, SSAStepper(), saveat=10.)
plot(sol, fmt=:svg)


rjs = regularjumps(repressilator)
lprob = JumpProblem(dprob, Direct(), rjs)
lsol = solve(lprob, SimpleTauLeaping(), dt=.1)
plot(lsol, plotdensity=1000, fmt=:svg)


bdp = @reaction_network begin
  c₁, X --> 2X
  c₂, X --> 0
  c₃, 0 --> X
end c₁ c₂ c₃
p = (1.0,2.0,50.)
u₀ = [5.]
tspan = (0.,4.);


# SDEProblem for CLE
sprob = SDEProblem(bdp, u₀, tspan, p)

# solve and plot, tstops is used to specify enough points 
# that the plot looks well-resolved
sol = solve(sprob, tstops=range(0., step=4e-3, length=1001))
plot(sol, fmt=:svg)


latexify(jacobianexprs(repressilator))


x = latexify(jacobianexprs(repressilator), starred=true);
display("text/latex", "$x");


using DiffEqTutorials
DiffEqTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file], remove_homedir=true)

