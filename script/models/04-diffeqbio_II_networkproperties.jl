
using DifferentialEquations, DiffEqBiological, Latexify, Plots
fmt = :svg
pyplot(fmt=fmt)
rn = @reaction_network begin
    hillr(D₂,α,K,n), ∅ --> m₁
    hillr(D₁,α,K,n), ∅ --> m₂
    (δ,γ), m₁ ↔ ∅
    (δ,γ), m₂ ↔ ∅
    β, m₁ --> m₁ + P₁
    β, m₂ --> m₂ + P₂
    μ, P₁ --> ∅
    μ, P₂ --> ∅
    (k₊,k₋), 2P₁ ↔ D₁ 
    (k₊,k₋), 2P₂ ↔ D₂
    (k₊,k₋), P₁+P₂ ↔ T
end α K n δ γ β μ k₊ k₋;


latexify(rn; env=:chemical)


x = latexify(rn; env=:chemical, starred=true, mathjax=false);
display("text/latex", "$x");


species(rn)


params(rn)


substratesymstoich(rn, 11)


substratesymstoich.(rn, 1:numreactions(rn))


netstoich.(rn, 1:numreactions(rn))


rxtospecies_depgraph(rn)


species(rn)[[3,4,7]]


speciestorx_depgraph(rn)[1]


findall(depset -> in(:m₁, depset), dependents.(rn, 1:numreactions(rn)))


rxtorx_depgraph(rn)


rnmin = @min_reaction_network begin
    (δ,γ), m₁ ↔ ∅
    (δ,γ), m₂ ↔ ∅
    β, m₁ --> m₁ + P₁
    β, m₂ --> m₂ + P₂
    μ, P₁ --> ∅
    μ, P₂ --> ∅
end δ γ β μ;


addspecies!(rnmin, :D₁)
addspecies!(rnmin, :D₂)
addspecies!(rnmin, :T)


addparam!(rnmin, :α)
addparam!(rnmin, :K)
addparam!(rnmin, :n)
addparam!(rnmin, :k₊)
addparam!(rnmin, :k₋)


addreaction!(rnmin, :(hillr(D₁,α,K,n)), :(∅ --> m₂))
addreaction!(rnmin, :((k₊,k₋)), :(2P₂ ↔ D₂))
addreaction!(rnmin, :k₊, :(2P₁ --> D₁))
addreaction!(rnmin, :k₋, :(D₁ --> 2P₁))


# signature is addreaction!(rnmin, paramexpr, substratestoich, productstoich)
addreaction!(rnmin, :(hillr(D₂,α,K,n)), (), (:m₁ => 1,))
addreaction!(rnmin, :k₊, (:P₁=>1, :P₂=>1), (:T=>1,))
addreaction!(rnmin, :k₋, (:T=>1,), (:P₁=>1, :P₂=>1))


setdiff(species(rn), species(rnmin))


setdiff(params(rn), params(rnmin))


rxidx = numreactions(rn)
setdiff(substrates(rn, rxidx), substrates(rnmin, rxidx))


setdiff(products(rn, rxidx), products(rnmin, rxidx))


rateexpr(rn, rxidx) == rateexpr(rnmin, rxidx)


addodes!(rnmin)


odeexprs(rnmin)


latexify(rnmin)


x = latexify(rnmin, starred=true);
display("text/latex", "$x");


latexify(jacobianexprs(rnmin))


x = latexify(jacobianexprs(rnmin), starred=true);
display("text/latex", "$x");


N = 64
h = 1 / N


rn = @empty_reaction_network

for i = 1:N
    addspecies!(rn, Symbol(:u, i))
end


addparam!(rn, :β)


for i = 1:N
    (i < N) && addreaction!(rn, :β, (Symbol(:u,i)=>1,), (Symbol(:u,i+1)=>1,))
    (i > 1) && addreaction!(rn, :β, (Symbol(:u,i)=>1,), (Symbol(:u,i-1)=>1,))
end


addodes!(rn)


u₀ = zeros(N)
u₀[div(N,2)] = 10000
p = [1/(h*h)]
tspan = (0.,.01)
oprob = ODEProblem(rn, u₀, tspan, p)


sol = solve(oprob, Rodas5())
times = [0., .0001, .001, .01]
plt = plot()
for time in times
    plot!(plt, 1:N, sol(time), fmt=fmt, xlabel="i", ylabel="uᵢ", label=string("t = ", time), lw=3)
end
plot(plt, ylims=(0.,10000.))


addjumps!(rn, build_regular_jumps=false, minimal_jumps=true)

# make the initial condition integer valued 
u₀ = zeros(Int, N)
u₀[div(N,2)] = 10000

# setup and solve the problem
dprob = DiscreteProblem(rn, u₀, tspan, p)
jprob = JumpProblem(dprob, DirectCR(), rn, save_positions=(false,false))
jsol = solve(jprob, SSAStepper(), saveat=times)


times = [0., .0001, .001, .01]
plts = []
for i = 1:4
    b = bar(1:N, jsol[i], legend=false, fmt=fmt, xlabel="i", ylabel="uᵢ", title=string("t = ", times[i]))
    plot!(b,sol(times[i]))
    push!(plts,b)
end
plot(plts...)


using DiffEqTutorials
DiffEqTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file], remove_homedir=true)

