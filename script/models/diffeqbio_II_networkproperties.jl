
using DiffEqBiological, Latexify
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


x = latexify(rn; env=:chemical, starred=true, mathjax=true);
display("text/latex", "$x");


species(rn)


params(rn)


substrates(rn, 11)


substrates.(rn, 1:numreactions(rn))


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


addreaction!(rnmin, :(hillr(D₂,α,K,n)), (), (:m₁ => 1,))
addreaction!(rnmin, :k₊, (:P₁=>1, :P₂=>1), (:T=>1,))
addreaction!(rnmin, :k₋, (:T=>1,), (:P₁=>1, :P₂=>1))


setdiff(species(rn), species(rnmin))


setdiff(params(rn), params(rnmin))


rxidx = numreactions(rn)
setdiff(substrates(rn, rxidx), substrates(rnmin, rxidx))


setdiff(products(rn, rxidx), products(rnmin, rxidx))


rateexpr(rn, rxidx) == rateexpr(rnmin, rxidx)

