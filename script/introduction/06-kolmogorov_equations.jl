
using Flux, StochasticDiffEq
using  NeuralNetDiffEq
using Distributions
using LinearAlgebra


function phi(xi)
    y = Float64[]
    K = 100
    for x in eachcol(xi)
        val = max(maximum(x) - K , 0.00)
        y = push!(y , val)
    end
    y = reshape(y , 1 , size(y)[1] )
    return y
end


d = 1
r = 0.06
sigma = 0.4
xspan = (95.00 , 110.0)
tspan = (0.0 , 1.0)
σ(du , u , p , t) = du .= sigma.*u
μ(du , u , p , t) = du .= r.*u


sdealg = EM()
ensemblealg = EnsembleThreads()
dt = 0.01
dx = 0.01
trajectories = 100000


m = Chain(Dense(d, 32, leakyrelu),Dense(32, 128, leakyrelu),Dense(128 , 32 , leakyrelu) , Dense(32 , 1))
prob = KolmogorovPDEProblem(μ , σ , phi , xspan , tspan, d)
opt = Flux.ADAM(0.01)


sol = solve(prob, NeuralNetDiffEq.NNKolmogorov(m, opt, sdealg, ensemblealg), verbose = true, dt = dt,
            dx = dx , trajectories = trajectories , abstol=1e-6, maxiters = 3000)
