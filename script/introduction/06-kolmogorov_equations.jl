
using Flux, StochasticDiffEq
using  NeuralNetDiffEq
using Distributions
using Plots
using LinearAlgebra


g(u , p , t) = 0.5
f(u , p , t) = -u
d = 1
sdealg = EM()
ensemblealg = EnsembleThreads()
u0 = Normal(0.00 , 0.1)
xspan = (-3.0 , 3.0)
tspan = (0.0 , 1.0)
prob = SDEProblem(f , g , u0 , tspan ; xspan = xspan , d = d)
opt = Flux.ADAM(0.05)
m = Chain(Dense(1, 16, elu),Dense(16, 64, elu), Dense(64 , 5 , elu) , Dense(5 , 1))
m = fmap(f64 , m)
sol = solve(prob, NNKolmogorov(m,opt , sdealg,EnsembleDistributed()) , verbose = true, dt = 0.1,
            abstol=1e-10, dx = 0.001 , trajectories = 50000 ,  maxiters = 500 )


T = tspan[2]
mu = 1
sig = 0.5
sig2 = (sig^2 /(2*mu))exp(-2*mu*T)*(exp(2*mu*T) - 1)
analytical(xi) = pdf.(Normal(0.00 , sqrt(sig2)) , xi)


xs = xspan[1]:0.01:xspan[2]
x_val = collect(xs)
x_val= reshape(x_val , 1 , size(x_val)[1])
m = testmode!(m)
y_val = m(x_val)
y_val = reshape(y_val , length(xs) , 1)
x_val = collect(xs)
plot(x_val , y_val,linewidth=3,title="Solution to the linear PDE with a thick line",
 xaxis="Time (t)",yaxis="u(t) (in μm)",label="NNKolmogorov")
plot!(x_val , analytical(x_val),linewidth=3,title="Solution to the linear PDE with a thick line",
      xaxis="Time (t)",yaxis="u(t) (in μm)",label="Analytical")

