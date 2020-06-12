
using Flux, StochasticDiffEq
using  NeuralNetDiffEq
using Plots
using LinearAlgebra


function phi(xi)
    y = Float64[]
    K = 100
    for x in eachcol(xi)
        val = max(K - maximum(x) , 0.00)
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
m = fmap(cu , m)
prob = KolmogorovPDEProblem(μ , σ , phi , xspan , tspan, d)
opt = Flux.ADAM(0.01)


sol = solve(prob, NeuralNetDiffEq.NNKolmogorov(m, opt, sdealg, ensemblealg), verbose = true, dt = dt,
            dx = dx , trajectories = trajectories , abstol=1e-6, maxiters = 3000)


monte_carlo_sol = []
x_out = collect(98:1.00:110.00)
for x in x_out
  u₀= [x]
  g_val(du , u , p , t) = 0.4.*u
  f_val(du , u , p , t) = du .= 0.04.*u
  dt = 0.01
  tspan = (0.0,1.0)
  prob = SDEProblem(f,g,u₀,tspan)
  output_func(sol,i) = (sol[end],false)
  ensembleprob_val = EnsembleProblem(prob , output_func = output_func )
  sim_val = solve(ensembleprob_val, EM(), EnsembleThreads() , dt=0.01, trajectories=100000,adaptive=false)
  s = reduce(hcat , sim_val.u)
  global monte_carlo_sol = push!(monte_carlo_sol , mean(phi(s)))
end


m = fmap(cpu , m)
x_model = reshape(x_out, 1 , size(x_out)[1])
y_out = m(x_model)
y_out = reshape(y_out , 13 , 1)


plot(x_out , y_out , lw = 3 ,  xaxis="Initial Stock Price", yaxis="Payoff")
plot!(x_out , monte_carlo_sol , lw = 3 ,  xaxis="Initial Stock Price", yaxis="Payoff")

