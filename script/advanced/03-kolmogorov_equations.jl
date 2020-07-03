
using Flux, StochasticDiffEq
using NeuralNetDiffEq
using Plots
using CuArrays
using CUDAnative


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
r = 0.04
sigma = 0.2
xspan = (80.00 , 115.0)
tspan = (0.0 , 1.0)
σ(du , u , p , t) = du .= sigma.*u
μ(du , u , p , t) = du .= r.*u
prob = KolmogorovPDEProblem(μ , σ , phi , xspan , tspan, d)


sdealg = EM()
ensemblealg = EnsembleThreads()
dt = 0.01
dx = 0.01
trajectories = 100000


m = Chain(Dense(d, 64, elu),Dense(64, 128, elu),Dense(128 , 16 , elu) , Dense(16 , 1))
use_gpu = false
if CUDAnative.functional() == true
  m = fmap(CuArrays.cu , m)
  use_gpu = true
end
opt = Flux.ADAM(0.0005)


@time sol = solve(prob, NeuralNetDiffEq.NNKolmogorov(m, opt, sdealg, ensemblealg), verbose = true, dt = dt,
            dx = dx , trajectories = trajectories , abstol=1e-6, maxiters = 1000 , use_gpu = use_gpu)


monte_carlo_sol = []
x_out = collect(85:2.00:110.00)
for x in x_out
  u₀= [x]
  g_val(du , u , p , t) = du .= 0.2.*u
  f_val(du , u , p , t) = du .= 0.04.*u
  dt = 0.01
  tspan = (0.0,1.0)
  prob = SDEProblem(f_val,g_val,u₀,tspan)
  output_func(sol,i) = (sol[end],false)
  ensembleprob_val = EnsembleProblem(prob , output_func = output_func )
  sim_val = solve(ensembleprob_val, EM(), EnsembleThreads() , dt=0.01, trajectories=100000,adaptive=false)
  s = reduce(hcat , sim_val.u)
  mean_phi = sum(phi(s))/length(phi(s))
  global monte_carlo_sol = push!(monte_carlo_sol , mean_phi)
end


x_model = reshape(x_out, 1 , size(x_out)[1])
if use_gpu == true
  m = fmap(cpu , m)
end
y_out = m(x_model)
y_out = reshape(y_out , 13 , 1)


plot(x_out , y_out , lw = 3 ,  xaxis="Initial Stock Price", yaxis="Payoff" , label = "NNKolmogorov")
plot!(x_out , monte_carlo_sol , lw = 3 ,  xaxis="Initial Stock Price", yaxis="Payoff" ,label = "Monte Carlo Solutions")

