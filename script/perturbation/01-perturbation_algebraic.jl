
using Symbolics, SymbolicUtils

function solve_newton(f, x, x₀; abstol=1e-8, maxiter=50)
    xₙ = Float64(x₀)
    fₙ₊₁ = x - f / Symbolics.derivative(f, x)

    for i = 1:maxiter
        xₙ₊₁ = substitute(fₙ₊₁, Dict(x => xₙ))
        if abs(xₙ₊₁ - xₙ) < abstol
            return xₙ₊₁
        else
            xₙ = xₙ₊₁
        end
    end
    return xₙ₊₁
end


n = 2
@variables ϵ a[1:n]


x = 1 + a[1]*ϵ + a[2]*ϵ^2


  eq = x^5 + ϵ*x - 1


expand(eq)


function collect_powers(eq, x, ns; max_power=100)
    eq = substitute(expand(eq), Dict(x^j => 0 for j=last(ns)+1:max_power))

    eqs = []
    for i in ns
        powers = Dict(x^j => (i==j ? 1 : 0) for j=1:last(ns))
        push!(eqs, substitute(eq, powers))
    end
    eqs
end


eqs = collect_powers(eq, ϵ, 1:2)


substitute(expand(eq), Dict(
  ϵ => 0,
  ϵ^2 => 1,
  ϵ^3 => 0,
  ϵ^4 => 0,
  ϵ^5 => 0,
  ϵ^6 => 0,
  ϵ^7 => 0,
  ϵ^8 => 0)
)


function solve_coef(eqs, ps)
    vals = Dict()

    for i = 1:length(ps)
        eq = substitute(eqs[i], vals)
        vals[ps[i]] = Symbolics.solve_for(eq ~ 0, ps[i])
    end
    vals
end


solve_coef(eqs, a)


X = 𝜀 -> 1 + a[1]*𝜀 + a[2]*𝜀^2


def_taylor(x, ps) = sum([a*x^i for (i,a) in enumerate(ps)])
def_taylor(x, ps, p₀) = p₀ + def_taylor(x, ps)


n = 3
@variables ϵ M a[1:n]
x = def_taylor(ϵ, a, M)


expand_sin(x, n) = sum([(isodd(k) ? -1 : 1)*(-x)^(2k-1)/factorial(2k-1) for k=1:n])


expand_sin(0.1, 10) ≈ sin(0.1)


eq = x - ϵ * expand_sin(x, n) - M


eqs = collect_powers(eq, ϵ, 1:n)


vals = solve_coef(eqs, a)


x′ = substitute(x, vals)
X = (𝜀, 𝑀) -> substitute(x′, Dict(ϵ => 𝜀, M => 𝑀))
X(0.01671, π/2)

