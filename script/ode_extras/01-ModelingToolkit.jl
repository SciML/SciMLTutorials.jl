
using ModelingToolkit

### Define a differential equation system

@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]
de = ODESystem(eqs)
ode_f = ODEFunction(de, [x,y,z], [σ,ρ,β])

### Use in DifferentialEquations.jl

using OrdinaryDiffEq
u₀ = ones(3)
tspan = (0.0,100.0)
p = [10.0,28.0,10/3]
prob = ODEProblem(ode_f,u₀,tspan,p)
sol = solve(prob,Tsit5())

using Plots
plot(sol,vars=(1,2,3))


generate_function(de, [x,y,z], [σ,ρ,β])


generate_function(de, [x,y,z], [σ,ρ,β]; version=ModelingToolkit.SArrayFunction)


jac = calculate_jacobian(de)


jac_expr = generate_jacobian(de)


ModelingToolkit.generate_factorized_W(de)[1]


@variables x y z
@parameters σ ρ β

# Define a nonlinear system
eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs, [x,y,z])
nlsys_func = generate_function(ns, [x,y,z], [σ,ρ,β])


nl_f = @eval eval(nlsys_func)
# Make a closure over the parameters for for NLsolve.jl
f2 = (du,u) -> nl_f(du,u,(10.0,26.0,2.33))

using NLsolve
nlsolve(f2,ones(3))


@derivatives D3'''~t
@derivatives D2''~t
@variables u(t), x(t)
eqs = [D3(u) ~ 2(D2(u)) + D(u) + D(x) + 1
       D2(x) ~ D(x) + 2]
de = ODESystem(eqs)
de1 = ode_order_lowering(de)


de1.eqs


@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t Dx'~x Dy'~y Dz'~z
eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]
J = [Dx(eqs[1].rhs) Dy(eqs[1].rhs) Dz(eqs[1].rhs)
 Dx(eqs[2].rhs) Dy(eqs[2].rhs) Dz(eqs[2].rhs)
 Dx(eqs[3].rhs) Dy(eqs[3].rhs) Dz(eqs[3].rhs)]


J = expand_derivatives.(J)


using LinearAlgebra
luJ = lu(J)


luJ.L


invJ = inv(J)


function lorenz(du,u,p,t)
 du[1] = p[1]*(u[2]-u[1])
 du[2] = u[1]*(p[2]-u[3]) - u[2]
 du[3] = u[1]*u[2] - p[3]*u[3]
end


u = [x,y,z]
du = similar(u)
p = [σ,ρ,β]
lorenz(du,u,p,t)
du


J = [Dx(du[1]) Dy(du[1]) Dz(du[1])
     Dx(du[2]) Dy(du[2]) Dz(du[2])
     Dx(du[3]) Dy(du[3]) Dz(du[3])]
J = expand_derivatives.(J)


using SparseArrays
function SparseArrays.SparseMatrixCSC(M::Matrix{T}) where {T<:ModelingToolkit.Expression}
    idxs = findall(!iszero, M)
    I = [i[1] for i in idxs]
    J = [i[2] for i in idxs]
    V = [M[i] for i in idxs]
    return SparseArrays.sparse_IJ_sorted!(I, J, V, size(M)...)
end
sJ = SparseMatrixCSC(J)


@parameters σ(..)
eqs = [D(x) ~ σ(t-1)*(y-x),
       D(y) ~ x*(σ(t^2)-z)-y,
       D(z) ~ x*y - β*z]


@derivatives Dₜ'~t
Dₜ(x*(σ(t^2)-z)-y)
expand_derivatives(Dₜ(x*(σ(t^2)-z)-y))


_f(x) = 2x + x^2
_f(x)


f(x) = 2x + x^2
@register f(x)


f(x)


function ModelingToolkit.derivative(::typeof(f), args::NTuple{1,Any}, ::Val{1})
    2 + 2args[1]
end
expand_derivatives(Dx(f(x)))


using PuMaS

subject = process_nmtran(example_nmtran_data("event_data/data23"),
                         [], [:ev1,:cp,:periph,:resp])[1]


mrm = @model begin
    @param   θ ∈ VectorDomain(12)
    @random  η ~ MvNormal(Matrix{Float64}(I, 11, 11))

    @pre begin
        Ka1     = θ[1]
        CL      = θ[2]*exp(η[1])
        Vc      = θ[3]*exp(η[2])
        Q       = θ[4]*exp(η[3])
        Vp      = θ[5]*exp(η[4])
        Kin     = θ[6]*exp(η[5])
        Kout    = θ[7]*exp(η[6])
        IC50    = θ[8]*exp(η[7])
        IMAX    = θ[9]*exp(η[8])
        γ       = θ[10]*exp(η[9])
        Vmax    = θ[11]*exp(η[10])
        Km      = θ[12]*exp(η[11])
    end

    @init begin
        Resp = θ[6]/θ[7]
    end

    @dynamics begin
        Ev1'    = -Ka1*Ev1
        Cent'   =  Ka1*Ev1 - (CL+Vmax/(Km+(Cent/Vc))+Q)*(Cent/Vc)  + Q*(Periph/Vp)
        Periph' =  Q*(Cent/Vc)  - Q*(Periph/Vp)
        Resp'   =  Kin*(1-(IMAX*(Cent/Vc)^γ/(IC50^γ+(Cent/Vc)^γ)))  - Kout*Resp
    end

    @derived begin
        ev1    = Ev1
        cp     = Cent / θ[3]
        periph = Periph
        resp   = Resp
    end
end


@dynamics begin
  Ev1'    = -Ka1*Ev1
  Cent'   =  Ka1*Ev1 - (CL+Vmax/(Km+(Cent/Vc))+Q)*(Cent/Vc)  + Q*(Periph/Vp)
  Periph' =  Q*(Cent/Vc)  - Q*(Periph/Vp)
  Resp'   =  Kin*(1-(IMAX*(Cent/Vc)^γ/(IC50^γ+(Cent/Vc)^γ)))  - Kout*Resp
end


param = (θ = [
              1, # Ka1  Absorption rate constant 1 (1/time)
              1, # CL   Clearance (volume/time)
              20, # Vc   Central volume (volume)
              2, # Q    Inter-compartmental clearance (volume/time)
              10, # Vp   Peripheral volume of distribution (volume)
              10, # Kin  Response in rate constant (1/time)
              2, # Kout Response out rate constant (1/time)
              2, # IC50 Concentration for 50% of max inhibition (mass/volume)
              1, # IMAX Maximum inhibition
              1, # γ    Emax model sigmoidicity
              0, # Vmax Maximum reaction velocity (mass/time)
              2  # Km   Michaelis constant (mass/volume)
              ],)
randeffs = (η = zeros(11),)
sim = simobs(mrm, subject, param, randeffs)


plot(sim)


using DiffEqTutorials
DiffEqTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

