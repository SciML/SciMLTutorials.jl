
const v0 = -84.624
const v1 = 10.0
const C_K1 = 1.0f0
const C_x1 = 1.0f0
const C_Na = 1.0f0
const C_s = 1.0f0
const D_Ca = 0.0f0
const D_Na = 0.0f0
const g_s = 0.09f0
const g_Na = 4.0f0
const g_NaC = 0.005f0
const ENa = 50.0f0 + D_Na
const γ = 0.5f0
const C_m = 1.0f0


mutable struct BeelerReuterCpu <: Function
    t::Float64              # the last timestep time to calculate Δt
    diff_coef::Float64      # the diffusion-coefficient (coupling strength)

    C::Array{Float32, 2}    # intracellular calcium concentration
    M::Array{Float32, 2}    # sodium current activation gate (m)
    H::Array{Float32, 2}    # sodium current inactivation gate (h)
    J::Array{Float32, 2}    # sodium current slow inactivaiton gate (j)
    D::Array{Float32, 2}    # calcium current activaiton gate (d)
    F::Array{Float32, 2}    # calcium current inactivation gate (f)
    XI::Array{Float32, 2}   # inward-rectifying potassium current (iK1)

    Δu::Array{Float64, 2}   # place-holder for the Laplacian

    function BeelerReuterCpu(u0, diff_coef)
        self = new()

        ny, nx = size(u0)
        self.t = 0.0
        self.diff_coef = diff_coef

        self.C = fill(0.0001f0, (ny,nx))
        self.M = fill(0.01f0, (ny,nx))
        self.H = fill(0.988f0, (ny,nx))
        self.J = fill(0.975f0, (ny,nx))
        self.D = fill(0.003f0, (ny,nx))
        self.F = fill(0.994f0, (ny,nx))
        self.XI = fill(0.0001f0, (ny,nx))

        self.Δu = zeros(ny,nx)

        return self
    end
end


# 5-point stencil
function laplacian(Δu, u)
    n1, n2 = size(u)

    # internal nodes
    for j = 2:n2-1
        for i = 2:n1-1
            @inbounds  Δu[i,j] = u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1] - 4*u[i,j]
        end
    end

    # left/right edges
    for i = 2:n1-1
        @inbounds Δu[i,1] = u[i+1,1] + u[i-1,1] + 2*u[i,2] - 4*u[i,1]
        @inbounds Δu[i,n2] = u[i+1,n2] + u[i-1,n2] + 2*u[i,n2-1] - 4*u[i,n2]
    end

    # top/bottom edges
    for j = 2:n2-1
        @inbounds Δu[1,j] = u[1,j+1] + u[1,j-1] + 2*u[2,j] - 4*u[1,j]
        @inbounds Δu[n1,j] = u[n1,j+1] + u[n1,j-1] + 2*u[n1-1,j] - 4*u[n1,j]
    end

    # corners
    @inbounds Δu[1,1] = 2*(u[2,1] + u[1,2]) - 4*u[1,1]
    @inbounds Δu[n1,1] = 2*(u[n1-1,1] + u[n1,2]) - 4*u[n1,1]
    @inbounds Δu[1,n2] = 2*(u[2,n2] + u[1,n2-1]) - 4*u[1,n2]
    @inbounds Δu[n1,n2] = 2*(u[n1-1,n2] + u[n1,n2-1]) - 4*u[n1,n2]
end


@inline function rush_larsen(g, α, β, Δt)
    inf = α/(α+β)
    τ = 1f0 / (α+β)
    return clamp(g + (g - inf) * expm1(-Δt/τ), 0f0, 1f0)
end


function update_M_cpu(g, v, Δt)
    # the condition is needed here to prevent NaN when v == 47.0
    α = isapprox(v, 47.0f0) ? 10.0f0 : -(v+47.0f0) / (exp(-0.1f0*(v+47.0f0)) - 1.0f0)
    β = (40.0f0 * exp(-0.056f0*(v+72.0f0)))
    return rush_larsen(g, α, β, Δt)
end

function update_H_cpu(g, v, Δt)
    α = 0.126f0 * exp(-0.25f0*(v+77.0f0))
    β = 1.7f0 / (exp(-0.082f0*(v+22.5f0)) + 1.0f0)
   return rush_larsen(g, α, β, Δt)
end

function update_J_cpu(g, v, Δt)
    α = (0.55f0 * exp(-0.25f0*(v+78.0f0))) / (exp(-0.2f0*(v+78.0f0)) + 1.0f0)
    β = 0.3f0 / (exp(-0.1f0*(v+32.0f0)) + 1.0f0)
    return rush_larsen(g, α, β, Δt)
end

function update_D_cpu(g, v, Δt)
    α = γ * (0.095f0 * exp(-0.01f0*(v-5.0f0))) / (exp(-0.072f0*(v-5.0f0)) + 1.0f0)
    β = γ * (0.07f0 * exp(-0.017f0*(v+44.0f0))) / (exp(0.05f0*(v+44.0f0)) + 1.0f0)
    return rush_larsen(g, α, β, Δt)
end

function update_F_cpu(g, v, Δt)
    α = γ * (0.012f0 * exp(-0.008f0*(v+28.0f0))) / (exp(0.15f0*(v+28.0f0)) + 1.0f0)
    β = γ * (0.0065f0 * exp(-0.02f0*(v+30.0f0))) / (exp(-0.2f0*(v+30.0f0)) + 1.0f0)
    return rush_larsen(g, α, β, Δt)
end

function update_XI_cpu(g, v, Δt)
    α = (0.0005f0 * exp(0.083f0*(v+50.0f0))) / (exp(0.057f0*(v+50.0f0)) + 1.0f0)
    β = (0.0013f0 * exp(-0.06f0*(v+20.0f0))) / (exp(-0.04f0*(v+20.0f0)) + 1.0f0)
    return rush_larsen(g, α, β, Δt)
end


function update_C_cpu(g, d, f, v, Δt)
    ECa = D_Ca - 82.3f0 - 13.0278f0 * log(g)
    kCa = C_s * g_s * d * f
    iCa = kCa * (v - ECa)
    inf = 1.0f-7 * (0.07f0 - g)
    τ = 1f0 / 0.07f0
    return g + (g - inf) * expm1(-Δt/τ)
end


function update_gates_cpu(u, XI, M, H, J, D, F, C, Δt)
    let Δt = Float32(Δt)
        n1, n2 = size(u)
        for j = 1:n2
            for i = 1:n1
                v = Float32(u[i,j])

                XI[i,j] = update_XI_cpu(XI[i,j], v, Δt)
                M[i,j] = update_M_cpu(M[i,j], v, Δt)
                H[i,j] = update_H_cpu(H[i,j], v, Δt)
                J[i,j] = update_J_cpu(J[i,j], v, Δt)
                D[i,j] = update_D_cpu(D[i,j], v, Δt)
                F[i,j] = update_F_cpu(F[i,j], v, Δt)

                C[i,j] = update_C_cpu(C[i,j], D[i,j], F[i,j], v, Δt)
            end
        end
    end
end


# iK1 is the inward-rectifying potassium current
function calc_iK1(v)
    ea = exp(0.04f0*(v+85f0))
    eb = exp(0.08f0*(v+53f0))
    ec = exp(0.04f0*(v+53f0))
    ed = exp(-0.04f0*(v+23f0))
    return 0.35f0 * (4f0*(ea-1f0)/(eb + ec)
            + 0.2f0 * (isapprox(v, -23f0) ? 25f0 : (v+23f0) / (1f0-ed)))
end

# ix1 is the time-independent background potassium current
function calc_ix1(v, xi)
    ea = exp(0.04f0*(v+77f0))
    eb = exp(0.04f0*(v+35f0))
    return xi * 0.8f0 * (ea-1f0) / eb
end

# iNa is the sodium current (similar to the classic Hodgkin-Huxley model)
function calc_iNa(v, m, h, j)
    return C_Na * (g_Na * m^3 * h * j + g_NaC) * (v - ENa)
end

# iCa is the calcium current
function calc_iCa(v, d, f, c)
    ECa = D_Ca - 82.3f0 - 13.0278f0 * log(c)    # ECa is the calcium reversal potential
    return C_s * g_s * d * f * (v - ECa)
end

function update_du_cpu(du, u, XI, M, H, J, D, F, C)
    n1, n2 = size(u)

    for j = 1:n2
        for i = 1:n1
            v = Float32(u[i,j])

            # calculating individual currents
            iK1 = calc_iK1(v)
            ix1 = calc_ix1(v, XI[i,j])
            iNa = calc_iNa(v, M[i,j], H[i,j], J[i,j])
            iCa = calc_iCa(v, D[i,j], F[i,j], C[i,j])

            # total current
            I_sum = iK1 + ix1 + iNa + iCa

            # the reaction part of the reaction-diffusion equation
            du[i,j] = -I_sum / C_m
        end
    end
end


function (f::BeelerReuterCpu)(du, u, p, t)
    Δt = t - f.t

    if Δt != 0 || t == 0
        update_gates_cpu(u, f.XI, f.M, f.H, f.J, f.D, f.F, f.C, Δt)
        f.t = t
    end

    laplacian(f.Δu, u)

    # calculate the reaction portion
    update_du_cpu(du, u, f.XI, f.M, f.H, f.J, f.D, f.F, f.C)

    # ...add the diffusion portion
    du .+= f.diff_coef .* f.Δu
end


const N = 192;
u0 = fill(v0, (N, N));
u0[90:102,90:102] .= v1;   # a small square in the middle of the domain


using Plots
heatmap(u0)


using DifferentialEquations, Sundials

deriv_cpu = BeelerReuterCpu(u0, 1.0);
prob = ODEProblem(deriv_cpu, u0, (0.0, 50.0));


@time sol = solve(prob, CVODE_BDF(linear_solver=:GMRES), saveat=100.0);


heatmap(sol.u[end])


using CUDAnative, CuArrays

mutable struct BeelerReuterGpu <: Function
    t::Float64                  # the last timestep time to calculate Δt
    diff_coef::Float64          # the diffusion-coefficient (coupling strength)

    d_C::CuArray{Float32, 2}    # intracellular calcium concentration
    d_M::CuArray{Float32, 2}    # sodium current activation gate (m)
    d_H::CuArray{Float32, 2}    # sodium current inactivation gate (h)
    d_J::CuArray{Float32, 2}    # sodium current slow inactivaiton gate (j)
    d_D::CuArray{Float32, 2}    # calcium current activaiton gate (d)
    d_F::CuArray{Float32, 2}    # calcium current inactivation gate (f)
    d_XI::CuArray{Float32, 2}   # inward-rectifying potassium current (iK1)

    d_u::CuArray{Float64, 2}    # place-holder for u in the device memory
    d_du::CuArray{Float64, 2}   # place-holder for d_u in the device memory

    Δv::Array{Float64, 2}       # place-holder for voltage gradient

    function BeelerReuterGpu(u0, diff_coef)
        self = new()

        ny, nx = size(u0)
        @assert (nx % 16 == 0) && (ny % 16 == 0)
        self.t = 0.0
        self.diff_coef = diff_coef

        self.d_C = CuArray(fill(0.0001f0, (ny,nx)))
        self.d_M = CuArray(fill(0.01f0, (ny,nx)))
        self.d_H = CuArray(fill(0.988f0, (ny,nx)))
        self.d_J = CuArray(fill(0.975f0, (ny,nx)))
        self.d_D = CuArray(fill(0.003f0, (ny,nx)))
        self.d_F = CuArray(fill(0.994f0, (ny,nx)))
        self.d_XI = CuArray(fill(0.0001f0, (ny,nx)))

        self.d_u = CuArray(u0)
        self.d_du = CuArray(zeros(ny,nx))

        self.Δv = zeros(ny,nx)

        return self
    end
end


function rush_larsen_gpu(g, α, β, Δt)
    inf = α/(α+β)
    τ = 1.0/(α+β)
    return clamp(g + (g - inf) * CUDAnative.expm1(-Δt/τ), 0f0, 1f0)
end

function update_M_gpu(g, v, Δt)
    # the condition is needed here to prevent NaN when v == 47.0
    α = isapprox(v, 47.0f0) ? 10.0f0 : -(v+47.0f0) / (CUDAnative.exp(-0.1f0*(v+47.0f0)) - 1.0f0)
    β = (40.0f0 * CUDAnative.exp(-0.056f0*(v+72.0f0)))
    return rush_larsen_gpu(g, α, β, Δt)
end

function update_H_gpu(g, v, Δt)
    α = 0.126f0 * CUDAnative.exp(-0.25f0*(v+77.0f0))
    β = 1.7f0 / (CUDAnative.exp(-0.082f0*(v+22.5f0)) + 1.0f0)
    return rush_larsen_gpu(g, α, β, Δt)
end

function update_J_gpu(g, v, Δt)
    α = (0.55f0 * CUDAnative.exp(-0.25f0*(v+78.0f0))) / (CUDAnative.exp(-0.2f0*(v+78.0f0)) + 1.0f0)
    β = 0.3f0 / (CUDAnative.exp(-0.1f0*(v+32.0f0)) + 1.0f0)
    return rush_larsen_gpu(g, α, β, Δt)
end

function update_D_gpu(g, v, Δt)
    α = γ * (0.095f0 * CUDAnative.exp(-0.01f0*(v-5.0f0))) / (CUDAnative.exp(-0.072f0*(v-5.0f0)) + 1.0f0)
    β = γ * (0.07f0 * CUDAnative.exp(-0.017f0*(v+44.0f0))) / (CUDAnative.exp(0.05f0*(v+44.0f0)) + 1.0f0)
    return rush_larsen_gpu(g, α, β, Δt)
end

function update_F_gpu(g, v, Δt)
    α = γ * (0.012f0 * CUDAnative.exp(-0.008f0*(v+28.0f0))) / (CUDAnative.exp(0.15f0*(v+28.0f0)) + 1.0f0)
    β = γ * (0.0065f0 * CUDAnative.exp(-0.02f0*(v+30.0f0))) / (CUDAnative.exp(-0.2f0*(v+30.0f0)) + 1.0f0)
    return rush_larsen_gpu(g, α, β, Δt)
end

function update_XI_gpu(g, v, Δt)
    α = (0.0005f0 * CUDAnative.exp(0.083f0*(v+50.0f0))) / (CUDAnative.exp(0.057f0*(v+50.0f0)) + 1.0f0)
    β = (0.0013f0 * CUDAnative.exp(-0.06f0*(v+20.0f0))) / (CUDAnative.exp(-0.04f0*(v+20.0f0)) + 1.0f0)
    return rush_larsen_gpu(g, α, β, Δt)
end

function update_C_gpu(c, d, f, v, Δt)
    ECa = D_Ca - 82.3f0 - 13.0278f0 * CUDAnative.log(c)
    kCa = C_s * g_s * d * f
    iCa = kCa * (v - ECa)
    inf = 1.0f-7 * (0.07f0 - c)
    τ = 1f0 / 0.07f0
    return c + (c - inf) * CUDAnative.expm1(-Δt/τ)
end


# iK1 is the inward-rectifying potassium current
function calc_iK1(v)
    ea = CUDAnative.exp(0.04f0*(v+85f0))
    eb = CUDAnative.exp(0.08f0*(v+53f0))
    ec = CUDAnative.exp(0.04f0*(v+53f0))
    ed = CUDAnative.exp(-0.04f0*(v+23f0))
    return 0.35f0 * (4f0*(ea-1f0)/(eb + ec)
            + 0.2f0 * (isapprox(v, -23f0) ? 25f0 : (v+23f0) / (1f0-ed)))
end

# ix1 is the time-independent background potassium current
function calc_ix1(v, xi)
    ea = CUDAnative.exp(0.04f0*(v+77f0))
    eb = CUDAnative.exp(0.04f0*(v+35f0))
    return xi * 0.8f0 * (ea-1f0) / eb
end

# iNa is the sodium current (similar to the classic Hodgkin-Huxley model)
function calc_iNa(v, m, h, j)
    return C_Na * (g_Na * m^3 * h * j + g_NaC) * (v - ENa)
end

# iCa is the calcium current
function calc_iCa(v, d, f, c)
    ECa = D_Ca - 82.3f0 - 13.0278f0 * CUDAnative.log(c)    # ECa is the calcium reversal potential
    return C_s * g_s * d * f * (v - ECa)
end


function update_gates_gpu(u, XI, M, H, J, D, F, C, Δt)
    i = (blockIdx().x-UInt32(1)) * blockDim().x + threadIdx().x
    j = (blockIdx().y-UInt32(1)) * blockDim().y + threadIdx().y

    v = Float32(u[i,j])

    let Δt = Float32(Δt)
        XI[i,j] = update_XI_gpu(XI[i,j], v, Δt)
        M[i,j] = update_M_gpu(M[i,j], v, Δt)
        H[i,j] = update_H_gpu(H[i,j], v, Δt)
        J[i,j] = update_J_gpu(J[i,j], v, Δt)
        D[i,j] = update_D_gpu(D[i,j], v, Δt)
        F[i,j] = update_F_gpu(F[i,j], v, Δt)

        C[i,j] = update_C_gpu(C[i,j], D[i,j], F[i,j], v, Δt)
    end
    nothing
end

function update_du_gpu(du, u, XI, M, H, J, D, F, C)
    i = (blockIdx().x-UInt32(1)) * blockDim().x + threadIdx().x
    j = (blockIdx().y-UInt32(1)) * blockDim().y + threadIdx().y

    v = Float32(u[i,j])

    # calculating individual currents
    iK1 = calc_iK1(v)
    ix1 = calc_ix1(v, XI[i,j])
    iNa = calc_iNa(v, M[i,j], H[i,j], J[i,j])
    iCa = calc_iCa(v, D[i,j], F[i,j], C[i,j])

    # total current
    I_sum = iK1 + ix1 + iNa + iCa

    # the reaction part of the reaction-diffusion equation
    du[i,j] = -I_sum / C_m
    nothing
end


function (f::BeelerReuterGpu)(du, u, p, t)
    L = 16   # block size
    Δt = t - f.t
    copyto!(f.d_u, u)
    ny, nx = size(u)

    if Δt != 0 || t == 0
        @cuda blocks=(ny÷L,nx÷L) threads=(L,L) update_gates_gpu(
            f.d_u, f.d_XI, f.d_M, f.d_H, f.d_J, f.d_D, f.d_F, f.d_C, Δt)
        f.t = t
    end

    laplacian(f.Δv, u)

    # calculate the reaction portion
    @cuda blocks=(ny÷L,nx÷L) threads=(L,L) update_du_gpu(
        f.d_du, f.d_u, f.d_XI, f.d_M, f.d_H, f.d_J, f.d_D, f.d_F, f.d_C)

    copyto!(du, f.d_du)

    # ...add the diffusion portion
    du .+= f.diff_coef .* f.Δv
end


using DifferentialEquations, Sundials

deriv_gpu = BeelerReuterGpu(u0, 1.0);
prob = ODEProblem(deriv_gpu, u0, (0.0, 50.0));
@time sol = solve(prob, CVODE_BDF(linear_solver=:GMRES), saveat=100.0);


heatmap(sol.u[end])


using DiffEqTutorials
DiffEqTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

