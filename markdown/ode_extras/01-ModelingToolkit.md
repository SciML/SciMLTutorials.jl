---
author: "Chris Rackauckas"
title: "ModelingToolkit.jl, An IR and Compiler for Scientific Models"
---


A lot of people are building modeling languages for their specific domains. However, while the syntax my vary greatly between these domain-specific languages (DSLs), the internals of modeling frameworks are surprisingly similar: building differential equations, calculating Jacobians, etc.

#### ModelingToolkit.jl is metamodeling systemitized

After building our third modeling interface, we realized that this problem can be better approached by having a reusable internal structure which DSLs can target. This internal is ModelingToolkit.jl: an Intermediate Representation (IR) with a well-defined interface for defining system transformations and compiling to Julia functions for use in numerical libraries. Now a DSL can easily be written by simply defining the translation to ModelingToolkit.jl's primatives and querying for the mathematical quantities one needs.

### Basic usage: defining differential equation systems, with performance!

Let's explore the IR itself. ModelingToolkit.jl is friendly to use, and can used as a symbolic DSL in its own right. Let's define and solve the Lorenz differential equation system using ModelingToolkit to generate the functions:

````julia
using ModelingToolkit

### Define a differential equation system

@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]
de = ODESystem(eqs, t, [x,y,z], [σ,ρ,β])
ode_f = ODEFunction(de)

### Use in DifferentialEquations.jl

using OrdinaryDiffEq
u₀ = ones(3)
tspan = (0.0,100.0)
p = [10.0,28.0,10/3]
prob = ODEProblem(ode_f,u₀,tspan,p)
sol = solve(prob,Tsit5())

using Plots
plot(sol,vars=(1,2,3))
````


![](figures/01-ModelingToolkit_1_1.png)



### ModelingToolkit is a compiler for mathematical systems

At its core, ModelingToolkit is a compiler. It's IR is its type system, and its output are Julia functions (it's a compiler for Julia code to Julia code, written in Julia).

DifferentialEquations.jl wants a function `f(u,p,t)` or `f(du,u,p,t)` for defining an ODE system,
so ModelingToolkit.jl builds both. First the out of place version:

````julia
generate_function(de)[1]
````


````
:((var"##MTKArg#345", var"##MTKArg#346", var"##MTKArg#347")->begin
          if var"##MTKArg#345" isa Array || !(typeof(var"##MTKArg#345") <: StaticArray) && false
              return @inbounds(begin
                          let (x, y, z, σ, ρ, β, t) = (var"##MTKArg#345"[1], var"##MTKArg#345"[2], var"##MTKArg#345"[3], var"##MTK
Arg#346"[1], var"##MTKArg#346"[2], var"##MTKArg#346"[3], var"##MTKArg#347")
                              [σ * (y - x), x * (ρ - z) - y, x * y - β * z]
                          end
                      end)
          else
              X = @inbounds(begin
                          let (x, y, z, σ, ρ, β, t) = (var"##MTKArg#345"[1], var"##MTKArg#345"[2], var"##MTKArg#345"[3], var"##MTK
Arg#346"[1], var"##MTKArg#346"[2], var"##MTKArg#346"[3], var"##MTKArg#347")
                              (σ * (y - x), x * (ρ - z) - y, x * y - β * z)
                          end
                      end)
              construct = if var"##MTKArg#345" isa ModelingToolkit.StaticArrays.StaticArray
                      ModelingToolkit.StaticArrays.similar_type(typeof(var"##MTKArg#345"), eltype(X))
                  else
                      x->begin
                              convert(typeof(var"##MTKArg#345"), x)
                          end
                  end
              return construct(X)
          end
      end)
````





and the in-place:

````julia
generate_function(de)[2]
````


````
:((var"##MTIIPVar#354", var"##MTKArg#350", var"##MTKArg#351", var"##MTKArg#352")->begin
          @inbounds begin
                  let (x, y, z, σ, ρ, β, t) = (var"##MTKArg#350"[1], var"##MTKArg#350"[2], var"##MTKArg#350"[3], var"##MTKArg#351"
[1], var"##MTKArg#351"[2], var"##MTKArg#351"[3], var"##MTKArg#352")
                      var"##MTIIPVar#354"[1] = σ * (y - x)
                      var"##MTIIPVar#354"[2] = x * (ρ - z) - y
                      var"##MTIIPVar#354"[3] = x * y - β * z
                  end
              end
          nothing
      end)
````





ModelingToolkit.jl can be used to calculate the Jacobian of the differential equation system:

````julia
jac = calculate_jacobian(de)
````


````
3×3 Array{ModelingToolkit.Expression,2}:
   σ * -1             σ  Constant(0)
 ρ - z(t)  Constant(-1)    x(t) * -1
     y(t)          x(t)          -1β
````





It will automatically generate functions for using this Jacobian within the stiff ODE solvers for faster solving:

````julia
jac_expr = generate_jacobian(de)
````


````
(:((var"##MTKArg#355", var"##MTKArg#356", var"##MTKArg#357")->begin
          if var"##MTKArg#355" isa Array || !(typeof(var"##MTKArg#355") <: StaticArray) && false
              return @inbounds(begin
                          let (x, y, z, σ, ρ, β, t) = (var"##MTKArg#355"[1], var"##MTKArg#355"[2], var"##MTKArg#355"[3], var"##MTK
Arg#356"[1], var"##MTKArg#356"[2], var"##MTKArg#356"[3], var"##MTKArg#357")
                              [σ * -1 σ 0; ρ - z -1 x * -1; y x -1β]
                          end
                      end)
          else
              X = @inbounds(begin
                          let (x, y, z, σ, ρ, β, t) = (var"##MTKArg#355"[1], var"##MTKArg#355"[2], var"##MTKArg#355"[3], var"##MTK
Arg#356"[1], var"##MTKArg#356"[2], var"##MTKArg#356"[3], var"##MTKArg#357")
                              (σ * -1, ρ - z, y, σ, -1, x, 0, x * -1, -1β)
                          end
                      end)
              construct = if var"##MTKArg#355" isa ModelingToolkit.StaticArrays.StaticArray
                      ModelingToolkit.StaticArrays.SMatrix{3, 3}
                  else
                      x->begin
                              out = similar(typeof(var"##MTKArg#355"), 3, 3)
                              out .= x
                          end
                  end
              return construct(X)
          end
      end), :((var"##MTIIPVar#359", var"##MTKArg#355", var"##MTKArg#356", var"##MTKArg#357")->begin
          @inbounds begin
                  let (x, y, z, σ, ρ, β, t) = (var"##MTKArg#355"[1], var"##MTKArg#355"[2], var"##MTKArg#355"[3], var"##MTKArg#356"
[1], var"##MTKArg#356"[2], var"##MTKArg#356"[3], var"##MTKArg#357")
                      var"##MTIIPVar#359"[1] = σ * -1
                      var"##MTIIPVar#359"[2] = ρ - z
                      var"##MTIIPVar#359"[3] = y
                      var"##MTIIPVar#359"[4] = σ
                      var"##MTIIPVar#359"[5] = -1
                      var"##MTIIPVar#359"[6] = x
                      var"##MTIIPVar#359"[7] = 0
                      var"##MTIIPVar#359"[8] = x * -1
                      var"##MTIIPVar#359"[9] = -1β
                  end
              end
          nothing
      end))
````





It can even do fancy linear algebra. Stiff ODE solvers need to perform an LU-factorization which is their most expensive part. But ModelingToolkit.jl can skip this operation and instead generate the analytical solution to a matrix factorization, and build a Julia function for directly computing the factorization, which is then optimized in LLVM compiler passes.

````julia
ModelingToolkit.generate_factorized_W(de)[1]
````


````
(:((var"##MTKArg#360", var"##MTKArg#361", var"##MTKArg#362", var"##MTKArg#363")->begin
          if var"##MTKArg#360" isa Array || !(typeof(var"##MTKArg#360") <: StaticArray) && false
              return @inbounds(begin
                          let (x, y, z, σ, ρ, β, __MTKWgamma, t) = (var"##MTKArg#360"[1], var"##MTKArg#360"[2], var"##MTKArg#360"[
3], var"##MTKArg#361"[1], var"##MTKArg#361"[2], var"##MTKArg#361"[3], var"##MTKArg#362", var"##MTKArg#363")
                              [σ * -1 * __MTKWgamma + -1 __MTKWgamma * σ 0; __MTKWgamma * (ρ - z) * inv(σ * -1 * __MTKWgamma + -1)
 (__MTKWgamma * -1 + -1) - __MTKWgamma * (ρ - z) * inv(σ * -1 * __MTKWgamma + -1) * __MTKWgamma * σ x * -1 * __MTKWgamma - 0; __MT
KWgamma * y * inv(σ * -1 * __MTKWgamma + -1) (__MTKWgamma * x - __MTKWgamma * y * inv(σ * -1 * __MTKWgamma + -1) * __MTKWgamma * σ
) * inv((__MTKWgamma * -1 + -1) - __MTKWgamma * (ρ - z) * inv(σ * -1 * __MTKWgamma + -1) * __MTKWgamma * σ) ((-1 * β * __MTKWgamma
 + -1) - 0) - (__MTKWgamma * x - __MTKWgamma * y * inv(σ * -1 * __MTKWgamma + -1) * __MTKWgamma * σ) * inv((__MTKWgamma * -1 + -1)
 - __MTKWgamma * (ρ - z) * inv(σ * -1 * __MTKWgamma + -1) * __MTKWgamma * σ) * (x * -1 * __MTKWgamma - 0)]
                          end
                      end)
          else
              X = @inbounds(begin
                          let (x, y, z, σ, ρ, β, __MTKWgamma, t) = (var"##MTKArg#360"[1], var"##MTKArg#360"[2], var"##MTKArg#360"[
3], var"##MTKArg#361"[1], var"##MTKArg#361"[2], var"##MTKArg#361"[3], var"##MTKArg#362", var"##MTKArg#363")
                              (σ * -1 * __MTKWgamma + -1, __MTKWgamma * (ρ - z) * inv(σ * -1 * __MTKWgamma + -1), __MTKWgamma * y 
* inv(σ * -1 * __MTKWgamma + -1), __MTKWgamma * σ, (__MTKWgamma * -1 + -1) - __MTKWgamma * (ρ - z) * inv(σ * -1 * __MTKWgamma + -1
) * __MTKWgamma * σ, (__MTKWgamma * x - __MTKWgamma * y * inv(σ * -1 * __MTKWgamma + -1) * __MTKWgamma * σ) * inv((__MTKWgamma * -
1 + -1) - __MTKWgamma * (ρ - z) * inv(σ * -1 * __MTKWgamma + -1) * __MTKWgamma * σ), 0, x * -1 * __MTKWgamma - 0, ((-1 * β * __MTK
Wgamma + -1) - 0) - (__MTKWgamma * x - __MTKWgamma * y * inv(σ * -1 * __MTKWgamma + -1) * __MTKWgamma * σ) * inv((__MTKWgamma * -1
 + -1) - __MTKWgamma * (ρ - z) * inv(σ * -1 * __MTKWgamma + -1) * __MTKWgamma * σ) * (x * -1 * __MTKWgamma - 0))
                          end
                      end)
              construct = (x->begin
                          A = SMatrix{(3, 3)...}(x)
                          StaticArrays.LU(LowerTriangular(SMatrix{(3, 3)...}(UnitLowerTriangular(A))), UpperTriangular(A), SVector
(ntuple((n->begin
                                              n
                                          end), max((3, 3)...))))
                      end)
              return construct(X)
          end
      end), :((var"##MTIIPVar#365", var"##MTKArg#360", var"##MTKArg#361", var"##MTKArg#362", var"##MTKArg#363")->begin
          @inbounds begin
                  let (x, y, z, σ, ρ, β, __MTKWgamma, t) = (var"##MTKArg#360"[1], var"##MTKArg#360"[2], var"##MTKArg#360"[3], var"
##MTKArg#361"[1], var"##MTKArg#361"[2], var"##MTKArg#361"[3], var"##MTKArg#362", var"##MTKArg#363")
                      var"##MTIIPVar#365"[1] = σ * -1 * __MTKWgamma + -1
                      var"##MTIIPVar#365"[2] = __MTKWgamma * (ρ - z) * inv(σ * -1 * __MTKWgamma + -1)
                      var"##MTIIPVar#365"[3] = __MTKWgamma * y * inv(σ * -1 * __MTKWgamma + -1)
                      var"##MTIIPVar#365"[4] = __MTKWgamma * σ
                      var"##MTIIPVar#365"[5] = (__MTKWgamma * -1 + -1) - __MTKWgamma * (ρ - z) * inv(σ * -1 * __MTKWgamma + -1) * 
__MTKWgamma * σ
                      var"##MTIIPVar#365"[6] = (__MTKWgamma * x - __MTKWgamma * y * inv(σ * -1 * __MTKWgamma + -1) * __MTKWgamma *
 σ) * inv((__MTKWgamma * -1 + -1) - __MTKWgamma * (ρ - z) * inv(σ * -1 * __MTKWgamma + -1) * __MTKWgamma * σ)
                      var"##MTIIPVar#365"[7] = 0
                      var"##MTIIPVar#365"[8] = x * -1 * __MTKWgamma - 0
                      var"##MTIIPVar#365"[9] = ((-1 * β * __MTKWgamma + -1) - 0) - (__MTKWgamma * x - __MTKWgamma * y * inv(σ * -1
 * __MTKWgamma + -1) * __MTKWgamma * σ) * inv((__MTKWgamma * -1 + -1) - __MTKWgamma * (ρ - z) * inv(σ * -1 * __MTKWgamma + -1) * _
_MTKWgamma * σ) * (x * -1 * __MTKWgamma - 0)
                  end
              end
          nothing
      end))
````





### Solving Nonlinear systems

ModelingToolkit.jl is not just for differential equations. It can be used for any mathematical target that is representable by its IR. For example, let's solve a rootfinding problem `F(x)=0`. What we do is define a nonlinear system and generate a function for use in NLsolve.jl

````julia
@variables x y z
@parameters σ ρ β

# Define a nonlinear system
eqs = [0 ~ σ*(y-x),
       0 ~ x*(ρ-z)-y,
       0 ~ x*y - β*z]
ns = NonlinearSystem(eqs, [x,y,z], [σ,ρ,β])
nlsys_func = generate_function(ns)
````


````
(:((var"##MTKArg#373", var"##MTKArg#374")->begin
          if var"##MTKArg#373" isa Array || !(typeof(var"##MTKArg#373") <: StaticArray) && false
              return @inbounds(begin
                          let (x, y, z, σ, ρ, β) = (var"##MTKArg#373"[1], var"##MTKArg#373"[2], var"##MTKArg#373"[3], var"##MTKArg
#374"[1], var"##MTKArg#374"[2], var"##MTKArg#374"[3])
                              [(*)(σ, (-)(y, x)), (-)((*)(x, (-)(ρ, z)), y), (-)((*)(x, y), (*)(β, z))]
                          end
                      end)
          else
              X = @inbounds(begin
                          let (x, y, z, σ, ρ, β) = (var"##MTKArg#373"[1], var"##MTKArg#373"[2], var"##MTKArg#373"[3], var"##MTKArg
#374"[1], var"##MTKArg#374"[2], var"##MTKArg#374"[3])
                              ((*)(σ, (-)(y, x)), (-)((*)(x, (-)(ρ, z)), y), (-)((*)(x, y), (*)(β, z)))
                          end
                      end)
              construct = if var"##MTKArg#373" isa ModelingToolkit.StaticArrays.StaticArray
                      ModelingToolkit.StaticArrays.similar_type(typeof(var"##MTKArg#373"), eltype(X))
                  else
                      x->begin
                              convert(typeof(var"##MTKArg#373"), x)
                          end
                  end
              return construct(X)
          end
      end), :((var"##MTIIPVar#376", var"##MTKArg#373", var"##MTKArg#374")->begin
          @inbounds begin
                  let (x, y, z, σ, ρ, β) = (var"##MTKArg#373"[1], var"##MTKArg#373"[2], var"##MTKArg#373"[3], var"##MTKArg#374"[1]
, var"##MTKArg#374"[2], var"##MTKArg#374"[3])
                      var"##MTIIPVar#376"[1] = (*)(σ, (-)(y, x))
                      var"##MTIIPVar#376"[2] = (-)((*)(x, (-)(ρ, z)), y)
                      var"##MTIIPVar#376"[3] = (-)((*)(x, y), (*)(β, z))
                  end
              end
          nothing
      end))
````





We can then tell ModelingToolkit.jl to compile this function for use in NLsolve.jl, and then numerically solve the rootfinding problem:

````julia
nl_f = @eval eval(nlsys_func[2])
# Make a closure over the parameters for for NLsolve.jl
f2 = (du,u) -> nl_f(du,u,(10.0,26.0,2.33))

using NLsolve
nlsolve(f2,ones(3))
````


````
Results of Nonlinear Solver Algorithm
 * Algorithm: Trust-region with dogleg and autoscaling
 * Starting Point: [1.0, 1.0, 1.0]
 * Zero: [2.2228042243306243e-10, 2.2228042243645056e-10, -9.990339599422887e-11]
 * Inf-norm of residuals: 0.000000
 * Iterations: 3
 * Convergence: true
   * |x - x'| < 0.0e+00: false
   * |f(x)| < 1.0e-08: true
 * Function Calls (f): 4
 * Jacobian Calls (df/dx): 4
````





### Library of transformations on mathematical systems

The reason for using ModelingToolkit is not just for defining performant Julia functions for solving systems, but also for performing mathematical transformations which may be required in order to numerically solve the system. For example, let's solve a third order ODE. The way this is done is by transforming the third order ODE into a first order ODE, and then solving the resulting ODE. This transformation is given by the `ode_order_lowering` function.

````julia
@derivatives D3'''~t
@derivatives D2''~t
@variables u(t), x(t)
eqs = [D3(u) ~ 2(D2(u)) + D(u) + D(x) + 1
       D2(x) ~ D(x) + 2]
de = ODESystem(eqs, t, [u,x], [])
de1 = ode_order_lowering(de)
````


````
ModelingToolkit.ODESystem(ModelingToolkit.Equation[ModelingToolkit.Equation(derivative(uˍtt(t), t), ((2 * uˍtt(t) + uˍt(t)) + xˍt(
t)) + 1), ModelingToolkit.Equation(derivative(xˍt(t), t), xˍt(t) + 2), ModelingToolkit.Equation(derivative(uˍt(t), t), uˍtt(t)), M
odelingToolkit.Equation(derivative(u(t), t), uˍt(t)), ModelingToolkit.Equation(derivative(x(t), t), xˍt(t))], t, ModelingToolkit.V
ariable[uˍtt, xˍt, uˍt, u, x], ModelingToolkit.Variable[], Base.RefValue{Array{ModelingToolkit.Expression,1}}(ModelingToolkit.Expr
ession[]), Base.RefValue{Array{ModelingToolkit.Expression,2}}(Array{ModelingToolkit.Expression}(undef,0,0)), Base.RefValue{Array{M
odelingToolkit.Expression,2}}(Array{ModelingToolkit.Expression}(undef,0,0)), Base.RefValue{Array{ModelingToolkit.Expression,2}}(Ar
ray{ModelingToolkit.Expression}(undef,0,0)), Symbol("##ODESystem#378"), ModelingToolkit.ODESystem[])
````



````julia
de1.eqs
````


````
5-element Array{ModelingToolkit.Equation,1}:
 ModelingToolkit.Equation(derivative(uˍtt(t), t), ((2 * uˍtt(t) + uˍt(t)) + xˍt(t)) + 1)
 ModelingToolkit.Equation(derivative(xˍt(t), t), xˍt(t) + 2)
 ModelingToolkit.Equation(derivative(uˍt(t), t), uˍtt(t))
 ModelingToolkit.Equation(derivative(u(t), t), uˍt(t))
 ModelingToolkit.Equation(derivative(x(t), t), xˍt(t))
````





This has generated a system of 5 first order ODE systems which can now be used in the ODE solvers.

### Linear Algebra... for free?

Let's take a look at how to extend ModelingToolkit.jl in new directions. Let's define a Jacobian just by using the derivative primatives by hand:

````julia
@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t Dx'~x Dy'~y Dz'~z
eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]
J = [Dx(eqs[1].rhs) Dy(eqs[1].rhs) Dz(eqs[1].rhs)
 Dx(eqs[2].rhs) Dy(eqs[2].rhs) Dz(eqs[2].rhs)
 Dx(eqs[3].rhs) Dy(eqs[3].rhs) Dz(eqs[3].rhs)]
````


````
3×3 Array{ModelingToolkit.Operation,2}:
        derivative(σ * (y(t) - x(t)), x(t))  …         derivative(σ * (y(t) - x(t)), z(t))
 derivative(x(t) * (ρ - z(t)) - y(t), x(t))     derivative(x(t) * (ρ - z(t)) - y(t), z(t))
   derivative(x(t) * y(t) - β * z(t), x(t))       derivative(x(t) * y(t) - β * z(t), z(t))
````





Notice that this writes the derivatives in a "lazy" manner. If we want to actually compute the derivatives, we can expand out those expressions:

````julia
J = expand_derivatives.(J)
````


````
3×3 Array{ModelingToolkit.Expression,2}:
   σ * -1             σ  Constant(0)
 ρ - z(t)  Constant(-1)    x(t) * -1
     y(t)          x(t)          -1β
````





Here's the magic of ModelingToolkit.jl: **Julia treats ModelingToolkit expressions like a Number, and so generic numerical functions are directly usable on ModelingToolkit expressions!** Let's compute the LU-factorization of this Jacobian we defined using Julia's Base linear algebra library.

````julia
using LinearAlgebra
luJ = lu(J,Val(false))
````


````
Error: MethodError: ModelingToolkit.Expression(::ModelingToolkit.Constant) is ambiguous. Candidates:
  (::Type{T})(x::T) where T<:Number in Core at boot.jl:715
  ModelingToolkit.Expression(ex; mod) in ModelingToolkit at /builds/JuliaGPU/DiffEqTutorials.jl/.julia/packages/ModelingToolkit/aK
a4S/src/utils.jl:3
Possible fix, define
  ModelingToolkit.Expression(::ModelingToolkit.Expression)
````



````julia
luJ.L
````


````
Error: UndefVarError: luJ not defined
````





and the inverse?

````julia
invJ = inv(luJ)
````


````
Error: UndefVarError: luJ not defined
````





#### Thus ModelingToolkit.jl can utilize existing numerical code on symbolic codes

Let's follow this thread a little deeper.

### Automatically convert numerical codes to symbolic

Let's take someone's code written to numerically solve the Lorenz equation:

````julia
function lorenz(du,u,p,t)
 du[1] = p[1]*(u[2]-u[1])
 du[2] = u[1]*(p[2]-u[3]) - u[2]
 du[3] = u[1]*u[2] - p[3]*u[3]
end
````


````
lorenz (generic function with 1 method)
````





Since ModelingToolkit can trace generic numerical functions in Julia, let's trace it with Operations. When we do this, it'll spit out a symbolic representation of their numerical code:

````julia
u = [x,y,z]
du = similar(u)
p = [σ,ρ,β]
lorenz(du,u,p,t)
du
````


````
3-element Array{ModelingToolkit.Operation,1}:
        σ * (y(t) - x(t))
 x(t) * (ρ - z(t)) - y(t)
   x(t) * y(t) - β * z(t)
````





We can then perform symbolic manipulations on their numerical code, and build a new numerical code that optimizes/fixes their original function!

````julia
J = [Dx(du[1]) Dy(du[1]) Dz(du[1])
     Dx(du[2]) Dy(du[2]) Dz(du[2])
     Dx(du[3]) Dy(du[3]) Dz(du[3])]
J = expand_derivatives.(J)
````


````
3×3 Array{ModelingToolkit.Expression,2}:
   σ * -1             σ  Constant(0)
 ρ - z(t)  Constant(-1)    x(t) * -1
     y(t)          x(t)          -1β
````





### Automated Sparsity Detection

In many cases one has to speed up large modeling frameworks by taking into account sparsity. While ModelingToolkit.jl can be used to compute Jacobians, we can write a standard Julia function in order to get a spase matrix of expressions which automatically detects and utilizes the sparsity of their function.

````julia
using SparseArrays
function SparseArrays.SparseMatrixCSC(M::Matrix{T}) where {T<:ModelingToolkit.Expression}
    idxs = findall(!iszero, M)
    I = [i[1] for i in idxs]
    J = [i[2] for i in idxs]
    V = [M[i] for i in idxs]
    return SparseArrays.sparse(I, J, V, size(M)...)
end
sJ = SparseMatrixCSC(J)
````


````
3×3 SparseArrays.SparseMatrixCSC{ModelingToolkit.Expression,Int64} with 8 stored entries:
  [1, 1]  =  σ * -1
  [2, 1]  =  ρ - z(t)
  [3, 1]  =  y(t)
  [1, 2]  =  σ
  [2, 2]  =  Constant(-1)
  [3, 2]  =  x(t)
  [2, 3]  =  x(t) * -1
  [3, 3]  =  -1β
````





### Dependent Variables, Functions, Chain Rule

"Variables" are overloaded. When you are solving a differential equation, the variable `u(t)` is actually a function of time. In order to handle these kinds of variables in a mathematically correct and extensible manner, the ModelingToolkit IR actually treats variables as functions, and constant variables are simply 0-ary functions (`t()`).

We can utilize this idea to have parameters that are also functions. For example, we can have a parameter σ which acts as a function of 1 argument, and then utilize this function within our differential equations:

````julia
@parameters σ(..)
eqs = [D(x) ~ σ(t-1)*(y-x),
       D(y) ~ x*(σ(t^2)-z)-y,
       D(z) ~ x*y - β*z]
````


````
3-element Array{ModelingToolkit.Equation,1}:
 ModelingToolkit.Equation(derivative(x(t), t), σ(t - 1) * (y(t) - x(t)))
 ModelingToolkit.Equation(derivative(y(t), t), x(t) * (σ(t ^ 2) - z(t)) - y(t))
 ModelingToolkit.Equation(derivative(z(t), t), x(t) * y(t) - β * z(t))
````





Notice that when we calculate the derivative with respect to `t`, the chain rule is automatically handled:

````julia
@derivatives Dₜ'~t
Dₜ(x*(σ(t^2)-z)-y)
expand_derivatives(Dₜ(x*(σ(t^2)-z)-y))
````


````
(σ(t ^ 2) - z(t)) * derivative(x(t), t) + x(t) * (derivative(σ(t ^ 2), t) + -1 * derivative(z(t), t)) + -1 * derivative(y(t), t)
````





### Hackability: Extend directly from the language

ModelingToolkit.jl is written in Julia, and thus it can be directly extended from Julia itself. Let's define a normal Julia function and call it with a variable:

````julia
_f(x) = 2x + x^2
_f(x)
````


````
2 * x(t) + x(t) ^ 2
````





Recall that when we do that, it will automatically trace this function and then build a symbolic expression. But what if we wanted our function to be a primative in the symbolic framework? This can be done by registering the function.

````julia
f(x) = 2x + x^2
@register f(x)
````


````
f (generic function with 2 methods)
````





Now this function is a new primitive:

````julia
f(x)
````


````
f(x(t))
````





and we can now define derivatives of our function:

````julia
function ModelingToolkit.derivative(::typeof(f), args::NTuple{1,Any}, ::Val{1})
    2 + 2args[1]
end
expand_derivatives(Dx(f(x)))
````


````
2 + 2 * x(t)
````




## Appendix

 This tutorial is part of the DiffEqTutorials.jl repository, found at: <https://github.com/JuliaDiffEq/DiffEqTutorials.jl>

To locally run this tutorial, do the following commands:
```
using DiffEqTutorials
DiffEqTutorials.weave_file("ode_extras","01-ModelingToolkit.jmd")
```

Computer Information:
```
Julia Version 1.4.2
Commit 44fa15b150* (2020-05-23 18:35 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU E5-2603 v4 @ 1.70GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-8.0.1 (ORCJIT, broadwell)
Environment:
  JULIA_DEPOT_PATH = /builds/JuliaGPU/DiffEqTutorials.jl/.julia
  JULIA_CUDA_MEMORY_LIMIT = 536870912
  JULIA_PROJECT = @.
  JULIA_NUM_THREADS = 4

```

Package Information:

```
Status `/builds/JuliaGPU/DiffEqTutorials.jl/Project.toml`
[2169fc97-5a83-5252-b627-83903c6c433c] AlgebraicMultigrid 0.3.0
[7e558dbc-694d-5a72-987c-6f4ebed21442] ArbNumerics 1.0.5
[6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf] BenchmarkTools 0.5.0
[be33ccc6-a3ff-5ff2-a52e-74243cff1e17] CUDAnative 3.1.0
[159f3aea-2a34-519c-b102-8c37f9878175] Cairo 1.0.3
[3a865a2d-5b23-5a0f-bc46-62713ec82fae] CuArrays 2.2.1
[55939f99-70c6-5e9b-8bb0-5071ed7d61fd] DecFP 0.4.10
[abce61dc-4473-55a0-ba07-351d65e31d42] Decimals 0.4.1
[ebbdde9d-f333-5424-9be2-dbf1e9acfb5e] DiffEqBayes 2.15.0
[eb300fae-53e8-50a0-950c-e21f52c2b7e0] DiffEqBiological 4.3.0
[459566f4-90b8-5000-8ac3-15dfb0a30def] DiffEqCallbacks 2.13.3
[f3b72e0c-5b89-59e1-b016-84e28bfd966d] DiffEqDevTools 2.21.0
[77a26b50-5914-5dd7-bc55-306e6241c503] DiffEqNoiseProcess 4.2.0
[9fdde737-9c7f-55bf-ade8-46b3f136cc48] DiffEqOperators 4.10.0
[1130ab10-4a5a-5621-a13d-e4788d82bd4c] DiffEqParamEstim 1.14.1
[055956cb-9e8b-5191-98cc-73ae4a59e68a] DiffEqPhysics 3.2.0
[0c46a032-eb83-5123-abaf-570d42b7fbaa] DifferentialEquations 6.14.0
[31c24e10-a181-5473-b8eb-7969acd0382f] Distributions 0.23.4
[497a8b3b-efae-58df-a0af-a86822472b78] DoubleFloats 1.1.12
[587475ba-b771-5e3f-ad9e-33799f191a9c] Flux 0.10.4
[f6369f11-7733-5829-9624-2563aa707210] ForwardDiff 0.10.10
[7073ff75-c697-5162-941a-fcdaad2a7d2a] IJulia 1.21.2
[23fbe1c1-3f47-55db-b15f-69d7ec21a316] Latexify 0.13.5
[c7f686f2-ff18-58e9-bc7b-31028e88f75d] MCMCChains 3.0.12
[eff96d63-e80a-5855-80a2-b1b0885c5ab7] Measurements 2.2.1
[961ee093-0014-501f-94e3-6117800e7a78] ModelingToolkit 3.1.1
[2774e3e8-f4cf-5e23-947b-6d7e65073b56] NLsolve 4.4.0
[8faf48c0-8b73-11e9-0e63-2155955bfa4d] NeuralNetDiffEq 1.5.0
[429524aa-4258-5aef-a3af-852621145aeb] Optim 0.21.0
[1dea7af3-3e70-54e6-95c3-0bf5283fa5ed] OrdinaryDiffEq 5.41.0
[65888b18-ceab-5e60-b2b9-181511a3b968] ParameterizedFunctions 5.3.0
[91a5bcdd-55d7-5caf-9e0b-520d859cae80] Plots 1.4.0
[d330b81b-6aea-500a-939a-2ce795aea3ee] PyPlot 2.9.0
[731186ca-8d62-57ce-b412-fbd966d074cd] RecursiveArrayTools 2.4.4
[47a9eef4-7e08-11e9-0b38-333d64bd3804] SparseDiffTools 1.8.0
[684fba80-ace3-11e9-3d08-3bc7ed6f96df] SparsityDetection 0.3.2
[90137ffa-7385-5640-81b9-e52037218182] StaticArrays 0.12.3
[f3b207a7-027a-5e70-b257-86293d7955fd] StatsPlots 0.14.6
[789caeaf-c7a9-5a7d-9973-96adeb23e2a0] StochasticDiffEq 6.23.1
[c3572dad-4567-51f8-b174-8c6c989267f4] Sundials 4.2.3
[1986cc42-f94f-5a68-af5c-568840ba703d] Unitful 1.2.1
[44d3d7a6-8a23-5bf8-98c5-b353f8df5ec9] Weave 0.10.2
[b77e0a4c-d291-57a0-90e8-8db25a27a240] InteractiveUtils
[37e2e46d-f89d-539d-b4ee-838fcccc9c8e] LinearAlgebra
[44cfe95a-1eb2-52ea-b672-e2afdf69b78f] Pkg
```
