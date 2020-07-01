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
:((var"##MTKArg#361", var"##MTKArg#362", var"##MTKArg#363")->begin
          if var"##MTKArg#361" isa Array || !(typeof(var"##MTKArg#361") <: StaticArray) && false
              return @inbounds(begin
                          let (x, y, z, σ, ρ, β, t) = (var"##MTKArg#361"[1], var"##MTKArg#361"[2], var"##MTKArg#361"[3], var"##MTK
Arg#362"[1], var"##MTKArg#362"[2], var"##MTKArg#362"[3], var"##MTKArg#363")
                              [(getproperty(Base, :*))(σ, (getproperty(Base, :-))(y, x)), (getproperty(Base, :-))((getproperty(Bas
e, :*))(x, (getproperty(Base, :-))(ρ, z)), y), (getproperty(Base, :-))((getproperty(Base, :*))(x, y), (getproperty(Base, :*))(β, z
))]
                          end
                      end)
          else
              X = @inbounds(begin
                          let (x, y, z, σ, ρ, β, t) = (var"##MTKArg#361"[1], var"##MTKArg#361"[2], var"##MTKArg#361"[3], var"##MTK
Arg#362"[1], var"##MTKArg#362"[2], var"##MTKArg#362"[3], var"##MTKArg#363")
                              ((getproperty(Base, :*))(σ, (getproperty(Base, :-))(y, x)), (getproperty(Base, :-))((getproperty(Bas
e, :*))(x, (getproperty(Base, :-))(ρ, z)), y), (getproperty(Base, :-))((getproperty(Base, :*))(x, y), (getproperty(Base, :*))(β, z
)))
                          end
                      end)
              construct = if var"##MTKArg#361" isa ModelingToolkit.StaticArrays.StaticArray
                      (getproperty(ModelingToolkit.StaticArrays, :similar_type))(typeof(var"##MTKArg#361"), eltype(X))
                  else
                      x->begin
                              convert(typeof(var"##MTKArg#361"), x)
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
:((var"##MTIIPVar#371", var"##MTKArg#367", var"##MTKArg#368", var"##MTKArg#369")->begin
          @inbounds begin
                  let (x, y, z, σ, ρ, β, t) = (var"##MTKArg#367"[1], var"##MTKArg#367"[2], var"##MTKArg#367"[3], var"##MTKArg#368"
[1], var"##MTKArg#368"[2], var"##MTKArg#368"[3], var"##MTKArg#369")
                      var"##MTIIPVar#371"[1] = (getproperty(Base, :*))(σ, (getproperty(Base, :-))(y, x))
                      var"##MTIIPVar#371"[2] = (getproperty(Base, :-))((getproperty(Base, :*))(x, (getproperty(Base, :-))(ρ, z)), 
y)
                      var"##MTIIPVar#371"[3] = (getproperty(Base, :-))((getproperty(Base, :*))(x, y), (getproperty(Base, :*))(β, z
))
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
           -1σ             σ  Constant(0)
 -1 * z(t) + ρ  Constant(-1)    -1 * x(t)
          y(t)          x(t)          -1β
````





It will automatically generate functions for using this Jacobian within the stiff ODE solvers for faster solving:

````julia
jac_expr = generate_jacobian(de)
````


````
(:((var"##MTKArg#373", var"##MTKArg#374", var"##MTKArg#375")->begin
          if var"##MTKArg#373" isa Array || !(typeof(var"##MTKArg#373") <: StaticArray) && false
              return @inbounds(begin
                          let (x, y, z, σ, ρ, β, t) = (var"##MTKArg#373"[1], var"##MTKArg#373"[2], var"##MTKArg#373"[3], var"##MTK
Arg#374"[1], var"##MTKArg#374"[2], var"##MTKArg#374"[3], var"##MTKArg#375")
                              [(getproperty(Base, :*))(-1, σ) σ 0; (getproperty(Base, :+))((getproperty(Base, :*))(-1, z), ρ) -1 (
getproperty(Base, :*))(-1, x); y x (getproperty(Base, :*))(-1, β)]
                          end
                      end)
          else
              X = @inbounds(begin
                          let (x, y, z, σ, ρ, β, t) = (var"##MTKArg#373"[1], var"##MTKArg#373"[2], var"##MTKArg#373"[3], var"##MTK
Arg#374"[1], var"##MTKArg#374"[2], var"##MTKArg#374"[3], var"##MTKArg#375")
                              ((getproperty(Base, :*))(-1, σ), (getproperty(Base, :+))((getproperty(Base, :*))(-1, z), ρ), y, σ, -
1, x, 0, (getproperty(Base, :*))(-1, x), (getproperty(Base, :*))(-1, β))
                          end
                      end)
              construct = if var"##MTKArg#373" isa ModelingToolkit.StaticArrays.StaticArray
                      ModelingToolkit.StaticArrays.SMatrix{3, 3}
                  else
                      x->begin
                              out = similar(typeof(var"##MTKArg#373"), 3, 3)
                              out .= x
                          end
                  end
              return construct(X)
          end
      end), :((var"##MTIIPVar#377", var"##MTKArg#373", var"##MTKArg#374", var"##MTKArg#375")->begin
          @inbounds begin
                  let (x, y, z, σ, ρ, β, t) = (var"##MTKArg#373"[1], var"##MTKArg#373"[2], var"##MTKArg#373"[3], var"##MTKArg#374"
[1], var"##MTKArg#374"[2], var"##MTKArg#374"[3], var"##MTKArg#375")
                      var"##MTIIPVar#377"[1] = (getproperty(Base, :*))(-1, σ)
                      var"##MTIIPVar#377"[2] = (getproperty(Base, :+))((getproperty(Base, :*))(-1, z), ρ)
                      var"##MTIIPVar#377"[3] = y
                      var"##MTIIPVar#377"[4] = σ
                      var"##MTIIPVar#377"[5] = -1
                      var"##MTIIPVar#377"[6] = x
                      var"##MTIIPVar#377"[7] = 0
                      var"##MTIIPVar#377"[8] = (getproperty(Base, :*))(-1, x)
                      var"##MTIIPVar#377"[9] = (getproperty(Base, :*))(-1, β)
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
(:((var"##MTKArg#379", var"##MTKArg#380", var"##MTKArg#381", var"##MTKArg#382")->begin
          if var"##MTKArg#379" isa Array || !(typeof(var"##MTKArg#379") <: StaticArray) && false
              return @inbounds(begin
                          let (x, y, z, σ, ρ, β, __MTKWgamma, t) = (var"##MTKArg#379"[1], var"##MTKArg#379"[2], var"##MTKArg#379"[
3], var"##MTKArg#380"[1], var"##MTKArg#380"[2], var"##MTKArg#380"[3], var"##MTKArg#381", var"##MTKArg#382")
                              [(getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __MTKWgamma, σ)) (getproperty(Base, :*))(__
MTKWgamma, σ) 0; (getproperty(Base, :*))((getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __MTKWg
amma, σ))), __MTKWgamma, (getproperty(Base, :+))((getproperty(Base, :*))(-1, z), ρ)) (getproperty(Base, :*))(-1, (getproperty(Base
, :+))(1, (getproperty(Base, :*))((getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __MTKWgamma, σ
))), __MTKWgamma ^ 2, σ, (getproperty(Base, :+))((getproperty(Base, :*))(-1, z), ρ)), __MTKWgamma)) (getproperty(Base, :*))(-1, x,
 __MTKWgamma); (getproperty(Base, :*))((getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __MTKWgam
ma, σ))), y, __MTKWgamma) (getproperty(Base, :*))((getproperty(Base, :inv))((getproperty(Base, :*))(-1, (getproperty(Base, :+))(1,
 (getproperty(Base, :*))((getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __MTKWgamma, σ))), __MT
KWgamma ^ 2, σ, (getproperty(Base, :+))((getproperty(Base, :*))(-1, z), ρ)), __MTKWgamma))), (getproperty(Base, :+))((getproperty(
Base, :*))(x, __MTKWgamma), (getproperty(Base, :*))(-1, (getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :
*))(-1, __MTKWgamma, σ))), y, __MTKWgamma ^ 2, σ))) (getproperty(Base, :+))((getproperty(Base, :*))(-1, (getproperty(Base, :+))(1,
 (getproperty(Base, :*))(-1, (getproperty(Base, :inv))((getproperty(Base, :*))(-1, (getproperty(Base, :+))(1, (getproperty(Base, :
*))((getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __MTKWgamma, σ))), __MTKWgamma ^ 2, σ, (getp
roperty(Base, :+))((getproperty(Base, :*))(-1, z), ρ)), __MTKWgamma))), x, __MTKWgamma, (getproperty(Base, :+))((getproperty(Base,
 :*))(x, __MTKWgamma), (getproperty(Base, :*))(-1, (getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-
1, __MTKWgamma, σ))), y, __MTKWgamma ^ 2, σ))))), (getproperty(Base, :*))(-1, __MTKWgamma, β))]
                          end
                      end)
          else
              X = @inbounds(begin
                          let (x, y, z, σ, ρ, β, __MTKWgamma, t) = (var"##MTKArg#379"[1], var"##MTKArg#379"[2], var"##MTKArg#379"[
3], var"##MTKArg#380"[1], var"##MTKArg#380"[2], var"##MTKArg#380"[3], var"##MTKArg#381", var"##MTKArg#382")
                              ((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __MTKWgamma, σ)), (getproperty(Base, :*))((
getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __MTKWgamma, σ))), __MTKWgamma, (getproperty(Base
, :+))((getproperty(Base, :*))(-1, z), ρ)), (getproperty(Base, :*))((getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getprop
erty(Base, :*))(-1, __MTKWgamma, σ))), y, __MTKWgamma), (getproperty(Base, :*))(__MTKWgamma, σ), (getproperty(Base, :*))(-1, (getp
roperty(Base, :+))(1, (getproperty(Base, :*))((getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __
MTKWgamma, σ))), __MTKWgamma ^ 2, σ, (getproperty(Base, :+))((getproperty(Base, :*))(-1, z), ρ)), __MTKWgamma)), (getproperty(Base
, :*))((getproperty(Base, :inv))((getproperty(Base, :*))(-1, (getproperty(Base, :+))(1, (getproperty(Base, :*))((getproperty(Base,
 :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __MTKWgamma, σ))), __MTKWgamma ^ 2, σ, (getproperty(Base, :+))((ge
tproperty(Base, :*))(-1, z), ρ)), __MTKWgamma))), (getproperty(Base, :+))((getproperty(Base, :*))(x, __MTKWgamma), (getproperty(Ba
se, :*))(-1, (getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __MTKWgamma, σ))), y, __MTKWgamma ^
 2, σ))), 0, (getproperty(Base, :*))(-1, x, __MTKWgamma), (getproperty(Base, :+))((getproperty(Base, :*))(-1, (getproperty(Base, :
+))(1, (getproperty(Base, :*))(-1, (getproperty(Base, :inv))((getproperty(Base, :*))(-1, (getproperty(Base, :+))(1, (getproperty(B
ase, :*))((getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __MTKWgamma, σ))), __MTKWgamma ^ 2, σ,
 (getproperty(Base, :+))((getproperty(Base, :*))(-1, z), ρ)), __MTKWgamma))), x, __MTKWgamma, (getproperty(Base, :+))((getproperty
(Base, :*))(x, __MTKWgamma), (getproperty(Base, :*))(-1, (getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, 
:*))(-1, __MTKWgamma, σ))), y, __MTKWgamma ^ 2, σ))))), (getproperty(Base, :*))(-1, __MTKWgamma, β)))
                          end
                      end)
              construct = (x->begin
                          A = SMatrix{(3, 3)...}(x)
                          (getproperty(StaticArrays, :LU))(LowerTriangular(SMatrix{(3, 3)...}(UnitLowerTriangular(A))), UpperTrian
gular(A), SVector(ntuple((n->begin
                                              n
                                          end), max((3, 3)...))))
                      end)
              return construct(X)
          end
      end), :((var"##MTIIPVar#384", var"##MTKArg#379", var"##MTKArg#380", var"##MTKArg#381", var"##MTKArg#382")->begin
          @inbounds begin
                  let (x, y, z, σ, ρ, β, __MTKWgamma, t) = (var"##MTKArg#379"[1], var"##MTKArg#379"[2], var"##MTKArg#379"[3], var"
##MTKArg#380"[1], var"##MTKArg#380"[2], var"##MTKArg#380"[3], var"##MTKArg#381", var"##MTKArg#382")
                      var"##MTIIPVar#384"[1] = (getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __MTKWgamma, σ))
                      var"##MTIIPVar#384"[2] = (getproperty(Base, :*))((getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getp
roperty(Base, :*))(-1, __MTKWgamma, σ))), __MTKWgamma, (getproperty(Base, :+))((getproperty(Base, :*))(-1, z), ρ))
                      var"##MTIIPVar#384"[3] = (getproperty(Base, :*))((getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getp
roperty(Base, :*))(-1, __MTKWgamma, σ))), y, __MTKWgamma)
                      var"##MTIIPVar#384"[4] = (getproperty(Base, :*))(__MTKWgamma, σ)
                      var"##MTIIPVar#384"[5] = (getproperty(Base, :*))(-1, (getproperty(Base, :+))(1, (getproperty(Base, :*))((get
property(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __MTKWgamma, σ))), __MTKWgamma ^ 2, σ, (getproperty(
Base, :+))((getproperty(Base, :*))(-1, z), ρ)), __MTKWgamma))
                      var"##MTIIPVar#384"[6] = (getproperty(Base, :*))((getproperty(Base, :inv))((getproperty(Base, :*))(-1, (getp
roperty(Base, :+))(1, (getproperty(Base, :*))((getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __
MTKWgamma, σ))), __MTKWgamma ^ 2, σ, (getproperty(Base, :+))((getproperty(Base, :*))(-1, z), ρ)), __MTKWgamma))), (getproperty(Bas
e, :+))((getproperty(Base, :*))(x, __MTKWgamma), (getproperty(Base, :*))(-1, (getproperty(Base, :inv))((getproperty(Base, :+))(-1,
 (getproperty(Base, :*))(-1, __MTKWgamma, σ))), y, __MTKWgamma ^ 2, σ)))
                      var"##MTIIPVar#384"[7] = 0
                      var"##MTIIPVar#384"[8] = (getproperty(Base, :*))(-1, x, __MTKWgamma)
                      var"##MTIIPVar#384"[9] = (getproperty(Base, :+))((getproperty(Base, :*))(-1, (getproperty(Base, :+))(1, (get
property(Base, :*))(-1, (getproperty(Base, :inv))((getproperty(Base, :*))(-1, (getproperty(Base, :+))(1, (getproperty(Base, :*))((
getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __MTKWgamma, σ))), __MTKWgamma ^ 2, σ, (getproper
ty(Base, :+))((getproperty(Base, :*))(-1, z), ρ)), __MTKWgamma))), x, __MTKWgamma, (getproperty(Base, :+))((getproperty(Base, :*))
(x, __MTKWgamma), (getproperty(Base, :*))(-1, (getproperty(Base, :inv))((getproperty(Base, :+))(-1, (getproperty(Base, :*))(-1, __
MTKWgamma, σ))), y, __MTKWgamma ^ 2, σ))))), (getproperty(Base, :*))(-1, __MTKWgamma, β))
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
(:((var"##MTKArg#394", var"##MTKArg#395")->begin
          if var"##MTKArg#394" isa Array || !(typeof(var"##MTKArg#394") <: StaticArray) && false
              return @inbounds(begin
                          let (x, y, z, σ, ρ, β) = (var"##MTKArg#394"[1], var"##MTKArg#394"[2], var"##MTKArg#394"[3], var"##MTKArg
#395"[1], var"##MTKArg#395"[2], var"##MTKArg#395"[3])
                              [(*)(σ, (-)(y, x)), (-)((*)(x, (-)(ρ, z)), y), (-)((*)(x, y), (*)(β, z))]
                          end
                      end)
          else
              X = @inbounds(begin
                          let (x, y, z, σ, ρ, β) = (var"##MTKArg#394"[1], var"##MTKArg#394"[2], var"##MTKArg#394"[3], var"##MTKArg
#395"[1], var"##MTKArg#395"[2], var"##MTKArg#395"[3])
                              ((*)(σ, (-)(y, x)), (-)((*)(x, (-)(ρ, z)), y), (-)((*)(x, y), (*)(β, z)))
                          end
                      end)
              construct = if var"##MTKArg#394" isa ModelingToolkit.StaticArrays.StaticArray
                      (getproperty(ModelingToolkit.StaticArrays, :similar_type))(typeof(var"##MTKArg#394"), eltype(X))
                  else
                      x->begin
                              convert(typeof(var"##MTKArg#394"), x)
                          end
                  end
              return construct(X)
          end
      end), :((var"##MTIIPVar#397", var"##MTKArg#394", var"##MTKArg#395")->begin
          @inbounds begin
                  let (x, y, z, σ, ρ, β) = (var"##MTKArg#394"[1], var"##MTKArg#394"[2], var"##MTKArg#394"[3], var"##MTKArg#395"[1]
, var"##MTKArg#395"[2], var"##MTKArg#395"[3])
                      var"##MTIIPVar#397"[1] = (*)(σ, (-)(y, x))
                      var"##MTIIPVar#397"[2] = (-)((*)(x, (-)(ρ, z)), y)
                      var"##MTIIPVar#397"[3] = (-)((*)(x, y), (*)(β, z))
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
ray{ModelingToolkit.Expression}(undef,0,0)), Symbol("##ODESystem#400"), ModelingToolkit.ODESystem[])
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
           -1σ             σ  Constant(0)
 -1 * z(t) + ρ  Constant(-1)    -1 * x(t)
          y(t)          x(t)          -1β
````





Here's the magic of ModelingToolkit.jl: **Julia treats ModelingToolkit expressions like a Number, and so generic numerical functions are directly usable on ModelingToolkit expressions!** Let's compute the LU-factorization of this Jacobian we defined using Julia's Base linear algebra library.

````julia
using LinearAlgebra
luJ = lu(J,Val(false))
````


````
LinearAlgebra.LU{ModelingToolkit.Expression,Array{ModelingToolkit.Expression,2}}
L factor:
3×3 Array{ModelingToolkit.Expression,2}:
                Constant(1)  …  Constant(0)
 (-1 * z(t) + ρ) * inv(-1σ)     Constant(0)
            y(t) * inv(-1σ)     Constant(1)
U factor:
3×3 Array{ModelingToolkit.Expression,2}:
         -1σ  …                                                                                                                   
                                  Constant(0)
 Constant(0)                                                                                                                      
 -1 * x(t) - ((-1 * z(t) + ρ) * inv(-1σ)) * 0
 Constant(0)     (-1β - (y(t) * inv(-1σ)) * 0) - ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1 - ((-1 * z(t) + ρ) * inv(-1σ)) * σ)) * (
-1 * x(t) - ((-1 * z(t) + ρ) * inv(-1σ)) * 0)
````



````julia
luJ.L
````


````
3×3 Array{ModelingToolkit.Expression,2}:
                Constant(1)  …  Constant(0)
 (-1 * z(t) + ρ) * inv(-1σ)     Constant(0)
            y(t) * inv(-1σ)     Constant(1)
````





and the inverse?

````julia
invJ = inv(luJ)
````


````
3×3 Array{ModelingToolkit.Expression,2}:
 (-1σ) \ ((true - 0 * (((-1β - (y(t) * inv(-1σ)) * 0) - ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1 - ((-1 * z(t) + ρ) * inv(-1σ)) * 
σ)) * (-1 * x(t) - ((-1 * z(t) + ρ) * inv(-1σ)) * 0)) \ ((0 - (y(t) * inv(-1σ)) * true) - ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1
 - ((-1 * z(t) + ρ) * inv(-1σ)) * σ)) * (0 - ((-1 * z(t) + ρ) * inv(-1σ)) * true)))) - σ * ((-1 - ((-1 * z(t) + ρ) * inv(-1σ)) * σ
) \ ((0 - ((-1 * z(t) + ρ) * inv(-1σ)) * true) - (-1 * x(t) - ((-1 * z(t) + ρ) * inv(-1σ)) * 0) * (((-1β - (y(t) * inv(-1σ)) * 0) 
- ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1 - ((-1 * z(t) + ρ) * inv(-1σ)) * σ)) * (-1 * x(t) - ((-1 * z(t) + ρ) * inv(-1σ)) * 0)) 
\ ((0 - (y(t) * inv(-1σ)) * true) - ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1 - ((-1 * z(t) + ρ) * inv(-1σ)) * σ)) * (0 - ((-1 * z(
t) + ρ) * inv(-1σ)) * true))))))  …  (-1σ) \ ((0 - 0 * (((-1β - (y(t) * inv(-1σ)) * 0) - ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1 
- ((-1 * z(t) + ρ) * inv(-1σ)) * σ)) * (-1 * x(t) - ((-1 * z(t) + ρ) * inv(-1σ)) * 0)) \ ((true - (y(t) * inv(-1σ)) * 0) - ((x(t) 
- (y(t) * inv(-1σ)) * σ) * inv(-1 - ((-1 * z(t) + ρ) * inv(-1σ)) * σ)) * (0 - ((-1 * z(t) + ρ) * inv(-1σ)) * 0)))) - σ * ((-1 - ((
-1 * z(t) + ρ) * inv(-1σ)) * σ) \ ((0 - ((-1 * z(t) + ρ) * inv(-1σ)) * 0) - (-1 * x(t) - ((-1 * z(t) + ρ) * inv(-1σ)) * 0) * (((-1
β - (y(t) * inv(-1σ)) * 0) - ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1 - ((-1 * z(t) + ρ) * inv(-1σ)) * σ)) * (-1 * x(t) - ((-1 * z
(t) + ρ) * inv(-1σ)) * 0)) \ ((true - (y(t) * inv(-1σ)) * 0) - ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1 - ((-1 * z(t) + ρ) * inv(-
1σ)) * σ)) * (0 - ((-1 * z(t) + ρ) * inv(-1σ)) * 0))))))
                                                                                                                                  
                                                                                                                                  
                                                                                              (-1 - ((-1 * z(t) + ρ) * inv(-1σ)) *
 σ) \ ((0 - ((-1 * z(t) + ρ) * inv(-1σ)) * true) - (-1 * x(t) - ((-1 * z(t) + ρ) * inv(-1σ)) * 0) * (((-1β - (y(t) * inv(-1σ)) * 0
) - ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1 - ((-1 * z(t) + ρ) * inv(-1σ)) * σ)) * (-1 * x(t) - ((-1 * z(t) + ρ) * inv(-1σ)) * 0)
) \ ((0 - (y(t) * inv(-1σ)) * true) - ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1 - ((-1 * z(t) + ρ) * inv(-1σ)) * σ)) * (0 - ((-1 * 
z(t) + ρ) * inv(-1σ)) * true))))                                                                                                  
                                                                                                                                  
                                                                                                                            (-1 - 
((-1 * z(t) + ρ) * inv(-1σ)) * σ) \ ((0 - ((-1 * z(t) + ρ) * inv(-1σ)) * 0) - (-1 * x(t) - ((-1 * z(t) + ρ) * inv(-1σ)) * 0) * (((
-1β - (y(t) * inv(-1σ)) * 0) - ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1 - ((-1 * z(t) + ρ) * inv(-1σ)) * σ)) * (-1 * x(t) - ((-1 *
 z(t) + ρ) * inv(-1σ)) * 0)) \ ((true - (y(t) * inv(-1σ)) * 0) - ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1 - ((-1 * z(t) + ρ) * inv
(-1σ)) * σ)) * (0 - ((-1 * z(t) + ρ) * inv(-1σ)) * 0))))
                                                                                                                                  
                                                                                                                                  
                                                                                                                                  
                                                                                                       ((-1β - (y(t) * inv(-1σ)) *
 0) - ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1 - ((-1 * z(t) + ρ) * inv(-1σ)) * σ)) * (-1 * x(t) - ((-1 * z(t) + ρ) * inv(-1σ)) * 
0)) \ ((0 - (y(t) * inv(-1σ)) * true) - ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1 - ((-1 * z(t) + ρ) * inv(-1σ)) * σ)) * (0 - ((-1 
* z(t) + ρ) * inv(-1σ)) * true))                                                                                                  
                                                                                                                                  
                                                                                                                                  
                                                                                                                                  
((-1β - (y(t) * inv(-1σ)) * 0) - ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1 - ((-1 * z(t) + ρ) * inv(-1σ)) * σ)) * (-1 * x(t) - ((-1
 * z(t) + ρ) * inv(-1σ)) * 0)) \ ((true - (y(t) * inv(-1σ)) * 0) - ((x(t) - (y(t) * inv(-1σ)) * σ) * inv(-1 - ((-1 * z(t) + ρ) * i
nv(-1σ)) * σ)) * (0 - ((-1 * z(t) + ρ) * inv(-1σ)) * 0))
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
           -1σ             σ  Constant(0)
 -1 * z(t) + ρ  Constant(-1)    -1 * x(t)
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
  [1, 1]  =  -1σ
  [2, 1]  =  -1 * z(t) + ρ
  [3, 1]  =  y(t)
  [1, 2]  =  σ
  [2, 2]  =  Constant(-1)
  [3, 2]  =  x(t)
  [2, 3]  =  -1 * x(t)
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
derivative(x(t), t) * (σ(t ^ 2) + -1 * z(t)) + -1 * derivative(y(t), t) + x(t) * (derivative(σ(t ^ 2), t) + -1 * derivative(z(t), 
t))
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
  CPU: Intel(R) Core(TM) i7-9700K CPU @ 3.60GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-8.0.1 (ORCJIT, skylake)
Environment:
  JULIA_DEPOT_PATH = /builds/JuliaGPU/DiffEqTutorials.jl/.julia
  JULIA_CUDA_MEMORY_LIMIT = 536870912
  JULIA_PROJECT = @.
  JULIA_NUM_THREADS = 4

```

Package Information:

```
Status `/builds/JuliaGPU/DiffEqTutorials.jl/tutorials/ode_extras/Project.toml`
[f3b72e0c-5b89-59e1-b016-84e28bfd966d] DiffEqDevTools 2.22.0
[961ee093-0014-501f-94e3-6117800e7a78] ModelingToolkit 3.11.0
[76087f3c-5699-56af-9a33-bf431cd00edd] NLopt 0.6.0
[2774e3e8-f4cf-5e23-947b-6d7e65073b56] NLsolve 4.4.0
[429524aa-4258-5aef-a3af-852621145aeb] Optim 0.22.0
[1dea7af3-3e70-54e6-95c3-0bf5283fa5ed] OrdinaryDiffEq 5.41.0
[91a5bcdd-55d7-5caf-9e0b-520d859cae80] Plots 1.5.1
[37e2e46d-f89d-539d-b4ee-838fcccc9c8e] LinearAlgebra
[2f01184e-e22b-5df5-ae63-d93ebab69eaf] SparseArrays
```
