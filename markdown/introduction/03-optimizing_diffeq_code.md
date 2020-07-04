---
author: "Chris Rackauckas"
title: "Optimizing DiffEq Code"
---


In this notebook we will walk through some of the main tools for optimizing your code in order to efficiently solve DifferentialEquations.jl. User-side optimizations are important because, for sufficiently difficult problems, most of the time will be spent inside of your `f` function, the function you are trying to solve. "Efficient" integrators are those that reduce the required number of `f` calls to hit the error tolerance. The main ideas for optimizing your DiffEq code, or any Julia function, are the following:

- Make it non-allocating
- Use StaticArrays for small arrays
- Use broadcast fusion
- Make it type-stable
- Reduce redundant calculations
- Make use of BLAS calls
- Optimize algorithm choice

We'll discuss these strategies in the context of small and large systems. Let's start with small systems.

## Optimizing Small Systems (<100 DEs)

Let's take the classic Lorenz system from before. Let's start by naively writing the system in its out-of-place form:

````julia
function lorenz(u,p,t)
 dx = 10.0*(u[2]-u[1])
 dy = u[1]*(28.0-u[3]) - u[2]
 dz = u[1]*u[2] - (8/3)*u[3]
 [dx,dy,dz]
end
````


````
lorenz (generic function with 1 method)
````





Here, `lorenz` returns an object, `[dx,dy,dz]`, which is created within the body of `lorenz`.

This is a common code pattern from high-level languages like MATLAB, SciPy, or R's deSolve. However, the issue with this form is that it allocates a vector, `[dx,dy,dz]`, at each step. Let's benchmark the solution process with this choice of function:

````julia
using DifferentialEquations, BenchmarkTools
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz,u0,tspan)
@benchmark solve(prob,Tsit5())
````


````
BenchmarkTools.Trial: 
  memory estimate:  10.83 MiB
  allocs estimate:  101435
  --------------
  minimum time:     3.330 ms (0.00% GC)
  median time:      3.378 ms (0.00% GC)
  mean time:        4.411 ms (17.67% GC)
  maximum time:     11.071 ms (49.26% GC)
  --------------
  samples:          1133
  evals/sample:     1
````





The BenchmarkTools package's `@benchmark` runs the code multiple times to get an accurate measurement. The minimum time is the time it takes when your OS and other background processes aren't getting in the way. Notice that in this case it takes about 5ms to solve and allocates around 11.11 MiB. However, if we were to use this inside of a real user code we'd see a lot of time spent doing garbage collection (GC) to clean up all of the arrays we made. Even if we turn off saving we have these allocations.

````julia
@benchmark solve(prob,Tsit5(),save_everystep=false)
````


````
BenchmarkTools.Trial: 
  memory estimate:  9.47 MiB
  allocs estimate:  88655
  --------------
  minimum time:     2.918 ms (0.00% GC)
  median time:      2.953 ms (0.00% GC)
  mean time:        3.995 ms (16.86% GC)
  maximum time:     10.556 ms (49.23% GC)
  --------------
  samples:          1251
  evals/sample:     1
````





The problem of course is that arrays are created every time our derivative function is called. This function is called multiple times per step and is thus the main source of memory usage. To fix this, we can use the in-place form to ***make our code non-allocating***:

````julia
function lorenz!(du,u,p,t)
 du[1] = 10.0*(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end
````


````
lorenz! (generic function with 1 method)
````





Here, instead of creating an array each time, we utilized the cache array `du`. When the inplace form is used, DifferentialEquations.jl takes a different internal route that minimizes the internal allocations as well. When we benchmark this function, we will see quite a difference.

````julia
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz!,u0,tspan)
@benchmark solve(prob,Tsit5())
````


````
BenchmarkTools.Trial: 
  memory estimate:  1.38 MiB
  allocs estimate:  12997
  --------------
  minimum time:     777.483 μs (0.00% GC)
  median time:      785.697 μs (0.00% GC)
  mean time:        913.323 μs (10.49% GC)
  maximum time:     6.689 ms (83.20% GC)
  --------------
  samples:          5463
  evals/sample:     1
````



````julia
@benchmark solve(prob,Tsit5(),save_everystep=false)
````


````
BenchmarkTools.Trial: 
  memory estimate:  5.27 KiB
  allocs estimate:  57
  --------------
  minimum time:     344.549 μs (0.00% GC)
  median time:      347.849 μs (0.00% GC)
  mean time:        348.663 μs (0.00% GC)
  maximum time:     382.426 μs (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1
````





There is a 4x time difference just from that change! Notice there are still some allocations and this is due to the construction of the integration cache. But this doesn't scale with the problem size:

````julia
tspan = (0.0,500.0) # 5x longer than before
prob = ODEProblem(lorenz!,u0,tspan)
@benchmark solve(prob,Tsit5(),save_everystep=false)
````


````
BenchmarkTools.Trial: 
  memory estimate:  5.27 KiB
  allocs estimate:  57
  --------------
  minimum time:     1.729 ms (0.00% GC)
  median time:      1.737 ms (0.00% GC)
  mean time:        1.741 ms (0.00% GC)
  maximum time:     1.934 ms (0.00% GC)
  --------------
  samples:          2870
  evals/sample:     1
````





since that's all just setup allocations.

#### But if the system is small we can optimize even more.

Allocations are only expensive if they are "heap allocations". For a more in-depth definition of heap allocations, [there are a lot of sources online](http://net-informations.com/faq/net/stack-heap.htm). But a good working definition is that heap allocations are variable-sized slabs of memory which have to be pointed to, and this pointer indirection costs time. Additionally, the heap has to be managed and the garbage controllers has to actively keep track of what's on the heap.

However, there's an alternative to heap allocations, known as stack allocations. The stack is statically-sized (known at compile time) and thus its accesses are quick. Additionally, the exact block of memory is known in advance by the compiler, and thus re-using the memory is cheap. This means that allocating on the stack has essentially no cost!

Arrays have to be heap allocated because their size (and thus the amount of memory they take up) is determined at runtime. But there are structures in Julia which are stack-allocated. `struct`s for example are stack-allocated "value-type"s. `Tuple`s are a stack-allocated collection. The most useful data structure for DiffEq though is the `StaticArray` from the package [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl). These arrays have their length determined at compile-time. They are created using macros attached to normal array expressions, for example:

````julia
using StaticArrays
A = @SVector [2.0,3.0,5.0]
````


````
3-element StaticArrays.SArray{Tuple{3},Float64,1,3} with indices SOneTo(3):
 2.0
 3.0
 5.0
````





Notice that the `3` after `SVector` gives the size of the `SVector`. It cannot be changed. Additionally, `SVector`s are immutable, so we have to create a new `SVector` to change values. But remember, we don't have to worry about allocations because this data structure is stack-allocated. `SArray`s have a lot of extra optimizations as well: they have fast matrix multiplication, fast QR factorizations, etc. which directly make use of the information about the size of the array. Thus, when possible they should be used.

Unfortunately static arrays can only be used for sufficiently small arrays. After a certain size, they are forced to heap allocate after some instructions and their compile time balloons. Thus static arrays shouldn't be used if your system has more than 100 variables. Additionally, only the native Julia algorithms can fully utilize static arrays.

Let's ***optimize `lorenz` using static arrays***. Note that in this case, we want to use the out-of-place allocating form, but this time we want to output a static array:

````julia
function lorenz_static(u,p,t)
 dx = 10.0*(u[2]-u[1])
 dy = u[1]*(28.0-u[3]) - u[2]
 dz = u[1]*u[2] - (8/3)*u[3]
 @SVector [dx,dy,dz]
end
````


````
lorenz_static (generic function with 1 method)
````





To make the solver internally use static arrays, we simply give it a static array as the initial condition:

````julia
u0 = @SVector [1.0,0.0,0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz_static,u0,tspan)
@benchmark solve(prob,Tsit5())
````


````
BenchmarkTools.Trial: 
  memory estimate:  466.30 KiB
  allocs estimate:  2596
  --------------
  minimum time:     336.320 μs (0.00% GC)
  median time:      342.997 μs (0.00% GC)
  mean time:        380.309 μs (6.07% GC)
  maximum time:     4.346 ms (86.60% GC)
  --------------
  samples:          10000
  evals/sample:     1
````



````julia
@benchmark solve(prob,Tsit5(),save_everystep=false)
````


````
BenchmarkTools.Trial: 
  memory estimate:  3.44 KiB
  allocs estimate:  32
  --------------
  minimum time:     235.919 μs (0.00% GC)
  median time:      240.566 μs (0.00% GC)
  mean time:        240.749 μs (0.00% GC)
  maximum time:     280.303 μs (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1
````





And that's pretty much all there is to it. With static arrays you don't have to worry about allocating, so use operations like `*` and don't worry about fusing operations (discussed in the next section). Do "the vectorized code" of R/MATLAB/Python and your code in this case will be fast, or directly use the numbers/values.

#### Exercise 1

Implement the out-of-place array, in-place array, and out-of-place static array forms for the [Henon-Heiles System](https://en.wikipedia.org/wiki/H%C3%A9non%E2%80%93Heiles_system) and time the results.

## Optimizing Large Systems

### Interlude: Managing Allocations with Broadcast Fusion

When your system is sufficiently large, or you have to make use of a non-native Julia algorithm, you have to make use of `Array`s. In order to use arrays in the most efficient manner, you need to be careful about temporary allocations. Vectorized calculations naturally have plenty of temporary array allocations. This is because a vectorized calculation outputs a vector. Thus:

````julia
A = rand(1000,1000); B = rand(1000,1000); C = rand(1000,1000)
test(A,B,C) = A + B + C
@benchmark test(A,B,C)
````


````
BenchmarkTools.Trial: 
  memory estimate:  7.63 MiB
  allocs estimate:  2
  --------------
  minimum time:     1.124 ms (0.00% GC)
  median time:      1.146 ms (0.00% GC)
  mean time:        1.313 ms (12.47% GC)
  maximum time:     3.252 ms (60.68% GC)
  --------------
  samples:          3791
  evals/sample:     1
````




That expression `A + B + C` creates 2 arrays. It first creates one for the output of `A + B`, then uses that result array to `+ C` to get the final result. 2 arrays! We don't want that! The first thing to do to fix this is to use broadcast fusion. [Broadcast fusion](https://julialang.org/blog/2017/01/moredots) puts expressions together. For example, instead of doing the `+` operations separately, if we were to add them all at the same time, then we would only have a single array that's created. For example:

````julia
test2(A,B,C) = map((a,b,c)->a+b+c,A,B,C)
@benchmark test2(A,B,C)
````


````
BenchmarkTools.Trial: 
  memory estimate:  7.63 MiB
  allocs estimate:  8
  --------------
  minimum time:     2.182 ms (0.00% GC)
  median time:      2.207 ms (0.00% GC)
  mean time:        2.360 ms (6.40% GC)
  maximum time:     4.155 ms (46.23% GC)
  --------------
  samples:          2116
  evals/sample:     1
````





Puts the whole expression into a single function call, and thus only one array is required to store output. This is the same as writing the loop:

````julia
function test3(A,B,C)
    D = similar(A)
    @inbounds for i in eachindex(A)
        D[i] = A[i] + B[i] + C[i]
    end
    D
end
@benchmark test3(A,B,C)
````


````
BenchmarkTools.Trial: 
  memory estimate:  7.63 MiB
  allocs estimate:  2
  --------------
  minimum time:     1.119 ms (0.00% GC)
  median time:      1.136 ms (0.00% GC)
  mean time:        1.304 ms (12.58% GC)
  maximum time:     3.225 ms (61.21% GC)
  --------------
  samples:          3816
  evals/sample:     1
````





However, Julia's broadcast is syntactic sugar for this. If multiple expressions have a `.`, then it will put those vectorized operations together. Thus:

````julia
test4(A,B,C) = A .+ B .+ C
@benchmark test4(A,B,C)
````


````
BenchmarkTools.Trial: 
  memory estimate:  7.63 MiB
  allocs estimate:  2
  --------------
  minimum time:     1.123 ms (0.00% GC)
  median time:      1.142 ms (0.00% GC)
  mean time:        1.309 ms (12.53% GC)
  maximum time:     3.245 ms (60.90% GC)
  --------------
  samples:          3800
  evals/sample:     1
````





is a version with only 1 array created (the output). Note that `.`s can be used with function calls as well:

````julia
sin.(A) .+ sin.(B)
````


````
1000×1000 Array{Float64,2}:
 1.29144   0.466547  0.23388   0.937569  …  0.971028  0.88962   1.15988
 0.854991  1.06985   1.22783   1.28579      0.698313  0.243468  0.82299
 0.522012  1.4456    0.858349  1.0117       0.913442  0.633198  0.783146
 0.902455  1.21767   1.08965   0.947601     0.762651  0.415332  0.3836
 0.618146  1.34909   0.437813  0.848666     1.2956    0.550397  1.32508
 1.2157    0.935543  1.49661   1.16645   …  0.539873  0.277857  0.244978
 0.881675  1.29354   0.91445   1.20212      0.402342  0.361929  1.38935
 1.43355   1.37723   0.749495  0.14644      0.928398  1.24622   0.692594
 1.08933   1.02454   1.05411   0.993281     0.151414  0.21478   0.910877
 0.496662  0.352621  0.805208  0.616445     0.783966  0.579161  0.670726
 ⋮                                       ⋱                      
 0.604933  0.897257  0.905171  0.548894     0.867156  0.679873  0.510797
 1.06376   1.31646   0.394503  0.292554     0.19358   0.950823  1.00226
 0.816534  0.987653  0.614967  0.870429     1.15783   0.48936   0.987409
 1.56149   1.3981    0.532621  1.4926       1.15511   1.4547    1.13767
 0.978015  0.845137  0.863285  1.30935   …  1.34296   0.987646  0.297642
 0.999255  0.369314  0.784103  0.975572     0.293934  0.673055  0.771749
 0.814539  1.20448   0.99619   0.838672     0.406095  0.924025  0.766225
 0.924922  1.44719   1.19294   0.272792     0.304482  0.568076  1.22903
 1.13388   1.18528   0.562323  1.1867       0.723689  0.9604    0.839514
````





Also, the `@.` macro applys a dot to every operator:

````julia
test5(A,B,C) = @. A + B + C #only one array allocated
@benchmark test5(A,B,C)
````


````
BenchmarkTools.Trial: 
  memory estimate:  7.63 MiB
  allocs estimate:  2
  --------------
  minimum time:     1.126 ms (0.00% GC)
  median time:      1.145 ms (0.00% GC)
  mean time:        1.313 ms (12.52% GC)
  maximum time:     3.270 ms (60.56% GC)
  --------------
  samples:          3789
  evals/sample:     1
````





Using these tools we can get rid of our intermediate array allocations for many vectorized function calls. But we are still allocating the output array. To get rid of that allocation, we can instead use mutation. Mutating broadcast is done via `.=`. For example, if we pre-allocate the output:

````julia
D = zeros(1000,1000);
````


````
1000×1000 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 ⋮                        ⋮              ⋱            ⋮                   
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0
.0
````





Then we can keep re-using this cache for subsequent calculations. The mutating broadcasting form is:

````julia
test6!(D,A,B,C) = D .= A .+ B .+ C #only one array allocated
@benchmark test6!(D,A,B,C)
````


````
BenchmarkTools.Trial: 
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     1.040 ms (0.00% GC)
  median time:      1.049 ms (0.00% GC)
  mean time:        1.051 ms (0.00% GC)
  maximum time:     1.324 ms (0.00% GC)
  --------------
  samples:          4729
  evals/sample:     1
````





If we use `@.` before the `=`, then it will turn it into `.=`:

````julia
test7!(D,A,B,C) = @. D = A + B + C #only one array allocated
@benchmark test7!(D,A,B,C)
````


````
BenchmarkTools.Trial: 
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     1.042 ms (0.00% GC)
  median time:      1.050 ms (0.00% GC)
  mean time:        1.053 ms (0.00% GC)
  maximum time:     1.326 ms (0.00% GC)
  --------------
  samples:          4721
  evals/sample:     1
````





Notice that in this case, there is no "output", and instead the values inside of `D` are what are changed (like with the DiffEq inplace function). Many Julia functions have a mutating form which is denoted with a `!`. For example, the mutating form of the `map` is `map!`:

````julia
test8!(D,A,B,C) = map!((a,b,c)->a+b+c,D,A,B,C)
@benchmark test8!(D,A,B,C)
````


````
BenchmarkTools.Trial: 
  memory estimate:  32 bytes
  allocs estimate:  1
  --------------
  minimum time:     2.306 ms (0.00% GC)
  median time:      2.319 ms (0.00% GC)
  mean time:        2.326 ms (0.00% GC)
  maximum time:     2.551 ms (0.00% GC)
  --------------
  samples:          2146
  evals/sample:     1
````





Some operations require using an alternate mutating form in order to be fast. For example, matrix multiplication via `*` allocates a temporary:

````julia
@benchmark A*B
````


````
BenchmarkTools.Trial: 
  memory estimate:  7.63 MiB
  allocs estimate:  2
  --------------
  minimum time:     6.330 ms (0.00% GC)
  median time:      6.456 ms (0.00% GC)
  mean time:        6.647 ms (2.35% GC)
  maximum time:     8.672 ms (23.39% GC)
  --------------
  samples:          752
  evals/sample:     1
````





Instead, we can use the mutating form `mul!` into a cache array to avoid allocating the output:

````julia
using LinearAlgebra
@benchmark mul!(D,A,B) # same as D = A * B
````


````
BenchmarkTools.Trial: 
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     5.963 ms (0.00% GC)
  median time:      6.003 ms (0.00% GC)
  mean time:        6.013 ms (0.00% GC)
  maximum time:     9.758 ms (0.00% GC)
  --------------
  samples:          831
  evals/sample:     1
````





For repeated calculations this reduced allocation can stop GC cycles and thus lead to more efficient code. Additionally, ***we can fuse together higher level linear algebra operations using BLAS***. The package [SugarBLAS.jl](https://github.com/lopezm94/SugarBLAS.jl) makes it easy to write higher level operations like `alpha*B*A + beta*C` as mutating BLAS calls.

### Example Optimization: Gierer-Meinhardt Reaction-Diffusion PDE Discretization

Let's optimize the solution of a Reaction-Diffusion PDE's discretization. In its discretized form, this is the ODE:

$$
\begin{align}
du &= D_1 (A_y u + u A_x) + \frac{au^2}{v} + \bar{u} - \alpha u\\
dv &= D_2 (A_y v + v A_x) + a u^2 + \beta v
\end{align}
$$

where $u$, $v$, and $A$ are matrices. Here, we will use the simplified version where $A$ is the tridiagonal stencil $[1,-2,1]$, i.e. it's the 2D discretization of the LaPlacian. The native code would be something along the lines of:

````julia
# Generate the constants
p = (1.0,1.0,1.0,10.0,0.001,100.0) # a,α,ubar,β,D1,D2
N = 100
Ax = Array(Tridiagonal([1.0 for i in 1:N-1],[-2.0 for i in 1:N],[1.0 for i in 1:N-1]))
Ay = copy(Ax)
Ax[2,1] = 2.0
Ax[end-1,end] = 2.0
Ay[1,2] = 2.0
Ay[end,end-1] = 2.0

function basic_version!(dr,r,p,t)
  a,α,ubar,β,D1,D2 = p
  u = r[:,:,1]
  v = r[:,:,2]
  Du = D1*(Ay*u + u*Ax)
  Dv = D2*(Ay*v + v*Ax)
  dr[:,:,1] = Du .+ a.*u.*u./v .+ ubar .- α*u
  dr[:,:,2] = Dv .+ a.*u.*u .- β*v
end

a,α,ubar,β,D1,D2 = p
uss = (ubar+β)/α
vss = (a/β)*uss^2
r0 = zeros(100,100,2)
r0[:,:,1] .= uss.+0.1.*rand.()
r0[:,:,2] .= vss

prob = ODEProblem(basic_version!,r0,(0.0,0.1),p)
````


````
ODEProblem with uType Array{Float64,3} and tType Float64. In-place: true
timespan: (0.0, 0.1)
u0: [11.055761709400883 11.024247930398587 … 11.014567891656775 11.02235927
9364153; 11.084473440945578 11.05875630087183 … 11.009734595896523 11.00007
2375751836; … ; 11.013784472099452 11.034456431973112 … 11.060212612466495 
11.039922075072242; 11.083391484651129 11.07991404732281 … 11.0967350008302
45 11.041394909310073]

[12.100000000000001 12.100000000000001 … 12.100000000000001 12.100000000000
001; 12.100000000000001 12.100000000000001 … 12.100000000000001 12.10000000
0000001; … ; 12.100000000000001 12.100000000000001 … 12.100000000000001 12.
100000000000001; 12.100000000000001 12.100000000000001 … 12.100000000000001
 12.100000000000001]
````





In this version we have encoded our initial condition to be a 3-dimensional array, with `u[:,:,1]` being the `A` part and `u[:,:,2]` being the `B` part.

````julia
@benchmark solve(prob,Tsit5())
````


````
BenchmarkTools.Trial: 
  memory estimate:  186.88 MiB
  allocs estimate:  8552
  --------------
  minimum time:     48.590 ms (4.94% GC)
  median time:      51.207 ms (9.27% GC)
  mean time:        50.505 ms (8.08% GC)
  maximum time:     51.454 ms (9.43% GC)
  --------------
  samples:          99
  evals/sample:     1
````





While this version isn't very efficient,

#### We recommend writing the "high-level" code first, and iteratively optimizing it!

The first thing that we can do is get rid of the slicing allocations. The operation `r[:,:,1]` creates a temporary array instead of a "view", i.e. a pointer to the already existing memory. To make it a view, add `@view`. Note that we have to be careful with views because they point to the same memory, and thus changing a view changes the original values:

````julia
A = rand(4)
@show A
````


````
A = [0.3656247108087518, 0.3933396785440917, 0.5893584038255757, 0.95140149
02416317]
````



````julia
B = @view A[1:3]
B[2] = 2
@show A
````


````
A = [0.3656247108087518, 2.0, 0.5893584038255757, 0.9514014902416317]
4-element Array{Float64,1}:
 0.3656247108087518
 2.0
 0.5893584038255757
 0.9514014902416317
````





Notice that changing `B` changed `A`. This is something to be careful of, but at the same time we want to use this since we want to modify the output `dr`. Additionally, the last statement is a purely element-wise operation, and thus we can make use of broadcast fusion there. Let's rewrite `basic_version!` to ***avoid slicing allocations*** and to ***use broadcast fusion***:

````julia
function gm2!(dr,r,p,t)
  a,α,ubar,β,D1,D2 = p
  u = @view r[:,:,1]
  v = @view r[:,:,2]
  du = @view dr[:,:,1]
  dv = @view dr[:,:,2]
  Du = D1*(Ay*u + u*Ax)
  Dv = D2*(Ay*v + v*Ax)
  @. du = Du + a.*u.*u./v + ubar - α*u
  @. dv = Dv + a.*u.*u - β*v
end
prob = ODEProblem(gm2!,r0,(0.0,0.1),p)
@benchmark solve(prob,Tsit5())
````


````
BenchmarkTools.Trial: 
  memory estimate:  119.55 MiB
  allocs estimate:  7082
  --------------
  minimum time:     38.816 ms (6.24% GC)
  median time:      38.999 ms (6.22% GC)
  mean time:        39.408 ms (7.16% GC)
  maximum time:     41.605 ms (11.54% GC)
  --------------
  samples:          127
  evals/sample:     1
````





Now, most of the allocations are taking place in `Du = D1*(Ay*u + u*Ax)` since those operations are vectorized and not mutating. We should instead replace the matrix multiplications with `mul!`. When doing so, we will need to have cache variables to write into. This looks like:

````julia
Ayu = zeros(N,N)
uAx = zeros(N,N)
Du = zeros(N,N)
Ayv = zeros(N,N)
vAx = zeros(N,N)
Dv = zeros(N,N)
function gm3!(dr,r,p,t)
  a,α,ubar,β,D1,D2 = p
  u = @view r[:,:,1]
  v = @view r[:,:,2]
  du = @view dr[:,:,1]
  dv = @view dr[:,:,2]
  mul!(Ayu,Ay,u)
  mul!(uAx,u,Ax)
  mul!(Ayv,Ay,v)
  mul!(vAx,v,Ax)
  @. Du = D1*(Ayu + uAx)
  @. Dv = D2*(Ayv + vAx)
  @. du = Du + a*u*u./v + ubar - α*u
  @. dv = Dv + a*u*u - β*v
end
prob = ODEProblem(gm3!,r0,(0.0,0.1),p)
@benchmark solve(prob,Tsit5())
````


````
BenchmarkTools.Trial: 
  memory estimate:  29.75 MiB
  allocs estimate:  5318
  --------------
  minimum time:     33.112 ms (0.00% GC)
  median time:      33.238 ms (0.00% GC)
  mean time:        33.927 ms (2.02% GC)
  maximum time:     35.752 ms (6.71% GC)
  --------------
  samples:          148
  evals/sample:     1
````





But our temporary variables are global variables. We need to either declare the caches as `const` or localize them. We can localize them by adding them to the parameters, `p`. It's easier for the compiler to reason about local variables than global variables. ***Localizing variables helps to ensure type stability***.

````julia
p = (1.0,1.0,1.0,10.0,0.001,100.0,Ayu,uAx,Du,Ayv,vAx,Dv) # a,α,ubar,β,D1,D2
function gm4!(dr,r,p,t)
  a,α,ubar,β,D1,D2,Ayu,uAx,Du,Ayv,vAx,Dv = p
  u = @view r[:,:,1]
  v = @view r[:,:,2]
  du = @view dr[:,:,1]
  dv = @view dr[:,:,2]
  mul!(Ayu,Ay,u)
  mul!(uAx,u,Ax)
  mul!(Ayv,Ay,v)
  mul!(vAx,v,Ax)
  @. Du = D1*(Ayu + uAx)
  @. Dv = D2*(Ayv + vAx)
  @. du = Du + a*u*u./v + ubar - α*u
  @. dv = Dv + a*u*u - β*v
end
prob = ODEProblem(gm4!,r0,(0.0,0.1),p)
@benchmark solve(prob,Tsit5())
````


````
BenchmarkTools.Trial: 
  memory estimate:  29.66 MiB
  allocs estimate:  1054
  --------------
  minimum time:     26.598 ms (0.00% GC)
  median time:      26.729 ms (0.00% GC)
  mean time:        27.419 ms (2.47% GC)
  maximum time:     29.288 ms (8.08% GC)
  --------------
  samples:          183
  evals/sample:     1
````





We could then use the BLAS `gemmv` to optimize the matrix multiplications some more, but instead let's devectorize the stencil.

````julia
p = (1.0,1.0,1.0,10.0,0.001,100.0,N)
function fast_gm!(du,u,p,t)
  a,α,ubar,β,D1,D2,N = p

  @inbounds for j in 2:N-1, i in 2:N-1
    du[i,j,1] = D1*(u[i-1,j,1] + u[i+1,j,1] + u[i,j+1,1] + u[i,j-1,1] - 4u[i,j,1]) +
              a*u[i,j,1]^2/u[i,j,2] + ubar - α*u[i,j,1]
  end

  @inbounds for j in 2:N-1, i in 2:N-1
    du[i,j,2] = D2*(u[i-1,j,2] + u[i+1,j,2] + u[i,j+1,2] + u[i,j-1,2] - 4u[i,j,2]) +
            a*u[i,j,1]^2 - β*u[i,j,2]
  end

  @inbounds for j in 2:N-1
    i = 1
    du[1,j,1] = D1*(2u[i+1,j,1] + u[i,j+1,1] + u[i,j-1,1] - 4u[i,j,1]) +
            a*u[i,j,1]^2/u[i,j,2] + ubar - α*u[i,j,1]
  end
  @inbounds for j in 2:N-1
    i = 1
    du[1,j,2] = D2*(2u[i+1,j,2] + u[i,j+1,2] + u[i,j-1,2] - 4u[i,j,2]) +
            a*u[i,j,1]^2 - β*u[i,j,2]
  end
  @inbounds for j in 2:N-1
    i = N
    du[end,j,1] = D1*(2u[i-1,j,1] + u[i,j+1,1] + u[i,j-1,1] - 4u[i,j,1]) +
           a*u[i,j,1]^2/u[i,j,2] + ubar - α*u[i,j,1]
  end
  @inbounds for j in 2:N-1
    i = N
    du[end,j,2] = D2*(2u[i-1,j,2] + u[i,j+1,2] + u[i,j-1,2] - 4u[i,j,2]) +
           a*u[i,j,1]^2 - β*u[i,j,2]
  end

  @inbounds for i in 2:N-1
    j = 1
    du[i,1,1] = D1*(u[i-1,j,1] + u[i+1,j,1] + 2u[i,j+1,1] - 4u[i,j,1]) +
              a*u[i,j,1]^2/u[i,j,2] + ubar - α*u[i,j,1]
  end
  @inbounds for i in 2:N-1
    j = 1
    du[i,1,2] = D2*(u[i-1,j,2] + u[i+1,j,2] + 2u[i,j+1,2] - 4u[i,j,2]) +
              a*u[i,j,1]^2 - β*u[i,j,2]
  end
  @inbounds for i in 2:N-1
    j = N
    du[i,end,1] = D1*(u[i-1,j,1] + u[i+1,j,1] + 2u[i,j-1,1] - 4u[i,j,1]) +
             a*u[i,j,1]^2/u[i,j,2] + ubar - α*u[i,j,1]
  end
  @inbounds for i in 2:N-1
    j = N
    du[i,end,2] = D2*(u[i-1,j,2] + u[i+1,j,2] + 2u[i,j-1,2] - 4u[i,j,2]) +
             a*u[i,j,1]^2 - β*u[i,j,2]
  end

  @inbounds begin
    i = 1; j = 1
    du[1,1,1] = D1*(2u[i+1,j,1] + 2u[i,j+1,1] - 4u[i,j,1]) +
              a*u[i,j,1]^2/u[i,j,2] + ubar - α*u[i,j,1]
    du[1,1,2] = D2*(2u[i+1,j,2] + 2u[i,j+1,2] - 4u[i,j,2]) +
              a*u[i,j,1]^2 - β*u[i,j,2]

    i = 1; j = N
    du[1,N,1] = D1*(2u[i+1,j,1] + 2u[i,j-1,1] - 4u[i,j,1]) +
             a*u[i,j,1]^2/u[i,j,2] + ubar - α*u[i,j,1]
    du[1,N,2] = D2*(2u[i+1,j,2] + 2u[i,j-1,2] - 4u[i,j,2]) +
             a*u[i,j,1]^2 - β*u[i,j,2]

    i = N; j = 1
    du[N,1,1] = D1*(2u[i-1,j,1] + 2u[i,j+1,1] - 4u[i,j,1]) +
             a*u[i,j,1]^2/u[i,j,2] + ubar - α*u[i,j,1]
    du[N,1,2] = D2*(2u[i-1,j,2] + 2u[i,j+1,2] - 4u[i,j,2]) +
             a*u[i,j,1]^2 - β*u[i,j,2]

    i = N; j = N
    du[end,end,1] = D1*(2u[i-1,j,1] + 2u[i,j-1,1] - 4u[i,j,1]) +
             a*u[i,j,1]^2/u[i,j,2] + ubar - α*u[i,j,1]
    du[end,end,2] = D2*(2u[i-1,j,2] + 2u[i,j-1,2] - 4u[i,j,2]) +
             a*u[i,j,1]^2 - β*u[i,j,2]
   end
end
prob = ODEProblem(fast_gm!,r0,(0.0,0.1),p)
@benchmark solve(prob,Tsit5())
````


````
BenchmarkTools.Trial: 
  memory estimate:  29.62 MiB
  allocs estimate:  467
  --------------
  minimum time:     8.083 ms (0.00% GC)
  median time:      8.183 ms (0.00% GC)
  mean time:        8.841 ms (7.42% GC)
  maximum time:     10.712 ms (22.16% GC)
  --------------
  samples:          566
  evals/sample:     1
````





Lastly, we can do other things like multithread the main loops, but these optimizations get the last 2x-3x out. The main optimizations which apply everywhere are the ones we just performed (though the last one only works if your matrix is a stencil. This is known as a matrix-free implementation of the PDE discretization).

This gets us to about 8x faster than our original MATLAB/SciPy/R vectorized style code!

The last thing to do is then ***optimize our algorithm choice***. We have been using `Tsit5()` as our test algorithm, but in reality this problem is a stiff PDE discretization and thus one recommendation is to use `CVODE_BDF()`. However, instead of using the default dense Jacobian, we should make use of the sparse Jacobian afforded by the problem. The Jacobian is the matrix $\frac{df_i}{dr_j}$, where $r$ is read by the linear index (i.e. down columns). But since the $u$ variables depend on the $v$, the band size here is large, and thus this will not do well with a Banded Jacobian solver. Instead, we utilize sparse Jacobian algorithms. `CVODE_BDF` allows us to use a sparse Newton-Krylov solver by setting `linear_solver = :GMRES` (see [the solver documentation](https://docs.juliadiffeq.org/dev/solvers/ode_solve/#ode_solve_sundials-1), and thus we can solve this problem efficiently. Let's see how this scales as we increase the integration time.

````julia
prob = ODEProblem(fast_gm!,r0,(0.0,10.0),p)
@benchmark solve(prob,Tsit5())
````


````
BenchmarkTools.Trial: 
  memory estimate:  2.76 GiB
  allocs estimate:  41651
  --------------
  minimum time:     1.224 s (37.03% GC)
  median time:      1.605 s (40.31% GC)
  mean time:        1.533 s (42.60% GC)
  maximum time:     1.700 s (50.93% GC)
  --------------
  samples:          4
  evals/sample:     1
````



````julia
using Sundials
@benchmark solve(prob,CVODE_BDF(linear_solver=:GMRES))
````


````
BenchmarkTools.Trial: 
  memory estimate:  117.05 MiB
  allocs estimate:  20087
  --------------
  minimum time:     793.838 ms (0.00% GC)
  median time:      797.491 ms (0.34% GC)
  mean time:        797.332 ms (0.28% GC)
  maximum time:     799.731 ms (0.30% GC)
  --------------
  samples:          7
  evals/sample:     1
````



````julia
prob = ODEProblem(fast_gm!,r0,(0.0,100.0),p)
# Will go out of memory if we don't turn off `save_everystep`!
@benchmark solve(prob,Tsit5(),save_everystep=false)
````


````
BenchmarkTools.Trial: 
  memory estimate:  2.90 MiB
  allocs estimate:  77
  --------------
  minimum time:     5.673 s (0.00% GC)
  median time:      5.673 s (0.00% GC)
  mean time:        5.673 s (0.00% GC)
  maximum time:     5.673 s (0.00% GC)
  --------------
  samples:          1
  evals/sample:     1
````



````julia
@benchmark solve(prob,CVODE_BDF(linear_solver=:GMRES))
````


````
BenchmarkTools.Trial: 
  memory estimate:  303.39 MiB
  allocs estimate:  53851
  --------------
  minimum time:     2.184 s (0.00% GC)
  median time:      2.215 s (1.42% GC)
  mean time:        2.243 s (2.56% GC)
  maximum time:     2.329 s (6.04% GC)
  --------------
  samples:          3
  evals/sample:     1
````





Now let's check the allocation growth.

````julia
@benchmark solve(prob,CVODE_BDF(linear_solver=:GMRES),save_everystep=false)
````


````
BenchmarkTools.Trial: 
  memory estimate:  3.47 MiB
  allocs estimate:  44669
  --------------
  minimum time:     2.157 s (0.00% GC)
  median time:      2.158 s (0.00% GC)
  mean time:        2.159 s (0.00% GC)
  maximum time:     2.163 s (0.00% GC)
  --------------
  samples:          3
  evals/sample:     1
````



````julia
prob = ODEProblem(fast_gm!,r0,(0.0,500.0),p)
@benchmark solve(prob,CVODE_BDF(linear_solver=:GMRES),save_everystep=false)
````


````
BenchmarkTools.Trial: 
  memory estimate:  5.74 MiB
  allocs estimate:  84041
  --------------
  minimum time:     4.087 s (0.00% GC)
  median time:      4.089 s (0.00% GC)
  mean time:        4.089 s (0.00% GC)
  maximum time:     4.091 s (0.00% GC)
  --------------
  samples:          2
  evals/sample:     1
````





Notice that we've elimated almost all allocations, allowing the code to grow without hitting garbage collection and slowing down.

Why is `CVODE_BDF` doing well? What's happening is that, because the problem is stiff, the number of steps required by the explicit Runge-Kutta method grows rapidly, whereas `CVODE_BDF` is taking large steps. Additionally, the `GMRES` linear solver form is quite an efficient way to solve the implicit system in this case. This is problem-dependent, and in many cases using a Krylov method effectively requires a preconditioner, so you need to play around with testing other algorithms and linear solvers to find out what works best with your problem.

## Conclusion

Julia gives you the tools to optimize the solver "all the way", but you need to make use of it. The main thing to avoid is temporary allocations. For small systems, this is effectively done via static arrays. For large systems, this is done via in-place operations and cache arrays. Either way, the resulting solution can be immensely sped up over vectorized formulations by using these principles.


## Appendix

 This tutorial is part of the DiffEqTutorials.jl repository, found at: <https://github.com/JuliaDiffEq/DiffEqTutorials.jl>

To locally run this tutorial, do the following commands:
```
using DiffEqTutorials
DiffEqTutorials.weave_file("introduction","03-optimizing_diffeq_code.jmd")
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
  JULIA_CUDA_MEMORY_LIMIT = 2147483648
  JULIA_PROJECT = @.
  JULIA_NUM_THREADS = 4

```

Package Information:

```
Status `/builds/JuliaGPU/DiffEqTutorials.jl/tutorials/introduction/Project.toml`
[6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf] BenchmarkTools 0.5.0
[0c46a032-eb83-5123-abaf-570d42b7fbaa] DifferentialEquations 6.15.0
[65888b18-ceab-5e60-b2b9-181511a3b968] ParameterizedFunctions 5.3.0
[91a5bcdd-55d7-5caf-9e0b-520d859cae80] Plots 1.5.2
[90137ffa-7385-5640-81b9-e52037218182] StaticArrays 0.12.3
[c3572dad-4567-51f8-b174-8c6c989267f4] Sundials 4.2.5
[37e2e46d-f89d-539d-b4ee-838fcccc9c8e] LinearAlgebra
```
