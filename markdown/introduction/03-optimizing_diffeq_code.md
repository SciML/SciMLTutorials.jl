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
  minimum time:     3.314 ms (0.00% GC)
  median time:      3.349 ms (0.00% GC)
  mean time:        4.227 ms (16.91% GC)
  maximum time:     9.396 ms (47.13% GC)
  --------------
  samples:          1184
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
  minimum time:     2.906 ms (0.00% GC)
  median time:      2.924 ms (0.00% GC)
  mean time:        3.786 ms (15.89% GC)
  maximum time:     8.883 ms (49.49% GC)
  --------------
  samples:          1321
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
  minimum time:     779.317 μs (0.00% GC)
  median time:      787.006 μs (0.00% GC)
  mean time:        890.525 μs (9.62% GC)
  maximum time:     6.121 ms (81.63% GC)
  --------------
  samples:          5602
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
  minimum time:     343.929 μs (0.00% GC)
  median time:      346.587 μs (0.00% GC)
  mean time:        346.911 μs (0.00% GC)
  maximum time:     381.913 μs (0.00% GC)
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
  minimum time:     1.724 ms (0.00% GC)
  median time:      1.730 ms (0.00% GC)
  mean time:        1.732 ms (0.00% GC)
  maximum time:     1.855 ms (0.00% GC)
  --------------
  samples:          2886
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
  minimum time:     341.344 μs (0.00% GC)
  median time:      346.586 μs (0.00% GC)
  mean time:        377.242 μs (5.40% GC)
  maximum time:     3.822 ms (87.09% GC)
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
  minimum time:     235.138 μs (0.00% GC)
  median time:      238.182 μs (0.00% GC)
  mean time:        238.641 μs (0.00% GC)
  maximum time:     268.689 μs (0.00% GC)
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
  minimum time:     1.129 ms (0.00% GC)
  median time:      1.153 ms (0.00% GC)
  mean time:        1.305 ms (11.52% GC)
  maximum time:     3.004 ms (60.22% GC)
  --------------
  samples:          3816
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
  minimum time:     2.785 ms (0.00% GC)
  median time:      2.820 ms (0.00% GC)
  mean time:        2.961 ms (4.71% GC)
  maximum time:     4.565 ms (37.82% GC)
  --------------
  samples:          1687
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
  minimum time:     1.118 ms (0.00% GC)
  median time:      1.142 ms (0.00% GC)
  mean time:        1.293 ms (11.61% GC)
  maximum time:     2.991 ms (61.21% GC)
  --------------
  samples:          3848
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
  minimum time:     1.124 ms (0.00% GC)
  median time:      1.146 ms (0.00% GC)
  mean time:        1.298 ms (11.60% GC)
  maximum time:     3.003 ms (60.76% GC)
  --------------
  samples:          3835
  evals/sample:     1
````





is a version with only 1 array created (the output). Note that `.`s can be used with function calls as well:

````julia

sin.(A) .+ sin.(B)
````


````
1000×1000 Array{Float64,2}:
 0.522869  0.762457  0.705335  1.28392    …  0.697115  0.586575  0.712355
 0.757357  1.45263   0.575005  0.839452      1.14611   0.926982  0.618261
 0.892628  0.961034  0.383037  1.13252       0.771699  0.890253  1.13169
 0.280007  1.24839   1.01319   0.894068      0.823251  1.28038   0.268521
 0.438326  0.681747  0.912342  0.628293      0.268904  0.84825   1.24767
 1.29953   1.3449    1.0407    0.69007    …  0.516017  0.84036   0.750896
 1.06155   1.11441   0.421359  1.16101       0.198726  0.779754  1.48158
 1.04738   0.696495  1.17574   1.58097       0.896606  0.803252  0.470437
 0.41159   1.33406   1.12273   0.480433      0.486254  1.1082    0.111666
 1.24669   1.15907   0.564737  0.26263       1.31374   1.42007   0.937338
 ⋮                                        ⋱                      
 1.31633   1.47947   0.804315  0.649385      0.467765  1.48512   0.842383
 0.648239  0.794633  1.05609   1.14878       1.14888   0.782322  0.94693
 0.616616  1.39635   1.34339   0.160978      0.782674  1.2463    0.844689
 0.987097  0.967679  1.0116    1.07016       0.565597  0.820391  1.39805
 0.786271  0.703329  1.11445   0.380049   …  0.611932  0.471022  0.747319
 0.858423  0.602896  0.870189  0.667962      1.50664   0.217334  1.17214
 0.228472  0.860427  0.259358  0.0224333     0.883119  1.12868   1.05434
 0.47911   0.303095  0.842932  0.803767      0.620023  0.765244  0.804916
 1.27327   1.01102   0.642231  1.07592       0.949659  0.811927  0.975322
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
  minimum time:     1.127 ms (0.00% GC)
  median time:      1.153 ms (0.00% GC)
  mean time:        1.305 ms (11.60% GC)
  maximum time:     3.010 ms (60.49% GC)
  --------------
  samples:          3814
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
  median time:      1.046 ms (0.00% GC)
  mean time:        1.048 ms (0.00% GC)
  maximum time:     1.297 ms (0.00% GC)
  --------------
  samples:          4746
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
  minimum time:     1.046 ms (0.00% GC)
  median time:      1.055 ms (0.00% GC)
  mean time:        1.056 ms (0.00% GC)
  maximum time:     1.309 ms (0.00% GC)
  --------------
  samples:          4710
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
  minimum time:     2.345 ms (0.00% GC)
  median time:      2.354 ms (0.00% GC)
  mean time:        2.358 ms (0.00% GC)
  maximum time:     2.538 ms (0.00% GC)
  --------------
  samples:          2117
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
  minimum time:     6.527 ms (0.00% GC)
  median time:      6.710 ms (0.00% GC)
  mean time:        6.827 ms (2.08% GC)
  maximum time:     8.602 ms (21.50% GC)
  --------------
  samples:          732
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
  minimum time:     6.059 ms (0.00% GC)
  median time:      6.102 ms (0.00% GC)
  mean time:        6.104 ms (0.00% GC)
  maximum time:     6.702 ms (0.00% GC)
  --------------
  samples:          818
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
u0: [11.074738640769288 11.067258154135137 … 11.081834930743078 11.09898091
6537101; 11.09248438940732 11.081536332572746 … 11.065826658606165 11.08996
2204167723; … ; 11.095783201556893 11.063543569268473 … 11.070578958831506 
11.057631627995505; 11.036811963628905 11.079533065539584 … 11.058553648429
877 11.029358884339073]

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
  minimum time:     48.874 ms (4.48% GC)
  median time:      51.295 ms (8.50% GC)
  mean time:        50.751 ms (7.59% GC)
  maximum time:     51.455 ms (8.60% GC)
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
B = @view A[1:3]
B[2] = 2
@show A
````


````
A = [0.5522111940313656, 0.8844994354841427, 0.8492692302056744, 0.19853348
991877362]
A = [0.5522111940313656, 2.0, 0.8492692302056744, 0.19853348991877362]
4-element Array{Float64,1}:
 0.5522111940313656
 2.0
 0.8492692302056744
 0.19853348991877362
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
  minimum time:     38.841 ms (5.42% GC)
  median time:      39.225 ms (5.68% GC)
  mean time:        39.644 ms (6.70% GC)
  maximum time:     41.672 ms (10.48% GC)
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
  minimum time:     33.729 ms (0.00% GC)
  median time:      33.885 ms (0.00% GC)
  mean time:        34.519 ms (1.90% GC)
  maximum time:     36.290 ms (6.21% GC)
  --------------
  samples:          145
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
  minimum time:     27.036 ms (0.00% GC)
  median time:      27.232 ms (0.00% GC)
  mean time:        27.889 ms (2.32% GC)
  maximum time:     29.694 ms (7.53% GC)
  --------------
  samples:          180
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
  minimum time:     8.334 ms (0.00% GC)
  median time:      8.459 ms (0.00% GC)
  mean time:        9.081 ms (6.83% GC)
  maximum time:     10.852 ms (19.87% GC)
  --------------
  samples:          551
  evals/sample:     1
````





Lastly, we can do other things like multithread the main loops, but these optimizations get the last 2x-3x out. The main optimizations which apply everywhere are the ones we just performed (though the last one only works if your matrix is a stencil. This is known as a matrix-free implementation of the PDE discretization).

This gets us to about 8x faster than our original MATLAB/SciPy/R vectorized style code!

The last thing to do is then ***optimize our algorithm choice***. We have been using `Tsit5()` as our test algorithm, but in reality this problem is a stiff PDE discretization and thus one recommendation is to use `CVODE_BDF()`. However, instead of using the default dense Jacobian, we should make use of the sparse Jacobian afforded by the problem. The Jacobian is the matrix $\frac{df_i}{dr_j}$, where $r$ is read by the linear index (i.e. down columns). But since the $u$ variables depend on the $v$, the band size here is large, and thus this will not do well with a Banded Jacobian solver. Instead, we utilize sparse Jacobian algorithms. `CVODE_BDF` allows us to use a sparse Newton-Krylov solver by setting `linear_solver = :GMRES` (see [the solver documentation](https://docs.sciml.ai/dev/solvers/ode_solve/#ode_solve_sundials-1), and thus we can solve this problem efficiently. Let's see how this scales as we increase the integration time.

````julia

prob = ODEProblem(fast_gm!,r0,(0.0,10.0),p)
@benchmark solve(prob,Tsit5())
````


````
BenchmarkTools.Trial: 
  memory estimate:  2.76 GiB
  allocs estimate:  41633
  --------------
  minimum time:     1.361 s (41.01% GC)
  median time:      1.513 s (40.92% GC)
  mean time:        1.502 s (40.83% GC)
  maximum time:     1.621 s (50.09% GC)
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
  memory estimate:  55.10 MiB
  allocs estimate:  8770
  --------------
  minimum time:     351.671 ms (0.00% GC)
  median time:      354.166 ms (0.63% GC)
  mean time:        354.143 ms (0.36% GC)
  maximum time:     356.730 ms (0.69% GC)
  --------------
  samples:          15
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
  minimum time:     5.971 s (0.00% GC)
  median time:      5.971 s (0.00% GC)
  mean time:        5.971 s (0.00% GC)
  maximum time:     5.971 s (0.00% GC)
  --------------
  samples:          1
  evals/sample:     1
````



````julia

@benchmark solve(prob,CVODE_BDF(linear_solver=:GMRES))
````


````
BenchmarkTools.Trial: 
  memory estimate:  240.88 MiB
  allocs estimate:  40993
  --------------
  minimum time:     1.690 s (0.27% GC)
  median time:      1.717 s (1.51% GC)
  mean time:        1.742 s (3.17% GC)
  maximum time:     1.820 s (7.43% GC)
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
  memory estimate:  2.84 MiB
  allocs estimate:  33701
  --------------
  minimum time:     1.659 s (0.00% GC)
  median time:      1.660 s (0.00% GC)
  mean time:        1.660 s (0.00% GC)
  maximum time:     1.662 s (0.00% GC)
  --------------
  samples:          4
  evals/sample:     1
````



````julia

prob = ODEProblem(fast_gm!,r0,(0.0,500.0),p)
@benchmark solve(prob,CVODE_BDF(linear_solver=:GMRES),save_everystep=false)
````


````
BenchmarkTools.Trial: 
  memory estimate:  3.48 MiB
  allocs estimate:  44963
  --------------
  minimum time:     2.210 s (0.00% GC)
  median time:      2.216 s (0.00% GC)
  mean time:        2.216 s (0.00% GC)
  maximum time:     2.222 s (0.00% GC)
  --------------
  samples:          3
  evals/sample:     1
````





Notice that we've elimated almost all allocations, allowing the code to grow without hitting garbage collection and slowing down.

Why is `CVODE_BDF` doing well? What's happening is that, because the problem is stiff, the number of steps required by the explicit Runge-Kutta method grows rapidly, whereas `CVODE_BDF` is taking large steps. Additionally, the `GMRES` linear solver form is quite an efficient way to solve the implicit system in this case. This is problem-dependent, and in many cases using a Krylov method effectively requires a preconditioner, so you need to play around with testing other algorithms and linear solvers to find out what works best with your problem.

## Conclusion

Julia gives you the tools to optimize the solver "all the way", but you need to make use of it. The main thing to avoid is temporary allocations. For small systems, this is effectively done via static arrays. For large systems, this is done via in-place operations and cache arrays. Either way, the resulting solution can be immensely sped up over vectorized formulations by using these principles.


## Appendix

 This tutorial is part of the SciMLTutorials.jl repository, found at: <https://github.com/SciML/SciMLTutorials.jl>.
 For more information on doing scientific machine learning (SciML) with open source software, check out <https://sciml.ai/>.

To locally run this tutorial, do the following commands:
```
using SciMLTutorials
SciMLTutorials.weave_file("introduction","03-optimizing_diffeq_code.jmd")
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
  JULIA_LOAD_PATH = /builds/JuliaGPU/DiffEqTutorials.jl:
  JULIA_DEPOT_PATH = /builds/JuliaGPU/DiffEqTutorials.jl/.julia
  JULIA_CUDA_MEMORY_LIMIT = 2147483648
  JULIA_NUM_THREADS = 8

```

Package Information:

```
Status `/builds/JuliaGPU/DiffEqTutorials.jl/tutorials/introduction/Project.toml`
[6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf] BenchmarkTools 0.5.0
[0c46a032-eb83-5123-abaf-570d42b7fbaa] DifferentialEquations 6.15.0
[65888b18-ceab-5e60-b2b9-181511a3b968] ParameterizedFunctions 5.4.0
[91a5bcdd-55d7-5caf-9e0b-520d859cae80] Plots 1.6.0
[90137ffa-7385-5640-81b9-e52037218182] StaticArrays 0.12.4
[c3572dad-4567-51f8-b174-8c6c989267f4] Sundials 4.2.6
[37e2e46d-f89d-539d-b4ee-838fcccc9c8e] LinearAlgebra
```
