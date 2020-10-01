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

```julia

function lorenz(u,p,t)
 dx = 10.0*(u[2]-u[1])
 dy = u[1]*(28.0-u[3]) - u[2]
 dz = u[1]*u[2] - (8/3)*u[3]
 [dx,dy,dz]
end
```


```
lorenz (generic function with 1 method)
```





Here, `lorenz` returns an object, `[dx,dy,dz]`, which is created within the body of `lorenz`.

This is a common code pattern from high-level languages like MATLAB, SciPy, or R's deSolve. However, the issue with this form is that it allocates a vector, `[dx,dy,dz]`, at each step. Let's benchmark the solution process with this choice of function:

```julia

using DifferentialEquations, BenchmarkTools
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz,u0,tspan)
@benchmark solve(prob,Tsit5())
```


```
BenchmarkTools.Trial: 
  memory estimate:  10.83 MiB
  allocs estimate:  101434
  --------------
  minimum time:     3.349 ms (0.00% GC)
  median time:      3.439 ms (0.00% GC)
  mean time:        4.139 ms (16.55% GC)
  maximum time:     8.662 ms (53.54% GC)
  --------------
  samples:          1207
  evals/sample:     1
```





The BenchmarkTools package's `@benchmark` runs the code multiple times to get an accurate measurement. The minimum time is the time it takes when your OS and other background processes aren't getting in the way. Notice that in this case it takes about 5ms to solve and allocates around 11.11 MiB. However, if we were to use this inside of a real user code we'd see a lot of time spent doing garbage collection (GC) to clean up all of the arrays we made. Even if we turn off saving we have these allocations.

```julia

@benchmark solve(prob,Tsit5(),save_everystep=false)
```


```
BenchmarkTools.Trial: 
  memory estimate:  9.47 MiB
  allocs estimate:  88652
  --------------
  minimum time:     2.911 ms (0.00% GC)
  median time:      2.919 ms (0.00% GC)
  mean time:        3.475 ms (14.76% GC)
  maximum time:     7.222 ms (56.45% GC)
  --------------
  samples:          1438
  evals/sample:     1
```





The problem of course is that arrays are created every time our derivative function is called. This function is called multiple times per step and is thus the main source of memory usage. To fix this, we can use the in-place form to ***make our code non-allocating***:

```julia

function lorenz!(du,u,p,t)
 du[1] = 10.0*(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end
```


```
lorenz! (generic function with 1 method)
```





Here, instead of creating an array each time, we utilized the cache array `du`. When the inplace form is used, DifferentialEquations.jl takes a different internal route that minimizes the internal allocations as well. When we benchmark this function, we will see quite a difference.

```julia

u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz!,u0,tspan)
@benchmark solve(prob,Tsit5())
```


```
BenchmarkTools.Trial: 
  memory estimate:  1.38 MiB
  allocs estimate:  12996
  --------------
  minimum time:     780.183 μs (0.00% GC)
  median time:      786.882 μs (0.00% GC)
  mean time:        876.801 μs (9.55% GC)
  maximum time:     5.902 ms (81.34% GC)
  --------------
  samples:          5690
  evals/sample:     1
```



```julia

@benchmark solve(prob,Tsit5(),save_everystep=false)
```


```
BenchmarkTools.Trial: 
  memory estimate:  5.22 KiB
  allocs estimate:  54
  --------------
  minimum time:     345.404 μs (0.00% GC)
  median time:      348.122 μs (0.00% GC)
  mean time:        348.437 μs (0.00% GC)
  maximum time:     374.122 μs (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1
```





There is a 4x time difference just from that change! Notice there are still some allocations and this is due to the construction of the integration cache. But this doesn't scale with the problem size:

```julia

tspan = (0.0,500.0) # 5x longer than before
prob = ODEProblem(lorenz!,u0,tspan)
@benchmark solve(prob,Tsit5(),save_everystep=false)
```


```
BenchmarkTools.Trial: 
  memory estimate:  5.22 KiB
  allocs estimate:  54
  --------------
  minimum time:     1.734 ms (0.00% GC)
  median time:      1.741 ms (0.00% GC)
  mean time:        1.743 ms (0.00% GC)
  maximum time:     1.896 ms (0.00% GC)
  --------------
  samples:          2867
  evals/sample:     1
```





since that's all just setup allocations.

#### But if the system is small we can optimize even more.

Allocations are only expensive if they are "heap allocations". For a more in-depth definition of heap allocations, [there are a lot of sources online](http://net-informations.com/faq/net/stack-heap.htm). But a good working definition is that heap allocations are variable-sized slabs of memory which have to be pointed to, and this pointer indirection costs time. Additionally, the heap has to be managed and the garbage controllers has to actively keep track of what's on the heap.

However, there's an alternative to heap allocations, known as stack allocations. The stack is statically-sized (known at compile time) and thus its accesses are quick. Additionally, the exact block of memory is known in advance by the compiler, and thus re-using the memory is cheap. This means that allocating on the stack has essentially no cost!

Arrays have to be heap allocated because their size (and thus the amount of memory they take up) is determined at runtime. But there are structures in Julia which are stack-allocated. `struct`s for example are stack-allocated "value-type"s. `Tuple`s are a stack-allocated collection. The most useful data structure for DiffEq though is the `StaticArray` from the package [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl). These arrays have their length determined at compile-time. They are created using macros attached to normal array expressions, for example:

```julia

using StaticArrays
A = @SVector [2.0,3.0,5.0]
```


```
3-element StaticArrays.SArray{Tuple{3},Float64,1,3} with indices SOneTo(3):
 2.0
 3.0
 5.0
```





Notice that the `3` after `SVector` gives the size of the `SVector`. It cannot be changed. Additionally, `SVector`s are immutable, so we have to create a new `SVector` to change values. But remember, we don't have to worry about allocations because this data structure is stack-allocated. `SArray`s have a lot of extra optimizations as well: they have fast matrix multiplication, fast QR factorizations, etc. which directly make use of the information about the size of the array. Thus, when possible they should be used.

Unfortunately static arrays can only be used for sufficiently small arrays. After a certain size, they are forced to heap allocate after some instructions and their compile time balloons. Thus static arrays shouldn't be used if your system has more than 100 variables. Additionally, only the native Julia algorithms can fully utilize static arrays.

Let's ***optimize `lorenz` using static arrays***. Note that in this case, we want to use the out-of-place allocating form, but this time we want to output a static array:

```julia

function lorenz_static(u,p,t)
 dx = 10.0*(u[2]-u[1])
 dy = u[1]*(28.0-u[3]) - u[2]
 dz = u[1]*u[2] - (8/3)*u[3]
 @SVector [dx,dy,dz]
end
```


```
lorenz_static (generic function with 1 method)
```





To make the solver internally use static arrays, we simply give it a static array as the initial condition:

```julia

u0 = @SVector [1.0,0.0,0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz_static,u0,tspan)
@benchmark solve(prob,Tsit5())
```


```
BenchmarkTools.Trial: 
  memory estimate:  466.28 KiB
  allocs estimate:  2595
  --------------
  minimum time:     338.227 μs (0.00% GC)
  median time:      342.915 μs (0.00% GC)
  mean time:        368.412 μs (4.75% GC)
  maximum time:     3.288 ms (86.15% GC)
  --------------
  samples:          10000
  evals/sample:     1
```



```julia

@benchmark solve(prob,Tsit5(),save_everystep=false)
```


```
BenchmarkTools.Trial: 
  memory estimate:  3.33 KiB
  allocs estimate:  28
  --------------
  minimum time:     236.345 μs (0.00% GC)
  median time:      239.987 μs (0.00% GC)
  mean time:        240.294 μs (0.00% GC)
  maximum time:     253.074 μs (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1
```





And that's pretty much all there is to it. With static arrays you don't have to worry about allocating, so use operations like `*` and don't worry about fusing operations (discussed in the next section). Do "the vectorized code" of R/MATLAB/Python and your code in this case will be fast, or directly use the numbers/values.

#### Exercise 1

Implement the out-of-place array, in-place array, and out-of-place static array forms for the [Henon-Heiles System](https://en.wikipedia.org/wiki/H%C3%A9non%E2%80%93Heiles_system) and time the results.

## Optimizing Large Systems

### Interlude: Managing Allocations with Broadcast Fusion

When your system is sufficiently large, or you have to make use of a non-native Julia algorithm, you have to make use of `Array`s. In order to use arrays in the most efficient manner, you need to be careful about temporary allocations. Vectorized calculations naturally have plenty of temporary array allocations. This is because a vectorized calculation outputs a vector. Thus:

```julia

A = rand(1000,1000); B = rand(1000,1000); C = rand(1000,1000)
test(A,B,C) = A + B + C
@benchmark test(A,B,C)
```


```
BenchmarkTools.Trial: 
  memory estimate:  7.63 MiB
  allocs estimate:  2
  --------------
  minimum time:     1.132 ms (0.00% GC)
  median time:      1.156 ms (0.00% GC)
  mean time:        1.305 ms (11.25% GC)
  maximum time:     2.961 ms (59.89% GC)
  --------------
  samples:          3814
  evals/sample:     1
```




That expression `A + B + C` creates 2 arrays. It first creates one for the output of `A + B`, then uses that result array to `+ C` to get the final result. 2 arrays! We don't want that! The first thing to do to fix this is to use broadcast fusion. [Broadcast fusion](https://julialang.org/blog/2017/01/moredots) puts expressions together. For example, instead of doing the `+` operations separately, if we were to add them all at the same time, then we would only have a single array that's created. For example:

```julia

test2(A,B,C) = map((a,b,c)->a+b+c,A,B,C)
@benchmark test2(A,B,C)
```


```
BenchmarkTools.Trial: 
  memory estimate:  7.63 MiB
  allocs estimate:  8
  --------------
  minimum time:     2.813 ms (0.00% GC)
  median time:      2.837 ms (0.00% GC)
  mean time:        2.973 ms (4.49% GC)
  maximum time:     4.549 ms (37.07% GC)
  --------------
  samples:          1680
  evals/sample:     1
```





Puts the whole expression into a single function call, and thus only one array is required to store output. This is the same as writing the loop:

```julia

function test3(A,B,C)
    D = similar(A)
    @inbounds for i in eachindex(A)
        D[i] = A[i] + B[i] + C[i]
    end
    D
end
@benchmark test3(A,B,C)
```


```
BenchmarkTools.Trial: 
  memory estimate:  7.63 MiB
  allocs estimate:  2
  --------------
  minimum time:     1.124 ms (0.00% GC)
  median time:      1.148 ms (0.00% GC)
  mean time:        1.297 ms (11.34% GC)
  maximum time:     3.000 ms (59.89% GC)
  --------------
  samples:          3837
  evals/sample:     1
```





However, Julia's broadcast is syntactic sugar for this. If multiple expressions have a `.`, then it will put those vectorized operations together. Thus:

```julia

test4(A,B,C) = A .+ B .+ C
@benchmark test4(A,B,C)
```


```
BenchmarkTools.Trial: 
  memory estimate:  7.63 MiB
  allocs estimate:  2
  --------------
  minimum time:     1.127 ms (0.00% GC)
  median time:      1.153 ms (0.00% GC)
  mean time:        1.302 ms (11.31% GC)
  maximum time:     2.967 ms (60.42% GC)
  --------------
  samples:          3821
  evals/sample:     1
```





is a version with only 1 array created (the output). Note that `.`s can be used with function calls as well:

```julia

sin.(A) .+ sin.(B)
```


```
1000×1000 Array{Float64,2}:
 1.14822   0.653106  0.400998  0.161702  …  0.931953  1.06832   0.900178
 1.29987   1.60977   0.59406   0.802195     0.696662  0.720217  0.151015
 0.724521  0.291946  1.15396   0.741696     1.11365   1.16862   0.713377
 0.261441  1.08603   0.838485  1.09141      0.542798  0.848069  0.780998
 1.28902   1.06844   0.896085  0.353321     1.49491   1.26847   0.775937
 0.931985  0.597447  1.54438   0.346974  …  0.571503  1.59248   0.920396
 0.64776   1.07941   1.22177   1.19306      1.23698   0.921395  0.91717
 1.36254   0.920541  1.1556    1.10763      0.308959  0.457159  0.971354
 1.56507   0.851284  0.786668  0.904077     0.638127  1.17275   1.33017
 0.802636  0.749744  0.780183  0.636848     0.750101  1.15993   0.953158
 ⋮                                       ⋱                      
 0.733594  1.36478   0.340724  0.340201     0.843932  0.680601  1.04149
 0.994858  0.887994  1.27176   1.42426      0.566199  0.829179  1.28069
 1.17642   1.14502   0.785152  1.24175      0.969501  1.47033   0.873622
 1.31984   0.857059  1.13175   1.38851      0.851331  0.486896  0.42976
 0.929718  0.832845  0.848255  0.978073  …  1.20457   1.1796    1.42254
 0.640584  0.627328  1.11305   1.19686      1.01449   0.659512  1.1893
 1.083     1.00496   1.17      0.190102     0.832462  1.18698   0.342249
 0.383132  1.20826   0.836149  0.807005     1.47502   1.10737   0.358481
 1.47624   1.11963   1.00009   1.23959      1.60455   0.426289  1.00275
```





Also, the `@.` macro applys a dot to every operator:

```julia

test5(A,B,C) = @. A + B + C #only one array allocated
@benchmark test5(A,B,C)
```


```
BenchmarkTools.Trial: 
  memory estimate:  7.63 MiB
  allocs estimate:  2
  --------------
  minimum time:     1.131 ms (0.00% GC)
  median time:      1.156 ms (0.00% GC)
  mean time:        1.305 ms (11.30% GC)
  maximum time:     2.983 ms (60.39% GC)
  --------------
  samples:          3813
  evals/sample:     1
```





Using these tools we can get rid of our intermediate array allocations for many vectorized function calls. But we are still allocating the output array. To get rid of that allocation, we can instead use mutation. Mutating broadcast is done via `.=`. For example, if we pre-allocate the output:

```julia
D = zeros(1000,1000);
```





Then we can keep re-using this cache for subsequent calculations. The mutating broadcasting form is:

```julia

test6!(D,A,B,C) = D .= A .+ B .+ C #only one array allocated
@benchmark test6!(D,A,B,C)
```


```
BenchmarkTools.Trial: 
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     1.060 ms (0.00% GC)
  median time:      1.067 ms (0.00% GC)
  mean time:        1.069 ms (0.00% GC)
  maximum time:     1.319 ms (0.00% GC)
  --------------
  samples:          4649
  evals/sample:     1
```





If we use `@.` before the `=`, then it will turn it into `.=`:

```julia

test7!(D,A,B,C) = @. D = A + B + C #only one array allocated
@benchmark test7!(D,A,B,C)
```


```
BenchmarkTools.Trial: 
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     1.062 ms (0.00% GC)
  median time:      1.069 ms (0.00% GC)
  mean time:        1.071 ms (0.00% GC)
  maximum time:     1.328 ms (0.00% GC)
  --------------
  samples:          4638
  evals/sample:     1
```





Notice that in this case, there is no "output", and instead the values inside of `D` are what are changed (like with the DiffEq inplace function). Many Julia functions have a mutating form which is denoted with a `!`. For example, the mutating form of the `map` is `map!`:

```julia

test8!(D,A,B,C) = map!((a,b,c)->a+b+c,D,A,B,C)
@benchmark test8!(D,A,B,C)
```


```
BenchmarkTools.Trial: 
  memory estimate:  32 bytes
  allocs estimate:  1
  --------------
  minimum time:     2.338 ms (0.00% GC)
  median time:      2.348 ms (0.00% GC)
  mean time:        2.353 ms (0.00% GC)
  maximum time:     3.030 ms (0.00% GC)
  --------------
  samples:          2121
  evals/sample:     1
```





Some operations require using an alternate mutating form in order to be fast. For example, matrix multiplication via `*` allocates a temporary:

```julia

@benchmark A*B
```


```
BenchmarkTools.Trial: 
  memory estimate:  7.63 MiB
  allocs estimate:  2
  --------------
  minimum time:     6.391 ms (0.00% GC)
  median time:      6.505 ms (0.00% GC)
  mean time:        6.669 ms (2.26% GC)
  maximum time:     8.716 ms (20.84% GC)
  --------------
  samples:          749
  evals/sample:     1
```





Instead, we can use the mutating form `mul!` into a cache array to avoid allocating the output:

```julia

using LinearAlgebra
@benchmark mul!(D,A,B) # same as D = A * B
```


```
BenchmarkTools.Trial: 
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     5.961 ms (0.00% GC)
  median time:      6.011 ms (0.00% GC)
  mean time:        6.015 ms (0.00% GC)
  maximum time:     6.647 ms (0.00% GC)
  --------------
  samples:          831
  evals/sample:     1
```





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

```julia

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
```


```
ODEProblem with uType Array{Float64,3} and tType Float64. In-place: true
timespan: (0.0, 0.1)
u0: [11.042756953924554 11.072858222058855 … 11.049539261990668 11.07123915
5670655; 11.026783271064723 11.00785249711271 … 11.061663284438994 11.02660
7567528439; … ; 11.055715295076121 11.04035895317539 … 11.079236661294487 1
1.055839497187407; 11.096194292032463 11.05546331683075 … 11.08225337132605
8 11.0034864004305]

[12.100000000000001 12.100000000000001 … 12.100000000000001 12.100000000000
001; 12.100000000000001 12.100000000000001 … 12.100000000000001 12.10000000
0000001; … ; 12.100000000000001 12.100000000000001 … 12.100000000000001 12.
100000000000001; 12.100000000000001 12.100000000000001 … 12.100000000000001
 12.100000000000001]
```





In this version we have encoded our initial condition to be a 3-dimensional array, with `u[:,:,1]` being the `A` part and `u[:,:,2]` being the `B` part.

```julia

@benchmark solve(prob,Tsit5())
```


```
BenchmarkTools.Trial: 
  memory estimate:  186.88 MiB
  allocs estimate:  8551
  --------------
  minimum time:     49.219 ms (4.41% GC)
  median time:      51.664 ms (8.32% GC)
  mean time:        51.201 ms (7.60% GC)
  maximum time:     51.823 ms (8.33% GC)
  --------------
  samples:          98
  evals/sample:     1
```





While this version isn't very efficient,

#### We recommend writing the "high-level" code first, and iteratively optimizing it!

The first thing that we can do is get rid of the slicing allocations. The operation `r[:,:,1]` creates a temporary array instead of a "view", i.e. a pointer to the already existing memory. To make it a view, add `@view`. Note that we have to be careful with views because they point to the same memory, and thus changing a view changes the original values:

```julia

A = rand(4)
@show A
B = @view A[1:3]
B[2] = 2
@show A
```


```
A = [0.547441153831203, 0.4889157576611225, 0.06710977826513576, 0.14941164
720627875]
A = [0.547441153831203, 2.0, 0.06710977826513576, 0.14941164720627875]
4-element Array{Float64,1}:
 0.547441153831203
 2.0
 0.06710977826513576
 0.14941164720627875
```





Notice that changing `B` changed `A`. This is something to be careful of, but at the same time we want to use this since we want to modify the output `dr`. Additionally, the last statement is a purely element-wise operation, and thus we can make use of broadcast fusion there. Let's rewrite `basic_version!` to ***avoid slicing allocations*** and to ***use broadcast fusion***:

```julia

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
```


```
BenchmarkTools.Trial: 
  memory estimate:  119.55 MiB
  allocs estimate:  7081
  --------------
  minimum time:     39.356 ms (5.58% GC)
  median time:      39.533 ms (5.57% GC)
  mean time:        40.047 ms (6.73% GC)
  maximum time:     41.897 ms (10.35% GC)
  --------------
  samples:          125
  evals/sample:     1
```





Now, most of the allocations are taking place in `Du = D1*(Ay*u + u*Ax)` since those operations are vectorized and not mutating. We should instead replace the matrix multiplications with `mul!`. When doing so, we will need to have cache variables to write into. This looks like:

```julia

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
```


```
BenchmarkTools.Trial: 
  memory estimate:  29.75 MiB
  allocs estimate:  5317
  --------------
  minimum time:     33.887 ms (0.00% GC)
  median time:      34.208 ms (0.00% GC)
  mean time:        34.773 ms (1.82% GC)
  maximum time:     36.550 ms (5.81% GC)
  --------------
  samples:          144
  evals/sample:     1
```





But our temporary variables are global variables. We need to either declare the caches as `const` or localize them. We can localize them by adding them to the parameters, `p`. It's easier for the compiler to reason about local variables than global variables. ***Localizing variables helps to ensure type stability***.

```julia

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
```


```
BenchmarkTools.Trial: 
  memory estimate:  29.66 MiB
  allocs estimate:  1053
  --------------
  minimum time:     27.394 ms (0.00% GC)
  median time:      27.715 ms (0.00% GC)
  mean time:        28.289 ms (2.23% GC)
  maximum time:     30.058 ms (7.15% GC)
  --------------
  samples:          177
  evals/sample:     1
```





We could then use the BLAS `gemmv` to optimize the matrix multiplications some more, but instead let's devectorize the stencil.

```julia

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
```


```
BenchmarkTools.Trial: 
  memory estimate:  29.62 MiB
  allocs estimate:  466
  --------------
  minimum time:     8.368 ms (0.00% GC)
  median time:      8.490 ms (0.00% GC)
  mean time:        9.109 ms (6.77% GC)
  maximum time:     10.770 ms (19.44% GC)
  --------------
  samples:          549
  evals/sample:     1
```





Lastly, we can do other things like multithread the main loops, but these optimizations get the last 2x-3x out. The main optimizations which apply everywhere are the ones we just performed (though the last one only works if your matrix is a stencil. This is known as a matrix-free implementation of the PDE discretization).

This gets us to about 8x faster than our original MATLAB/SciPy/R vectorized style code!

The last thing to do is then ***optimize our algorithm choice***. We have been using `Tsit5()` as our test algorithm, but in reality this problem is a stiff PDE discretization and thus one recommendation is to use `CVODE_BDF()`. However, instead of using the default dense Jacobian, we should make use of the sparse Jacobian afforded by the problem. The Jacobian is the matrix $\frac{df_i}{dr_j}$, where $r$ is read by the linear index (i.e. down columns). But since the $u$ variables depend on the $v$, the band size here is large, and thus this will not do well with a Banded Jacobian solver. Instead, we utilize sparse Jacobian algorithms. `CVODE_BDF` allows us to use a sparse Newton-Krylov solver by setting `linear_solver = :GMRES` (see [the solver documentation](https://docs.sciml.ai/dev/solvers/ode_solve/#ode_solve_sundials-1), and thus we can solve this problem efficiently. Let's see how this scales as we increase the integration time.

```julia

prob = ODEProblem(fast_gm!,r0,(0.0,10.0),p)
@benchmark solve(prob,Tsit5())
```


```
BenchmarkTools.Trial: 
  memory estimate:  2.76 GiB
  allocs estimate:  41632
  --------------
  minimum time:     1.342 s (41.01% GC)
  median time:      1.550 s (39.28% GC)
  mean time:        1.571 s (41.38% GC)
  maximum time:     1.840 s (46.42% GC)
  --------------
  samples:          4
  evals/sample:     1
```



```julia

using Sundials
@benchmark solve(prob,CVODE_BDF(linear_solver=:GMRES))
```


```
BenchmarkTools.Trial: 
  memory estimate:  51.87 MiB
  allocs estimate:  8371
  --------------
  minimum time:     335.117 ms (0.00% GC)
  median time:      337.471 ms (0.64% GC)
  mean time:        337.130 ms (0.37% GC)
  maximum time:     338.900 ms (0.68% GC)
  --------------
  samples:          15
  evals/sample:     1
```



```julia

prob = ODEProblem(fast_gm!,r0,(0.0,100.0),p)
# Will go out of memory if we don't turn off `save_everystep`!
@benchmark solve(prob,Tsit5(),save_everystep=false)
```


```
BenchmarkTools.Trial: 
  memory estimate:  2.90 MiB
  allocs estimate:  74
  --------------
  minimum time:     5.886 s (0.00% GC)
  median time:      5.886 s (0.00% GC)
  mean time:        5.886 s (0.00% GC)
  maximum time:     5.886 s (0.00% GC)
  --------------
  samples:          1
  evals/sample:     1
```



```julia

@benchmark solve(prob,CVODE_BDF(linear_solver=:GMRES))
```


```
BenchmarkTools.Trial: 
  memory estimate:  302.02 MiB
  allocs estimate:  53934
  --------------
  minimum time:     2.211 s (0.00% GC)
  median time:      2.253 s (1.36% GC)
  mean time:        2.272 s (2.36% GC)
  maximum time:     2.353 s (5.53% GC)
  --------------
  samples:          3
  evals/sample:     1
```





Now let's check the allocation growth.

```julia

@benchmark solve(prob,CVODE_BDF(linear_solver=:GMRES),save_everystep=false)
```


```
BenchmarkTools.Trial: 
  memory estimate:  3.48 MiB
  allocs estimate:  44791
  --------------
  minimum time:     2.185 s (0.00% GC)
  median time:      2.194 s (0.00% GC)
  mean time:        2.192 s (0.00% GC)
  maximum time:     2.198 s (0.00% GC)
  --------------
  samples:          3
  evals/sample:     1
```



```julia

prob = ODEProblem(fast_gm!,r0,(0.0,500.0),p)
@benchmark solve(prob,CVODE_BDF(linear_solver=:GMRES),save_everystep=false)
```


```
BenchmarkTools.Trial: 
  memory estimate:  4.93 MiB
  allocs estimate:  69991
  --------------
  minimum time:     3.421 s (0.00% GC)
  median time:      3.422 s (0.00% GC)
  mean time:        3.422 s (0.00% GC)
  maximum time:     3.423 s (0.00% GC)
  --------------
  samples:          2
  evals/sample:     1
```





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
[65888b18-ceab-5e60-b2b9-181511a3b968] ParameterizedFunctions 5.6.0
[91a5bcdd-55d7-5caf-9e0b-520d859cae80] Plots 1.6.7
[90137ffa-7385-5640-81b9-e52037218182] StaticArrays 0.12.4
[c3572dad-4567-51f8-b174-8c6c989267f4] Sundials 4.3.0
[37e2e46d-f89d-539d-b4ee-838fcccc9c8e] LinearAlgebra
```
