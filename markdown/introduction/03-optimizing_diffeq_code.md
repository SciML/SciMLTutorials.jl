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
  allocs estimate:  101447
  --------------
  minimum time:     3.860 ms (0.00% GC)
  median time:      8.464 ms (0.00% GC)
  mean time:        8.830 ms (19.87% GC)
  maximum time:     50.326 ms (83.04% GC)
  --------------
  samples:          566
  evals/sample:     1
````





The BenchmarkTools package's `@benchmark` runs the code multiple times to get an accurate measurement. The minimum time is the time it takes when your OS and other background processes aren't getting in the way. Notice that in this case it takes about 5ms to solve and allocates around 11.11 MiB. However, if we were to use this inside of a real user code we'd see a lot of time spent doing garbage collection (GC) to clean up all of the arrays we made. Even if we turn off saving we have these allocations.

````julia
@benchmark solve(prob,Tsit5(),save_everystep=false)
````


````
BenchmarkTools.Trial: 
  memory estimate:  9.47 MiB
  allocs estimate:  88665
  --------------
  minimum time:     3.162 ms (0.00% GC)
  median time:      7.178 ms (0.00% GC)
  mean time:        7.653 ms (20.27% GC)
  maximum time:     49.797 ms (85.01% GC)
  --------------
  samples:          658
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
  allocs estimate:  13009
  --------------
  minimum time:     718.400 μs (0.00% GC)
  median time:      1.261 ms (0.00% GC)
  mean time:        1.396 ms (14.73% GC)
  maximum time:     43.266 ms (96.46% GC)
  --------------
  samples:          3575
  evals/sample:     1
````



````julia
@benchmark solve(prob,Tsit5(),save_everystep=false)
````


````
BenchmarkTools.Trial: 
  memory estimate:  6.28 KiB
  allocs estimate:  67
  --------------
  minimum time:     310.300 μs (0.00% GC)
  median time:      363.450 μs (0.00% GC)
  mean time:        378.548 μs (0.00% GC)
  maximum time:     688.700 μs (0.00% GC)
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
  memory estimate:  6.28 KiB
  allocs estimate:  67
  --------------
  minimum time:     1.549 ms (0.00% GC)
  median time:      1.828 ms (0.00% GC)
  mean time:        1.779 ms (0.00% GC)
  maximum time:     3.835 ms (0.00% GC)
  --------------
  samples:          2807
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
  memory estimate:  467.31 KiB
  allocs estimate:  2608
  --------------
  minimum time:     318.199 μs (0.00% GC)
  median time:      493.000 μs (0.00% GC)
  mean time:        537.263 μs (11.66% GC)
  maximum time:     36.188 ms (98.46% GC)
  --------------
  samples:          9284
  evals/sample:     1
````



````julia
@benchmark solve(prob,Tsit5(),save_everystep=false)
````


````
BenchmarkTools.Trial: 
  memory estimate:  4.36 KiB
  allocs estimate:  41
  --------------
  minimum time:     205.701 μs (0.00% GC)
  median time:      224.499 μs (0.00% GC)
  mean time:        227.612 μs (0.00% GC)
  maximum time:     510.401 μs (0.00% GC)
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
  minimum time:     3.551 ms (0.00% GC)
  median time:      4.203 ms (0.00% GC)
  mean time:        5.097 ms (17.47% GC)
  maximum time:     38.441 ms (88.83% GC)
  --------------
  samples:          979
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
  allocs estimate:  5
  --------------
  minimum time:     4.374 ms (0.00% GC)
  median time:      5.117 ms (0.00% GC)
  mean time:        5.923 ms (14.60% GC)
  maximum time:     37.805 ms (86.14% GC)
  --------------
  samples:          843
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
  minimum time:     3.500 ms (0.00% GC)
  median time:      4.170 ms (0.00% GC)
  mean time:        5.035 ms (17.39% GC)
  maximum time:     38.428 ms (88.59% GC)
  --------------
  samples:          991
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
  minimum time:     3.575 ms (0.00% GC)
  median time:      4.178 ms (0.00% GC)
  mean time:        5.064 ms (17.32% GC)
  maximum time:     37.581 ms (88.06% GC)
  --------------
  samples:          985
  evals/sample:     1
````





is a version with only 1 array created (the output). Note that `.`s can be used with function calls as well:

````julia
sin.(A) .+ sin.(B)
````


````
1000×1000 Array{Float64,2}:
 0.629199  1.10191   0.908789  0.601001  …  1.07944   1.54387   0.998294
 0.63178   1.27845   1.25595   0.545865     0.950411  1.10605   1.14804 
 0.74961   0.338372  0.776325  0.843301     0.621682  1.4814    0.387066
 1.12688   0.766303  1.33449   0.766908     1.02491   1.31279   1.35948 
 1.0273    1.28566   0.796759  1.06933      1.11355   0.941753  0.809312
 0.847441  0.703465  1.12598   0.664446  …  0.353486  0.113834  0.882723
 1.11164   0.909755  1.08992   0.886073     1.00769   1.02535   0.986124
 0.600466  1.22593   0.381924  1.1188       1.06073   0.928734  0.860063
 1.18456   0.930343  1.13676   1.67413      1.29347   0.376713  0.947774
 0.54174   1.39467   0.491944  0.438106     0.723928  1.194     1.09252 
 ⋮                                       ⋱                              
 1.02091   1.36771   0.921309  0.753442     0.991678  0.579141  1.00315 
 0.919283  0.839523  0.635499  1.11653      1.09411   0.921141  0.805901
 0.357252  0.897618  1.13678   0.834228     1.63192   1.56706   0.72399 
 1.28705   0.729184  1.32264   0.942499     1.25033   0.479953  1.08578 
 1.39789   0.817921  0.627723  1.21801   …  1.05692   0.851231  0.578276
 0.904645  0.833058  1.41274   0.696929     0.775396  1.47136   1.06926 
 0.984909  1.47736   1.08829   1.61864      0.511166  1.36079   1.17073 
 1.22154   1.25806   0.848001  1.21803      1.47232   1.08072   1.37376 
 0.946084  0.741297  0.43729   0.852772     1.22975   0.891545  0.84901
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
  minimum time:     3.606 ms (0.00% GC)
  median time:      4.221 ms (0.00% GC)
  mean time:        5.107 ms (17.55% GC)
  maximum time:     39.207 ms (88.49% GC)
  --------------
  samples:          977
  evals/sample:     1
````





Using these tools we can get rid of our intermediate array allocations for many vectorized function calls. But we are still allocating the output array. To get rid of that allocation, we can instead use mutation. Mutating broadcast is done via `.=`. For example, if we pre-allocate the output:

````julia
D = zeros(1000,1000);
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
  minimum time:     1.536 ms (0.00% GC)
  median time:      1.636 ms (0.00% GC)
  mean time:        1.655 ms (0.00% GC)
  maximum time:     2.217 ms (0.00% GC)
  --------------
  samples:          3001
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
  minimum time:     1.548 ms (0.00% GC)
  median time:      1.646 ms (0.00% GC)
  mean time:        1.664 ms (0.00% GC)
  maximum time:     2.256 ms (0.00% GC)
  --------------
  samples:          2984
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
  minimum time:     2.645 ms (0.00% GC)
  median time:      2.790 ms (0.00% GC)
  mean time:        2.845 ms (0.00% GC)
  maximum time:     3.788 ms (0.00% GC)
  --------------
  samples:          1752
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
  minimum time:     16.939 ms (0.00% GC)
  median time:      27.255 ms (0.00% GC)
  mean time:        27.735 ms (4.88% GC)
  maximum time:     88.387 ms (60.92% GC)
  --------------
  samples:          181
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
  minimum time:     14.639 ms (0.00% GC)
  median time:      24.124 ms (0.00% GC)
  mean time:        23.555 ms (0.00% GC)
  maximum time:     45.502 ms (0.00% GC)
  --------------
  samples:          212
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
u0: [11.088901276293319 11.023925325874055 … 11.046600137596931 11.06229787
4592385; 11.071035374548634 11.03407650194175 … 11.007487343004286 11.08108
2781677436; … ; 11.007445436693185 11.034901131743307 … 11.090654974956067 
11.006764746744533; 11.008229200769804 11.028737575584314 … 11.039895849449
023 11.09610574895775]

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
  allocs estimate:  8564
  --------------
  minimum time:     96.121 ms (0.00% GC)
  median time:      150.949 ms (23.11% GC)
  mean time:        141.973 ms (16.57% GC)
  maximum time:     207.772 ms (27.91% GC)
  --------------
  samples:          36
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
A = [0.8310220901369554, 0.3727296355218075, 0.0835910683892136, 0.56188607
68016369]
````



````julia
B = @view A[1:3]
B[2] = 2
@show A
````


````
A = [0.8310220901369554, 2.0, 0.0835910683892136, 0.5618860768016369]
4-element Array{Float64,1}:
 0.8310220901369554
 2.0               
 0.0835910683892136
 0.5618860768016369
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
  allocs estimate:  7094
  --------------
  minimum time:     76.812 ms (0.00% GC)
  median time:      105.176 ms (0.00% GC)
  mean time:        115.817 ms (13.81% GC)
  maximum time:     201.290 ms (32.51% GC)
  --------------
  samples:          44
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
  memory estimate:  29.76 MiB
  allocs estimate:  5330
  --------------
  minimum time:     54.567 ms (0.00% GC)
  median time:      63.950 ms (0.00% GC)
  mean time:        66.450 ms (4.17% GC)
  maximum time:     101.438 ms (27.05% GC)
  --------------
  samples:          76
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
  allocs estimate:  1066
  --------------
  minimum time:     46.456 ms (0.00% GC)
  median time:      55.021 ms (0.00% GC)
  mean time:        57.857 ms (4.80% GC)
  maximum time:     96.047 ms (34.60% GC)
  --------------
  samples:          87
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
  memory estimate:  29.63 MiB
  allocs estimate:  479
  --------------
  minimum time:     11.207 ms (0.00% GC)
  median time:      12.684 ms (0.00% GC)
  mean time:        16.075 ms (14.75% GC)
  maximum time:     70.401 ms (67.62% GC)
  --------------
  samples:          311
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
  allocs estimate:  41645
  --------------
  minimum time:     2.323 s (16.87% GC)
  median time:      2.644 s (25.79% GC)
  mean time:        2.644 s (25.79% GC)
  maximum time:     2.964 s (32.78% GC)
  --------------
  samples:          2
  evals/sample:     1
````



````julia
using Sundials
@benchmark solve(prob,CVODE_BDF(linear_solver=:GMRES))
````


````
BenchmarkTools.Trial: 
  memory estimate:  120.75 MiB
  allocs estimate:  20296
  --------------
  minimum time:     634.511 ms (0.00% GC)
  median time:      708.681 ms (0.00% GC)
  mean time:        699.203 ms (1.16% GC)
  maximum time:     751.103 ms (3.79% GC)
  --------------
  samples:          8
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
  allocs estimate:  87
  --------------
  minimum time:     6.519 s (0.00% GC)
  median time:      6.519 s (0.00% GC)
  mean time:        6.519 s (0.00% GC)
  maximum time:     6.519 s (0.00% GC)
  --------------
  samples:          1
  evals/sample:     1
````



````julia
@benchmark solve(prob,CVODE_BDF(linear_solver=:GMRES))
````


````
BenchmarkTools.Trial: 
  memory estimate:  338.28 MiB
  allocs estimate:  66015
  --------------
  minimum time:     2.028 s (0.50% GC)
  median time:      2.236 s (0.51% GC)
  mean time:        2.175 s (1.14% GC)
  maximum time:     2.262 s (2.33% GC)
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
  memory estimate:  3.82 MiB
  allocs estimate:  50172
  --------------
  minimum time:     1.887 s (0.00% GC)
  median time:      2.082 s (0.00% GC)
  mean time:        2.023 s (0.00% GC)
  maximum time:     2.100 s (0.00% GC)
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
  memory estimate:  5.46 MiB
  allocs estimate:  78553
  --------------
  minimum time:     3.096 s (0.00% GC)
  median time:      3.192 s (0.00% GC)
  mean time:        3.192 s (0.00% GC)
  maximum time:     3.289 s (0.00% GC)
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
Julia Version 1.3.0
Commit 46ce4d7933 (2019-11-26 06:09 UTC)
Platform Info:
  OS: Windows (x86_64-w64-mingw32)
  CPU: Intel(R) Core(TM) i7-8550U CPU @ 1.80GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-6.0.1 (ORCJIT, skylake)
Environment:
  JULIA_EDITOR = "C:\Users\accou\AppData\Local\atom\app-1.42.0\atom.exe"  -a
  JULIA_NUM_THREADS = 4

```

Package Information:

```
Status `~\.julia\dev\DiffEqTutorials\Project.toml`
```
