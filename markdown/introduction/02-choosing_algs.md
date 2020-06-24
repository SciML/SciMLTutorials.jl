---
author: "Chris Rackauckas"
title: "Choosing an ODE Algorithm"
---


While the default algorithms, along with `alg_hints = [:stiff]`, will suffice in most cases, there are times when you may need to exert more control. The purpose of this part of the tutorial is to introduce you to some of the most widely used algorithm choices and when they should be used. The corresponding page of the documentation is the [ODE Solvers](https://docs.juliadiffeq.org/dev/solvers/ode_solve/) page which goes into more depth.

## Diagnosing Stiffness

One of the key things to know for algorithm choices is whether your problem is stiff. Let's take for example the driven Van Der Pol equation:

````julia
using DifferentialEquations, ParameterizedFunctions
van! = @ode_def VanDerPol begin
  dy = μ*((1-x^2)*y - x)
  dx = 1*y
end μ

prob = ODEProblem(van!,[0.0,2.0],(0.0,6.3),1e6)
````


````
ODEProblem with uType Array{Float64,1} and tType Float64. In-place: true
timespan: (0.0, 6.3)
u0: [0.0, 2.0]
````





One indicating factor that should alert you to the fact that this model may be stiff is the fact that the parameter is `1e6`: large parameters generally mean stiff models. If we try to solve this with the default method:

````julia
sol = solve(prob,Tsit5())
````


````
retcode: MaxIters
Interpolation: specialized 4th order "free" interpolation
t: 999977-element Array{Float64,1}:
 0.0
 4.997501249375313e-10
 5.4972513743128435e-9
 3.289919594544218e-8
 9.055581394883546e-8
 1.7309428803584187e-7
 2.79375393394586e-7
 4.149527171475212e-7
 5.807919390815544e-7
 7.81280701490125e-7
 ⋮
 1.8457012081010522
 1.845702696026691
 1.8457041839548325
 1.8457056718857727
 1.845707159819413
 1.8457086477557534
 1.8457101356946952
 1.8457116236362385
 1.8457131115805805
u: 999977-element Array{Array{Float64,1},1}:
 [0.0, 2.0]
 [-0.0009987513736106552, 1.9999999999997504]
 [-0.010904339759596433, 1.9999999999699458]
 [-0.06265556194129239, 1.9999999989523902]
 [-0.1585948892562767, 1.9999999924944207]
 [-0.2700352862461109, 1.9999999746155703]
 [-0.3783197963325601, 1.9999999398563364]
 [-0.47467864703912216, 1.9999998815910678]
 [-0.5499302545937235, 1.999999796115446]
 [-0.6026934372089534, 1.9999996800439757]
 ⋮
 [-0.7770871866226842, 1.8321769350351387]
 [-0.7770880934309836, 1.8321757783565626]
 [-0.7770890004563554, 1.832174621674691]
 [-0.7770899073362528, 1.832173464989294]
 [-0.7770908141915421, 1.832172308300448]
 [-0.777091721022237, 1.8321711516081531]
 [-0.7770926279492066, 1.832169994912486]
 [-0.7770935349724621, 1.8321688382134467]
 [-0.7770944418503102, 1.8321676815108816]
````





Here it shows that maximum iterations were reached. Another thing that can happen is that the solution can return that the solver was unstable (exploded to infinity) or that `dt` became too small. If these happen, the first thing to do is to check that your model is correct. It could very well be that you made an error that causes the model to be unstable!

If the model is the problem, then stiffness could be the reason. We can thus hint to the solver to use an appropriate method:

````julia
sol = solve(prob,alg_hints = [:stiff])
````


````
retcode: Success
Interpolation: specialized 3rd order "free" stiffness-aware interpolation
t: 694-element Array{Float64,1}:
 0.0
 4.997501249375313e-10
 5.454105825317844e-9
 1.8954286226539402e-8
 4.149674379723465e-8
 7.308080698498873e-8
 1.1714649583268228e-7
 1.7481330012839592e-7
 2.4862371241875205e-7
 3.4025555832660864e-7
 ⋮
 5.69767949075004
 5.748994165486137
 5.811844321155623
 5.886853430367259
 5.969502584336209
 6.05645855489726
 6.143414525458311
 6.230370496019361
 6.3
u: 694-element Array{Array{Float64,1},1}:
 [0.0, 2.0]
 [-0.0009987513736106515, 1.9999999999997504]
 [-0.010819454588930516, 1.9999999999704143]
 [-0.036850919195836614, 1.999999999647449]
 [-0.07803540833511657, 1.999999998347307]
 [-0.13124866317048456, 1.9999999950290166]
 [-0.19755036285072544, 1.9999999877524572]
 [-0.272075415352352, 1.999999974149605]
 [-0.35045254633113243, 1.9999999510683737]
 [-0.4264538643666248, 1.9999999153142478]
 ⋮
 [0.6849948021041035, -1.9679959070285558]
 [0.7068255882516246, -1.9322949901761672]
 [0.7369247908644185, -1.8869463160135977]
 [0.7789756893010558, -1.8301403490903476]
 [0.8358041460218815, -1.7634992825515126]
 [0.9131711695745722, -1.6876171241799416]
 [1.0200095610067244, -1.6038403486733988]
 [1.182122272069454, -1.5086434776790882]
 [1.3982811580024197, -1.4194614700844543]
````





Or we can use the default algorithm. By default, DifferentialEquations.jl uses algorithms like `AutoTsit5(Rodas5())` which automatically detect stiffness and switch to an appropriate method once stiffness is known.

````julia
sol = solve(prob)
````


````
retcode: Success
Interpolation: Automatic order switching interpolation
t: 1927-element Array{Float64,1}:
 0.0
 4.997501249375313e-10
 5.4972513743128435e-9
 3.289919594544218e-8
 9.055581394883546e-8
 1.7309428803584187e-7
 2.79375393394586e-7
 4.149527171475212e-7
 5.807919390815544e-7
 7.81280701490125e-7
 ⋮
 6.204648226459174
 6.219556296846657
 6.233842035405188
 6.247504825256131
 6.260547615587699
 6.2729765212211905
 6.2848005774657025
 6.296031197826328
 6.3
u: 1927-element Array{Array{Float64,1},1}:
 [0.0, 2.0]
 [-0.0009987513736106552, 1.9999999999997504]
 [-0.010904339759596433, 1.9999999999699458]
 [-0.06265556194129239, 1.9999999989523902]
 [-0.1585948892562767, 1.9999999924944207]
 [-0.2700352862461109, 1.9999999746155703]
 [-0.3783197963325601, 1.9999999398563364]
 [-0.47467864703912216, 1.9999998815910678]
 [-0.5499302545937235, 1.999999796115446]
 [-0.6026934372089534, 1.9999996800439757]
 ⋮
 [1.1173199781760204, -1.5429755449491689]
 [1.1481789818936488, -1.526090225258629]
 [1.1805062564165951, -1.5094585545821166]
 [1.2143604506793437, -1.4931001271860491]
 [1.2498033047017658, -1.477032218841621]
 [1.286899775928202, -1.461269894062796]
 [1.3257186497126714, -1.4458259322647566]
 [1.3663322749653013, -1.4307111536240857]
 [1.3818948926590964, -1.4252572688140688]
````





Another way to understand stiffness is to look at the solution.

````julia
using Plots; gr()
sol = solve(prob,alg_hints = [:stiff],reltol=1e-6)
plot(sol,denseplot=false)
````


![](figures/02-choosing_algs_5_1.png)



Let's zoom in on the y-axis to see what's going on:

````julia
plot(sol,ylims = (-10.0,10.0))
````


![](figures/02-choosing_algs_6_1.png)



Notice how there are some extreme vertical shifts that occur. These vertical shifts are places where the derivative term is very large, and this is indicative of stiffness. This is an extreme example to highlight the behavior, but this general idea can be carried over to your problem. When in doubt, simply try timing using both a stiff solver and a non-stiff solver and see which is more efficient.

To try this out, let's use BenchmarkTools, a package that let's us relatively reliably time code blocks.

````julia
function lorenz!(du,u,p,t)
    σ,ρ,β = p
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end
u0 = [1.0,0.0,0.0]
p = (10,28,8/3)
tspan = (0.0,100.0)
prob = ODEProblem(lorenz!,u0,tspan,p)
````


````
ODEProblem with uType Array{Float64,1} and tType Float64. In-place: true
timespan: (0.0, 100.0)
u0: [1.0, 0.0, 0.0]
````





And now, let's use the `@btime` macro from benchmark tools to compare the use of non-stiff and stiff solvers on this problem.

````julia
using BenchmarkTools
@btime solve(prob);
````


````
872.884 μs (13109 allocations: 1.42 MiB)
retcode: Success
Interpolation: Automatic order switching interpolation
t: 1294-element Array{Float64,1}:
   0.0
   3.5678604836301404e-5
   0.0003924646531993154
   0.0032624077544510573
   0.009058075635317072
   0.01695646895607931
   0.0276899566248403
   0.041856345938267966
   0.06024040228733675
   0.08368539694547242
   ⋮
  99.39403070915297
  99.47001147494375
  99.54379656909015
  99.614651558349
  99.69093823148101
  99.78733023233721
  99.86114450046736
  99.96115759510786
 100.0
u: 1294-element Array{Array{Float64,1},1}:
 [1.0, 0.0, 0.0]
 [0.9996434557625105, 0.0009988049817849058, 1.781434788799208e-8]
 [0.9961045497425811, 0.010965399721242457, 2.146955365838907e-6]
 [0.9693591634199452, 0.08977060667778931, 0.0001438018342266937]
 [0.9242043615038835, 0.24228912482984957, 0.0010461623302512404]
 [0.8800455868998046, 0.43873645009348244, 0.0034242593451028745]
 [0.8483309877783048, 0.69156288756671, 0.008487623500490047]
 [0.8495036595681027, 1.0145425335433382, 0.01821208597613427]
 [0.9139069079152129, 1.4425597546855036, 0.03669381053327124]
 [1.0888636764765296, 2.052326153029042, 0.07402570506414284]
 ⋮
 [12.999157033749652, 14.10699925404482, 31.74244844521858]
 [11.646131422021162, 7.2855792145502845, 35.365000488215486]
 [7.777555445486692, 2.5166095828739574, 32.030953593541675]
 [4.739741627223412, 1.5919220588229062, 27.249779003951755]
 [3.2351668945618774, 2.3121727966182695, 22.724936101772805]
 [3.310411964698304, 4.28106626744641, 18.435441144016366]
 [4.527117863517627, 6.895878639772805, 16.58544600757436]
 [8.043672261487556, 12.711555298531689, 18.12537420595938]
 [9.97537965430362, 15.143884806010783, 21.00643286956427]
````



````julia
@btime solve(prob,alg_hints = [:stiff]);
````


````
9.639 ms (63542 allocations: 2.78 MiB)
retcode: Success
Interpolation: specialized 3rd order "free" stiffness-aware interpolation
t: 2353-element Array{Float64,1}:
   0.0
   3.5678604836301404e-5
   0.0003924646531993154
   0.0015292937978620778
   0.003162305724282332
   0.005285973552599222
   0.008023688976652094
   0.01144101871159114
   0.0156693459080261
   0.020826710718121224
   ⋮
  99.64643805711006
  99.69715002361579
  99.75331122793506
  99.79981678120402
  99.84132760768253
  99.88283843416104
  99.92434926063956
  99.96697093084641
 100.0
u: 2353-element Array{Array{Float64,1},1}:
 [1.0, 0.0, 0.0]
 [0.9996434557625105, 0.0009988049817849047, 1.7814347887985208e-8]
 [0.9961045497425969, 0.0109653997212298, 2.146955365112677e-6]
 [0.9851473616483439, 0.04246652425810003, 3.21927130421189e-5]
 [0.9702414465725462, 0.08706126023105658, 0.00013525574346441506]
 [0.9522854546465404, 0.14396668240059424, 0.00036967772967708476]
 [0.9314326271963391, 0.21565548792970338, 0.0008289976215629184]
 [0.9088641247467208, 0.3027798780735744, 0.001632872227230375]
 [0.8859643488922837, 0.4074631051549474, 0.002954345942806542]
 [0.8650679703454024, 0.5313629257625292, 0.0050183974807311285]
 ⋮
 [12.806695275278894, 12.494656408454457, 33.00345553979429]
 [11.743080430598853, 8.043506580216379, 34.89282227510988]
 [9.11982468203129, 3.968709691943006, 33.23949044247742]
 [6.836119847459929, 2.3677711311840715, 30.424290289803416]
 [5.229488956523552, 1.9977361551169743, 27.73491536637819]
 [4.156928665122608, 2.2016007548657157, 25.204702582122852]
 [3.5806565170126663, 2.724678757457538, 22.930442334447942]
 [3.4202075264870198, 3.4950829296305566, 20.89909215126914]
 [3.5509784304871257, 4.260058718117392, 19.561783277650832]
````





In this particular case, we can see that non-stiff solvers get us to the solution much more quickly.

## The Recommended Methods

When picking a method, the general rules are as follows:

- Higher order is more efficient at lower tolerances, lower order is more efficient at higher tolerances
- Adaptivity is essential in most real-world scenarios
- Runge-Kutta methods do well with non-stiff equations, Rosenbrock methods do well with small stiff equations, BDF methods do well with large stiff equations

While there are always exceptions to the rule, those are good guiding principles. Based on those, a simple way to choose methods is:

- The default is `Tsit5()`, a non-stiff Runge-Kutta method of Order 5
- If you use low tolerances (`1e-8`), try `Vern7()` or `Vern9()`
- If you use high tolerances, try `BS3()`
- If the problem is stiff, try `Rosenbrock23()`, `Rodas5()`, or `CVODE_BDF()`
- If you don't know, use `AutoTsit5(Rosenbrock23())` or `AutoVern9(Rodas5())`.

(This is a simplified version of the default algorithm chooser)

## Comparison to other Software

If you are familiar with MATLAB, SciPy, or R's DESolve, here's a quick translation start to have transfer your knowledge over.

- `ode23` -> `BS3()`
- `ode45`/`dopri5` -> `DP5()`, though in most cases `Tsit5()` is more efficient
- `ode23s` -> `Rosenbrock23()`, though in most cases `Rodas4()` is more efficient
- `ode113` -> `VCABM()`, though in many cases `Vern7()` is more efficient
- `dop853` -> `DP8()`, though in most cases `Vern7()` is more efficient
- `ode15s`/`vode` -> `QNDF()`, though in many cases `CVODE_BDF()`, `Rodas4()`
  or `radau()` are more efficient
- `ode23t` -> `Trapezoid()` for efficiency and `GenericTrapezoid()` for robustness
- `ode23tb` -> `TRBDF2`
- `lsoda` -> `lsoda()` (requires `]add LSODA; using LSODA`)
- `ode15i` -> `IDA()`, though in many cases `Rodas4()` can handle the DAE and is
  significantly more efficient


## Appendix

 This tutorial is part of the DiffEqTutorials.jl repository, found at: <https://github.com/JuliaDiffEq/DiffEqTutorials.jl>

To locally run this tutorial, do the following commands:
```
using DiffEqTutorials
DiffEqTutorials.weave_file("introduction","02-choosing_algs.jmd")
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
Status `/builds/JuliaGPU/DiffEqTutorials.jl/tutorials/introduction/Project.toml`
[6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf] BenchmarkTools 0.5.0
[0c46a032-eb83-5123-abaf-570d42b7fbaa] DifferentialEquations 6.14.0
[65888b18-ceab-5e60-b2b9-181511a3b968] ParameterizedFunctions 5.3.0
[91a5bcdd-55d7-5caf-9e0b-520d859cae80] Plots 1.4.3
[90137ffa-7385-5640-81b9-e52037218182] StaticArrays 0.12.3
[c3572dad-4567-51f8-b174-8c6c989267f4] Sundials 4.2.3
[37e2e46d-f89d-539d-b4ee-838fcccc9c8e] LinearAlgebra
```
