---
author: "Chris Rackauckas"
title: "An Intro to DifferentialEquations.jl"
---


## Basic Introduction Via Ordinary Differential Equations

This notebook will get you started with DifferentialEquations.jl by introducing you to the functionality for solving ordinary differential equations (ODEs). The corresponding documentation page is the [ODE tutorial](https://docs.juliadiffeq.org/dev/tutorials/ode_example/). While some of the syntax may be different for other types of equations, the same general principles hold in each case. Our goal is to give a gentle and thorough introduction that highlights these principles in a way that will help you generalize what you have learned.

### Background

If you are new to the study of differential equations, it can be helpful to do a quick background read on [the definition of ordinary differential equations](https://en.wikipedia.org/wiki/Ordinary_differential_equation). We define an ordinary differential equation as an equation which describes the way that a variable $u$ changes, that is

$$u' = f(u,p,t)$$

where $p$ are the parameters of the model, $t$ is the time variable, and $f$ is the nonlinear model of how $u$ changes. The initial value problem also includes the information about the starting value:

$$u(t_0) = u_0$$

Together, if you know the starting value and you know how the value will change with time, then you know what the value will be at any time point in the future. This is the intuitive definition of a differential equation.

### First Model: Exponential Growth

Our first model will be the canonical exponential growth model. This model says that the rate of change is proportional to the current value, and is this:

$$u' = au$$

where we have a starting value $u(0)=u_0$. Let's say we put 1 dollar into Bitcoin which is increasing at a rate of $98\%$ per year. Then calling now $t=0$ and measuring time in years, our model is:

$$u' = 0.98u$$

and $u(0) = 1.0$. We encode this into Julia by noticing that, in this setup, we match the general form when

````julia
f(u,p,t) = 0.98u
````


````
f (generic function with 1 method)
````





with $ u_0 = 1.0 $. If we want to solve this model on a time span from `t=0.0` to `t=1.0`, then we define an `ODEProblem` by specifying this function `f`, this initial condition `u0`, and this time span as follows:

````julia
using DifferentialEquations
f(u,p,t) = 0.98u
u0 = 1.0
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
````


````
ODEProblem with uType Float64 and tType Float64. In-place: false
timespan: (0.0, 1.0)
u0: 1.0
````





To solve our `ODEProblem` we use the command `solve`.

````julia
sol = solve(prob)
````


````
retcode: Success
Interpolation: Automatic order switching interpolation
t: 5-element Array{Float64,1}:
 0.0
 0.10042494449239292
 0.35218603951893646
 0.6934436028208104
 1.0
u: 5-element Array{Float64,1}:
 1.0
 1.1034222047865465
 1.4121908848175448
 1.9730384275622996
 2.664456142481451
````





and that's it: we have succesfully solved our first ODE!

#### Analyzing the Solution

Of course, the solution type is not interesting in and of itself. We want to understand the solution! The documentation page which explains in detail the functions for analyzing the solution is the [Solution Handling](https://docs.juliadiffeq.org/dev/basics/solution/) page. Here we will describe some of the basics. You can plot the solution using the plot recipe provided by [Plots.jl](http://docs.juliaplots.org/dev/):

````julia
using Plots; gr()
plot(sol)
````


![](figures/01-ode_introduction_4_1.png)



From the picture we see that the solution is an exponential curve, which matches our intuition. As a plot recipe, we can annotate the result using any of the [Plots.jl attributes](http://docs.juliaplots.org/dev/attributes/). For example:

````julia
plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
````


![](figures/01-ode_introduction_5_1.png)



Using the mutating `plot!` command we can add other pieces to our plot. For this ODE we know that the true solution is $u(t) = u_0 exp(at)$, so let's add some of the true solution to our plot:

````julia
plot!(sol.t, t->1.0*exp(0.98t),lw=3,ls=:dash,label="True Solution!")
````


![](figures/01-ode_introduction_6_1.png)



In the previous command I demonstrated `sol.t`, which grabs the array of time points that the solution was saved at:

````julia
sol.t
````


````
5-element Array{Float64,1}:
 0.0
 0.10042494449239292
 0.35218603951893646
 0.6934436028208104
 1.0
````





We can get the array of solution values using `sol.u`:

````julia
sol.u
````


````
5-element Array{Float64,1}:
 1.0
 1.1034222047865465
 1.4121908848175448
 1.9730384275622996
 2.664456142481451
````





`sol.u[i]` is the value of the solution at time `sol.t[i]`. We can compute arrays of functions of the solution values using standard comprehensions, like:

````julia
[t+u for (u,t) in tuples(sol)]
````


````
5-element Array{Float64,1}:
 1.0
 1.2038471492789395
 1.7643769243364813
 2.66648203038311
 3.664456142481451
````





However, one interesting feature is that, by default, the solution is a continuous function. If we check the print out again:

````julia
sol
````


````
retcode: Success
Interpolation: Automatic order switching interpolation
t: 5-element Array{Float64,1}:
 0.0
 0.10042494449239292
 0.35218603951893646
 0.6934436028208104
 1.0
u: 5-element Array{Float64,1}:
 1.0
 1.1034222047865465
 1.4121908848175448
 1.9730384275622996
 2.664456142481451
````





you see that it says that the solution has a order changing interpolation. The default algorithm automatically switches between methods in order to handle all types of problems. For non-stiff equations (like the one we are solving), it is a continuous function of 4th order accuracy. We can call the solution as a function of time `sol(t)`. For example, to get the value at `t=0.45`, we can use the command:

````julia
sol(0.45)
````


````
1.554261048055312
````





#### Controlling the Solver

DifferentialEquations.jl has a common set of solver controls among its algorithms which can be found [at the Common Solver Options](https://docs.juliadiffeq.org/dev/basics/common_solver_opts/) page. We will detail some of the most widely used options.

The most useful options are the tolerances `abstol` and `reltol`. These tell the internal adaptive time stepping engine how precise of a solution you want. Generally, `reltol` is the relative accuracy while `abstol` is the accuracy when `u` is near zero. These tolerances are local tolerances and thus are not global guarantees. However, a good rule of thumb is that the total solution accuracy is 1-2 digits less than the relative tolerances. Thus for the defaults `abstol=1e-6` and `reltol=1e-3`, you can expect a global accuracy of about 1-2 digits. If we want to get around 6 digits of accuracy, we can use the commands:

````julia
sol = solve(prob,abstol=1e-8,reltol=1e-8)
````


````
retcode: Success
Interpolation: Automatic order switching interpolation
t: 9-element Array{Float64,1}:
 0.0
 0.04127492324135852
 0.14679917846877366
 0.28631546412766684
 0.4381941361169628
 0.6118924302028597
 0.7985659100883337
 0.9993516479536952
 1.0
u: 9-element Array{Float64,1}:
 1.0
 1.0412786454705882
 1.1547261252949712
 1.3239095703537043
 1.5363819257509728
 1.8214895157178692
 2.1871396448296223
 2.662763824115295
 2.664456241933517
````





Now we can see no visible difference against the true solution:


````julia
plot(sol)
plot!(sol.t, t->1.0*exp(0.98t),lw=3,ls=:dash,label="True Solution!")
````


![](figures/01-ode_introduction_13_1.png)



Notice that by decreasing the tolerance, the number of steps the solver had to take was `9` instead of the previous `5`. There is a trade off between accuracy and speed, and it is up to you to determine what is the right balance for your problem.

Another common option is to use `saveat` to make the solver save at specific time points. For example, if we want the solution at an even grid of `t=0.1k` for integers `k`, we would use the command:

````julia
sol = solve(prob,saveat=0.1)
````


````
retcode: Success
Interpolation: 1st order linear
t: 11-element Array{Float64,1}:
 0.0
 0.1
 0.2
 0.3
 0.4
 0.5
 0.6
 0.7
 0.8
 0.9
 1.0
u: 11-element Array{Float64,1}:
 1.0
 1.102962785129292
 1.2165269512238264
 1.341783821227542
 1.4799379510586077
 1.632316207054161
 1.8003833264983584
 1.9857565541588758
 2.1902158127997695
 2.415725742084496
 2.664456142481451
````





Notice that when `saveat` is used the continuous output variables are no longer saved and thus `sol(t)`, the interpolation, is only first order. We can save at an uneven grid of points by passing a collection of values to `saveat`. For example:

````julia
sol = solve(prob,saveat=[0.2,0.7,0.9])
````


````
retcode: Success
Interpolation: 1st order linear
t: 3-element Array{Float64,1}:
 0.2
 0.7
 0.9
u: 3-element Array{Float64,1}:
 1.2165269512238264
 1.9857565541588758
 2.415725742084496
````





If we need to reduce the amount of saving, we can also turn off the continuous output directly via `dense=false`:

````julia
sol = solve(prob,dense=false)
````


````
retcode: Success
Interpolation: 1st order linear
t: 5-element Array{Float64,1}:
 0.0
 0.10042494449239292
 0.35218603951893646
 0.6934436028208104
 1.0
u: 5-element Array{Float64,1}:
 1.0
 1.1034222047865465
 1.4121908848175448
 1.9730384275622996
 2.664456142481451
````





and to turn off all intermediate saving we can use `save_everystep=false`:

````julia
sol = solve(prob,save_everystep=false)
````


````
retcode: Success
Interpolation: 1st order linear
t: 2-element Array{Float64,1}:
 0.0
 1.0
u: 2-element Array{Float64,1}:
 1.0
 2.664456142481451
````





If we want to solve and only save the final value, we can even set `save_start=false`.

````julia
sol = solve(prob,save_everystep=false,save_start = false)
````


````
retcode: Success
Interpolation: 1st order linear
t: 1-element Array{Float64,1}:
 1.0
u: 1-element Array{Float64,1}:
 2.664456142481451
````





Note that similarly on the other side there is `save_end=false`.

More advanced saving behaviors, such as saving functionals of the solution, are handled via the `SavingCallback` in the [Callback Library](https://docs.juliadiffeq.org/dev/features/callback_library/#saving_callback-1) which will be addressed later in the tutorial.

#### Choosing Solver Algorithms

There is no best algorithm for numerically solving a differential equation. When you call `solve(prob)`, DifferentialEquations.jl makes a guess at a good algorithm for your problem, given the properties that you ask for (the tolerances, the saving information, etc.). However, in many cases you may want more direct control. A later notebook will help introduce the various *algorithms* in DifferentialEquations.jl, but for now let's introduce the *syntax*.

The most crucial determining factor in choosing a numerical method is the stiffness of the model. Stiffness is roughly characterized by a Jacobian `f` with large eigenvalues. That's quite mathematical, and we can think of it more intuitively: if you have big numbers in `f` (like parameters of order `1e5`), then it's probably stiff. Or, as the creator of the MATLAB ODE Suite, Lawrence Shampine, likes to define it, if the standard algorithms are slow, then it's stiff. We will go into more depth about diagnosing stiffness in a later tutorial, but for now note that if you believe your model may be stiff, you can hint this to the algorithm chooser via `alg_hints = [:stiff]`.

````julia
sol = solve(prob,alg_hints=[:stiff])
````


````
retcode: Success
Interpolation: specialized 3rd order "free" stiffness-aware interpolation
t: 8-element Array{Float64,1}:
 0.0
 0.05653299582822294
 0.17270731152826024
 0.3164602871490142
 0.5057500163821153
 0.7292241858994543
 0.9912975001018789
 1.0
u: 8-element Array{Float64,1}:
 1.0
 1.0569657840332976
 1.1844199383303913
 1.3636037723365293
 1.6415399686182572
 2.0434491434754793
 2.641825616057761
 2.6644526430553817
````





Stiff algorithms have to solve implicit equations and linear systems at each step so they should only be used when required.

If we want to choose an algorithm directly, you can pass the algorithm type after the problem as `solve(prob,alg)`. For example, let's solve this problem using the `Tsit5()` algorithm, and just for show let's change the relative tolerance to `1e-6` at the same time:

````julia
sol = solve(prob,Tsit5(),reltol=1e-6)
````


````
retcode: Success
Interpolation: specialized 4th order "free" interpolation
t: 10-element Array{Float64,1}:
 0.0
 0.028970819746309166
 0.10049147151547619
 0.19458908698515082
 0.3071725081673423
 0.43945421453622546
 0.5883434923759523
 0.7524873357619015
 0.9293021330536031
 1.0
u: 10-element Array{Float64,1}:
 1.0
 1.0287982807225062
 1.1034941463604806
 1.2100931078233779
 1.351248605624241
 1.538280340326815
 1.7799346012651116
 2.090571742234628
 2.486102171447025
 2.6644562434913377
````





### Systems of ODEs: The Lorenz Equation

Now let's move to a system of ODEs. The [Lorenz equation](https://en.wikipedia.org/wiki/Lorenz_system) is the famous "butterfly attractor" that spawned chaos theory. It is defined by the system of ODEs:

$$
\begin{align}
\frac{dx}{dt} &= \sigma (y - x)\\
\frac{dy}{dt} &= x (\rho - z) -y\\
\frac{dz}{dt} &= xy - \beta z
\end{align}
$$

To define a system of differential equations in DifferentialEquations.jl, we define our `f` as a vector function with a vector initial condition. Thus, for the vector `u = [x,y,z]'`, we have the derivative function:

````julia
function lorenz!(du,u,p,t)
    σ,ρ,β = p
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end
````


````
lorenz! (generic function with 1 method)
````





Notice here we used the in-place format which writes the output to the preallocated vector `du`. For systems of equations the in-place format is faster. We use the initial condition $u_0 = [1.0,0.0,0.0]$ as follows:

````julia
u0 = [1.0,0.0,0.0]
````


````
3-element Array{Float64,1}:
 1.0
 0.0
 0.0
````





Lastly, for this model we made use of the parameters `p`. We need to set this value in the `ODEProblem` as well. For our model we want to solve using the parameters $\sigma = 10$, $\rho = 28$, and $\beta = 8/3$, and thus we build the parameter collection:

````julia
p = (10,28,8/3) # we could also make this an array, or any other type!
````


````
(10, 28, 2.6666666666666665)
````





Now we generate the `ODEProblem` type. In this case, since we have parameters, we add the parameter values to the end of the constructor call. Let's solve this on a time span of `t=0` to `t=100`:

````julia
tspan = (0.0,100.0)
prob = ODEProblem(lorenz!,u0,tspan,p)
````


````
ODEProblem with uType Array{Float64,1} and tType Float64. In-place: true
timespan: (0.0, 100.0)
u0: [1.0, 0.0, 0.0]
````





Now, just as before, we solve the problem:

````julia
sol = solve(prob)
````


````
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





The same solution handling features apply to this case. Thus `sol.t` stores the time points and `sol.u` is an array storing the solution at the corresponding time points.

However, there are a few extra features which are good to know when dealing with systems of equations. First of all, `sol` also acts like an array. `sol[i]` returns the solution at the `i`th time point.

````julia
sol.t[10],sol[10]
````


````
(0.08368539694547242, [1.0888636764765296, 2.052326153029042, 0.07402570506
414284])
````





Additionally, the solution acts like a matrix where `sol[j,i]` is the value of the `j`th variable at time `i`:

````julia
sol[2,10]
````


````
2.052326153029042
````





We can get a real matrix by performing a conversion:

````julia
A = Array(sol)
````


````
3×1294 Array{Float64,2}:
 1.0  0.999643     0.996105    0.969359     …   4.52712   8.04367   9.97538
 0.0  0.000998805  0.0109654   0.0897706        6.89588  12.7116   15.1439
 0.0  1.78143e-8   2.14696e-6  0.000143802     16.5854   18.1254   21.0064
````





This is the same as sol, i.e. `sol[i,j] = A[i,j]`, but now it's a true matrix. Plotting will by default show the time series for each variable:

````julia
plot(sol)
````


![](figures/01-ode_introduction_29_1.png)



If we instead want to plot values against each other, we can use the `vars` command. Let's plot variable `1` against variable `2` against variable `3`:

````julia
plot(sol,vars=(1,2,3))
````


![](figures/01-ode_introduction_30_1.png)



This is the classic Lorenz attractor plot, where the `x` axis is `u[1]`, the `y` axis is `u[2]`, and the `z` axis is `u[3]`. Note that the plot recipe by default uses the interpolation, but we can turn this off:

````julia
plot(sol,vars=(1,2,3),denseplot=false)
````


![](figures/01-ode_introduction_31_1.png)



Yikes! This shows how calculating the continuous solution has saved a lot of computational effort by computing only a sparse solution and filling in the values! Note that in vars, `0=time`, and thus we can plot the time series of a single component like:

````julia
plot(sol,vars=(0,2))
````


![](figures/01-ode_introduction_32_1.png)



## Internal Types

The last basic user-interface feature to explore is the choice of types. DifferentialEquations.jl respects your input types to determine the internal types that are used. Thus since in the previous cases, when we used `Float64` values for the initial condition, this meant that the internal values would be solved using `Float64`. We made sure that time was specified via `Float64` values, meaning that time steps would utilize 64-bit floats as well. But, by simply changing these types we can change what is used internally.

As a quick example, let's say we want to solve an ODE defined by a matrix. To do this, we can simply use a matrix as input.

````julia
A  = [1. 0  0 -5
      4 -2  4 -3
     -4  0  0  1
      5 -2  2  3]
u0 = rand(4,2)
tspan = (0.0,1.0)
f(u,p,t) = A*u
prob = ODEProblem(f,u0,tspan)
sol = solve(prob)
````


````
retcode: Success
Interpolation: Automatic order switching interpolation
t: 10-element Array{Float64,1}:
 0.0
 0.04362622558308152
 0.1160883077551463
 0.198706031655623
 0.302460390063917
 0.422394977611238
 0.5667359726754638
 0.7214657872172093
 0.8827242841053362
 1.0
u: 10-element Array{Array{Float64,2},1}:
 [0.87711279044473 0.648155124773695; 0.4185251824823173 0.4795544197424175
; 0.8588907285543594 0.7681300386777359; 0.7492603809206482 0.9365850837133
614]
 [0.7133009857038476 0.43680898925630374; 0.5372566450939396 0.518983631535
743; 0.7587546344777232 0.7194919100412043; 1.0712454648804757 1.2177182869
285645]
 [0.2703542199110978 -0.06803702196602579; 0.5254513469705496 0.40670688578
53301; 0.7070459449839606 0.7653541586971416; 1.5679725141308694 1.63770155
94729114]
 [-0.4867282070513183 -0.8636888718534027; 0.23559398189709813 0.0569115177
56903285; 0.8856450759843708 1.0649433813690432; 2.0432809573730717 2.01420
74463982182]
 [-1.774502474990292 -2.1407369809782906; -0.419653363200437 -0.58236826028
83294; 1.577846962488419 1.903862093611193; 2.436088704714329 2.27118015473
4261]
 [-3.597388117824066 -3.8544402213409383; -1.2928122016065546 -1.3077355725
858306; 3.156782405783615 3.6071955140306886; 2.4971713194954352 2.16981765
7701055]
 [-5.902080517788761 -5.874777466099962; -1.8825558763292398 -1.55802002101
11457; 6.2293517151241184 6.6937167156578425; 1.8361398961595583 1.32940217
61866412]
 [-7.7797665277364 -7.264466183762919; -0.9721424208246652 -0.1660506116421
5185; 10.680554468548427 10.896484684237814; 0.04684471870831475 -0.5901939
4672959]
 [-7.9384715596770326 -6.761557663922413; 2.8181988381633634 4.094791191140
576; 15.65606642663034 15.226865654511265; -3.0919976794633977 -3.732717451
0842998]
 [-6.112849595971283 -4.422604134522122; 7.888503767206559 9.35331142548942
7; 18.495323068968812 17.324100940184252; -6.0962867491929025 -6.6113770450
00404]
````





There is no real difference from what we did before, but now in this case `u0` is a `4x2` matrix. Because of that, the solution at each time point is matrix:

````julia
sol[3]
````


````
4×2 Array{Float64,2}:
 0.270354  -0.068037
 0.525451   0.406707
 0.707046   0.765354
 1.56797    1.6377
````





In DifferentialEquations.jl, you can use any type that defines `+`, `-`, `*`, `/`, and has an appropriate `norm`. For example, if we want arbitrary precision floating point numbers, we can change the input to be a matrix of `BigFloat`:

````julia
big_u0 = big.(u0)
````


````
4×2 Array{BigFloat,2}:
 0.877113  0.648155
 0.418525  0.479554
 0.858891  0.76813
 0.74926   0.936585
````





and we can solve the `ODEProblem` with arbitrary precision numbers by using that initial condition:

````julia
prob = ODEProblem(f,big_u0,tspan)
sol = solve(prob)
````


````
retcode: Success
Interpolation: Automatic order switching interpolation
t: 5-element Array{Float64,1}:
 0.0
 0.1591988729994351
 0.42145884936170375
 0.7265465144089278
 1.0
u: 5-element Array{Array{BigFloat,2},1}:
 [0.877112790444730006100826358306221663951873779296875 0.64815512477369496
8354220691253431141376495361328125; 0.4185251824823172928091707944986410439
014434814453125 0.479554419742417525895916696754284203052520751953125; 0.85
88907285543594394283672954770736396312713623046875 0.7681300386777358824019
75686661899089813232421875; 0.749260380920648172420328592124860733747482299
8046875 0.9365850837133613548957100647385232150554656982421875]
 [-0.0923675895118209726858559064285547097937917784845450810670805094498519
4896493393 -0.4555880976277627736520821237667745555164157309620538508777274
886943193034761816; 0.40678231500074278436156127119239115622189474797022885
51905204158435953278469669 0.2490568380860821843113046822081906510941680513
740596141777668737598626271887532; 0.76402451173215624838314041001513944660
06296584664396442917563513653529936126167 0.8849008997505133075599969394301
084819910753970699141289063253860806569094754233; 1.83069363297601593375842
5215972343229600593099824484033897026198793863407744896 1.85017854027339616
6965148864745662051617288778756781117380023824095037718493349]
 [-3.5823365048149769227574277854304837419015430380595116481740726291653071
93319024 -3.840676613605395318512603044797434953326537442887131642617122746
949069423823231; -1.2865600436511164015295701639508001515068654660280563396
67600761287545727690198 -1.303150853867293543243964706398604556678605123620
522973615073053591919747594387; 3.14100309531303263419691880602697698879283
424118196564913962618710268599463537 3.590757113340202210487760458127062370
407226791075173043955645599432917020631252; 2.49864895604633757270346546183
435354426649157357548799623507433576015774806633 2.172546795641131989336281
078779949542430891838632625475527493586327546061875557]
 [-7.8195889118803138341392808881925757232633922057699167923446709808344137
15145138 -7.285405277362908432524610568267914992620156820040618253300242329
688062400473336; -0.9025641552293431534888679514928789482784663626125975079
085350483769141027875954 -0.08010426992162881863730294193821232219275235022
558494019074974384471491678645119; 10.8391244854581860158323435143376465561
7610090729439949239262261017902218099341 11.0411465664948638679480093508717
005442849361976267711488857882137483836659603; -0.0323293928668488912281797
3170768368034037724815747165095293163049681809166471707 -0.6719185528054829
245800886815952865085865592535451555086444394937735125518096503]
 [-6.1128405641260526683797771872769747417066412710775771651223370544121817
57944277 -4.422590595604050898737626800778981378609671465013765837598765648
605902588746339; 7.88854818599766731500468489349022579423626585322494137260
7729031073926203588125 9.35335739385370920785040150812319439409722769019170
4526372416968166100207805471; 18.495354222771303036145571872931312218347021
43587475889142007183791441372520976 17.324126203338848201734640663266807359
57083486455966514339234912184938889358224; -6.09631052112275135198511024486
4169843315187903250496019935005666328183144359227 -6.6114005323231647323642
42379400540485622344865810075655164017959589735434769853]
````



````julia
sol[1,3]
````


````
-3.582336504814976922757427785430483741901543038059511648174072629165307193
319024
````





To really make use of this, we would want to change `abstol` and `reltol` to be small! Notice that the type for "time" is different than the type for the dependent variables, and this can be used to optimize the algorithm via keeping multiple precisions. We can convert time to be arbitrary precision as well by defining our time span with `BigFloat` variables:

````julia
prob = ODEProblem(f,big_u0,big.(tspan))
sol = solve(prob)
````


````
retcode: Success
Interpolation: Automatic order switching interpolation
t: 5-element Array{BigFloat,1}:
 0.0
 0.159198872999435114729916544714062297724799128281135417471843904364004476
7541098
 0.421458849361703734820659565440821911260632026512155717037377648458182961
8098184
 0.726546514408927875213921471666014963503247797752019598359924152024986153
471474
 1.0
u: 5-element Array{Array{BigFloat,2},1}:
 [0.877112790444730006100826358306221663951873779296875 0.64815512477369496
8354220691253431141376495361328125; 0.4185251824823172928091707944986410439
014434814453125 0.479554419742417525895916696754284203052520751953125; 0.85
88907285543594394283672954770736396312713623046875 0.7681300386777358824019
75686661899089813232421875; 0.749260380920648172420328592124860733747482299
8046875 0.9365850837133613548957100647385232150554656982421875]
 [-0.0923675895118210231055011453676776838211672640999315512494469056879190
1097988876 -0.4555880976277628265837295597611388961911393996654414866355482
046497869295957339; 0.40678231500074276462623389208406258407614039122057954
97141569155623856639377487 0.2490568380860821606912057838341378757231183125
909646878292175104302915059150027; 0.76402451173215626038113506231406449818
43994977516280543757793783797494205607762 0.8849008997505133275871433884869
757947716373694704293308673732633410026265981635; 1.83069363297601596508573
0006561857004497591588343678324914712317311164359507366 1.85017854027339619
1746159009854044703485896462959132091709015883958156881826858]
 [-3.5823365048149771764627265693404464175701888420874418696805751324005640
0597125 -3.8406766136053955505622631493991727353619561642155650502334609622
95985385929521; -1.28656004365111650708231258334173355010735849950925999519
3951463370061295661023 -1.3031508538672936207497721481564715184790423639144
04504027096220157609411329999; 3.141003095313032899776834657270664002695490
86938356500800052949656931919309261 3.5907571133402024872296470203126786680
13934649381698414296696099301423899993325; 2.498648956046337548073881029741
516101987499974871958709280833061756141567495227 2.172546795641131943601114
33750177397388193841649110774733843774015281626495788]
 [-7.8195889118803141005684586148707143571824602209302298256686372240505300
02553878 -7.285405277362908569108417647150961625498578839788882756002407314
100731400433594; -0.9025641552293426670984401011485822693313597656022589244
852032747045119557979202 -0.08010426992162822026537369804494449106835185818
289339671524520605373532631194854; 10.8391244854581871029199799321783353812
7735500282249480081612731239739284042043 11.0411465664948648584439275684863
3370754631283353913643448996522695639261968432; -0.032329392866849437852517
19038667285576741291958191912506136601138006413476437061 -0.671918552805483
4882077647927053171590480608849694235190277839109590675617442996]
 [-6.1128405641260533447474032674983477034445153280319954613661204021350939
59925346 -4.422590595604051693502211921378968322038450383850098594971118091
580362582749123; 7.88854818599766587055208500658714277444266448611781133602
3001893034398631442701 9.35335739385370774419817594590723711660314634369537
874905228756987779426284309; 18.4953542227713025266905480533282235651606111
2270373743217512524132563939116307 17.3241262033388478942316788030684905382
8375660156117945968014758418409724985166; -6.096310521122750584834262379131
216524430537477311758513870662938226329599561962 -6.61140053232316401056337
793089207392651245288755421015400742852181146606307065]
````





Let's end by showing a more complicated use of types. For small arrays, it's usually faster to do operations on static arrays via the package [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl). The syntax is similar to that of normal arrays, but for these special arrays we utilize the `@SMatrix` macro to indicate we want to create a static array.

````julia
using StaticArrays
A  = @SMatrix [ 1.0  0.0 0.0 -5.0
                4.0 -2.0 4.0 -3.0
               -4.0  0.0 0.0  1.0
                5.0 -2.0 2.0  3.0]
u0 = @SMatrix rand(4,2)
tspan = (0.0,1.0)
f(u,p,t) = A*u
prob = ODEProblem(f,u0,tspan)
sol = solve(prob)
````


````
retcode: Success
Interpolation: Automatic order switching interpolation
t: 10-element Array{Float64,1}:
 0.0
 0.04143510784247056
 0.11037727247289317
 0.19442114911076813
 0.294154149040499
 0.4194089557424268
 0.5557734361924471
 0.7022784962178806
 0.8596294137108398
 1.0
u: 10-element Array{StaticArrays.SArray{Tuple{4,2},Float64,2,8},1}:
 [0.14302874818208688 0.09656371121132623; 0.3180117177445976 0.23059672616
447568; 0.6693883576625506 0.5126094242676411; 0.976403745862638 0.15714287
317992515]
 [-0.07579022363146071 0.060953213615606364; 0.2803398025770371 0.283452839
6246962; 0.7072731293919325 0.5071463659068303; 1.1479757279826535 0.217822
1721187863]
 [-0.5350882134737153 -0.028585044522144457; 0.12156004860835715 0.33218717
440487794; 0.8766521116749152 0.5200114749761936; 1.3879211124192634 0.3064
614957012348]
 [-1.2366290472943684 -0.18402881071451532; -0.18681033159545637 0.33660728
16434097; 1.2961414355642724 0.5838056336016412; 1.5806927952870558 0.38653
76385846703]
 [-2.2165091595904505 -0.41975457543982486; -0.6048196391110293 0.292341937
56220095; 2.1426504831986932 0.7439476129613601; 1.6255387087556 0.42894008
77639607]
 [-3.523100995503733 -0.7528405502774808; -0.9535096505021803 0.23096500455
880228; 3.7708035611323485 1.089175060571188; 1.325679075031437 0.379727985
13848]
 [-4.733612877146546 -1.0763765269952008; -0.6799723079418734 0.27918729720
585617; 6.170454719382499 1.632017716645011; 0.4595120877997112 0.171698272
4641927]
 [-5.267247474429787 -1.2340358111824612; 0.9394370482731442 0.630102795929
8159; 9.112335848190057 2.3197932634988367; -1.1591177012767833 -0.24628006
96939955]
 [-4.196747963144789 -0.9760666707628697; 4.76667626964911 1.52041282418377
42; 11.833330613760527 2.9566654611274843; -3.6392391749746755 -0.901273960
0148229]
 [-1.12976078496356 -0.1901029571172682; 10.185902867690352 2.8266373535690
663; 12.742225420859654 3.137299258974818; -6.2802826464502335 -1.601805205
639891]
````



````julia
sol[3]
````


````
4×2 StaticArrays.SArray{Tuple{4,2},Float64,2,8} with indices SOneTo(4)×SOne
To(2):
 -0.535088  -0.028585
  0.12156    0.332187
  0.876652   0.520011
  1.38792    0.306461
````





## Conclusion

These are the basic controls in DifferentialEquations.jl. All equations are defined via a problem type, and the `solve` command is used with an algorithm choice (or the default) to get a solution. Every solution acts the same, like an array `sol[i]` with `sol.t[i]`, and also like a continuous function `sol(t)` with a nice plot command `plot(sol)`. The Common Solver Options can be used to control the solver for any equation type. Lastly, the types used in the numerical solving are determined by the input types, and this can be used to solve with arbitrary precision and add additional optimizations (this can be used to solve via GPUs for example!). While this was shown on ODEs, these techniques generalize to other types of equations as well.

````
Error: ArgumentError: Package DiffEqTutorials not found in current path:
- Run `import Pkg; Pkg.add("DiffEqTutorials")` to install the DiffEqTutoria
ls package.
````



````
Error: UndefVarError: DiffEqTutorials not defined
````


