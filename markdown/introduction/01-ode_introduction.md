---
author: "Chris Rackauckas"
title: "An Intro to DifferentialEquations.jl"
---


## Basic Introduction Via Ordinary Differential Equations

This notebook will get you started with DifferentialEquations.jl by introducing you to the functionality for solving ordinary differential equations (ODEs). The corresponding documentation page is the [ODE tutorial](https://docs.sciml.ai/dev/tutorials/ode_example/). While some of the syntax may be different for other types of equations, the same general principles hold in each case. Our goal is to give a gentle and thorough introduction that highlights these principles in a way that will help you generalize what you have learned.

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
Interpolation: automatic order switching interpolation
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

Of course, the solution type is not interesting in and of itself. We want to understand the solution! The documentation page which explains in detail the functions for analyzing the solution is the [Solution Handling](https://docs.sciml.ai/dev/basics/solution/) page. Here we will describe some of the basics. You can plot the solution using the plot recipe provided by [Plots.jl](http://docs.juliaplots.org/dev/):

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
Interpolation: automatic order switching interpolation
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

DifferentialEquations.jl has a common set of solver controls among its algorithms which can be found [at the Common Solver Options](https://docs.sciml.ai/dev/basics/common_solver_opts/) page. We will detail some of the most widely used options.

The most useful options are the tolerances `abstol` and `reltol`. These tell the internal adaptive time stepping engine how precise of a solution you want. Generally, `reltol` is the relative accuracy while `abstol` is the accuracy when `u` is near zero. These tolerances are local tolerances and thus are not global guarantees. However, a good rule of thumb is that the total solution accuracy is 1-2 digits less than the relative tolerances. Thus for the defaults `abstol=1e-6` and `reltol=1e-3`, you can expect a global accuracy of about 1-2 digits. If we want to get around 6 digits of accuracy, we can use the commands:

````julia

sol = solve(prob,abstol=1e-8,reltol=1e-8)
````


````
retcode: Success
Interpolation: automatic order switching interpolation
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

More advanced saving behaviors, such as saving functionals of the solution, are handled via the `SavingCallback` in the [Callback Library](https://docs.sciml.ai/dev/features/callback_library/#saving_callback-1) which will be addressed later in the tutorial.

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
Interpolation: automatic order switching interpolation
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
Interpolation: automatic order switching interpolation
t: 10-element Array{Float64,1}:
 0.0
 0.04297991782675459
 0.11895240940449249
 0.20804668088620082
 0.31572426897446076
 0.44241711288613594
 0.5849882423692904
 0.7423235141702504
 0.9131601569809136
 1.0
u: 10-element Array{Array{Float64,2},1}:
 [0.16644874626312411 0.7648556689386914; 0.7574214464300553 0.249424456632
42974; 0.8384677441652046 0.7903670489396488; 0.4184056814230097 0.35676734
57715661]
 [0.07175649710715612 0.6922706532497522; 0.7951327566627149 0.410617506396
8266; 0.8376325440746346 0.6851057024982075; 0.5093015838043589 0.611159786
2810678]
 [-0.14995119929672543 0.4212263220844022; 0.7925911119471691 0.51124018260
49353; 0.89168892770078 0.5741797952908354; 0.6367806766436428 1.0384233546
432857]
 [-0.4825974280989417 -0.1265681508498277; 0.7087852137321302 0.35892861655
997654; 1.0634219968651892 0.6270134718543703; 0.718925639280158 1.47276275
98042963]
 [-0.947123684241147 -1.0909952730580932; 0.5659541789318394 -0.11716176880
380064; 1.4474937852601266 1.0586407826937836; 0.6981474427980019 1.8424403
394066964]
 [-1.485463136044281 -2.544029276891878; 0.4865346418407177 -0.843705993169
7983; 2.1436547203895264 2.213805155694563; 0.474246704743588 1.96031688433
40755]
 [-1.8950760540352147 -4.320776694728493; 0.7477904610283553 -1.39337474819
3556; 3.1559526711772996 4.430673039275123; -0.06338096208105759 1.54048408
46754344]
 [-1.7991960988478894 -5.896202282427029; 1.7623051465675168 -0.88598461160
54957; 4.2785997598146865 7.838497459181548; -1.000040017443705 0.238352961
1513815]
 [-0.6238206223751326 -6.172487945714782; 3.951579234225358 1.9997936905090
024; 4.901813525025712 11.923916171717064; -2.315624269273623 -2.2600758503
64929]
 [0.5284337383988529 -5.34356411220458; 5.483763586694698 4.670235909860765
; 4.698298312846905 13.680115397976639; -3.0220912203180212 -3.926768838506
1534]
````





There is no real difference from what we did before, but now in this case `u0` is a `4x2` matrix. Because of that, the solution at each time point is matrix:

````julia

sol[3]
````


````
4×2 Array{Float64,2}:
 -0.149951  0.421226
  0.792591  0.51124
  0.891689  0.57418
  0.636781  1.03842
````





In DifferentialEquations.jl, you can use any type that defines `+`, `-`, `*`, `/`, and has an appropriate `norm`. For example, if we want arbitrary precision floating point numbers, we can change the input to be a matrix of `BigFloat`:

````julia

big_u0 = big.(u0)
````


````
4×2 Array{BigFloat,2}:
 0.166449  0.764856
 0.757421  0.249424
 0.838468  0.790367
 0.418406  0.356767
````





and we can solve the `ODEProblem` with arbitrary precision numbers by using that initial condition:

````julia

prob = ODEProblem(f,big_u0,tspan)
sol = solve(prob)
````


````
retcode: Success
Interpolation: automatic order switching interpolation
t: 6-element Array{Float64,1}:
 0.0
 0.09966372028780551
 0.3430294439172079
 0.6404424417661176
 0.9602673462329638
 1.0
u: 6-element Array{Array{BigFloat,2},1}:
 [0.16644874626312411436401816899888217449188232421875 0.764855668938691390
934536684653721749782562255859375; 0.75742144643005526738477328763110563158
98895263671875 0.2494244566324297363024697915534488856792449951171875; 0.83
846774416520464257018829812295734882354736328125 0.790367048939648819683156
9439382292330265045166015625; 0.4184056814230097121054541275952942669391632
080078125 0.356767345771566102285987653885968029499053955078125]
 [-0.0875706938487612410420016981558766485104539240018417440384239045229021
3815951498 0.50737623740836119225566556752854775603554511468411661872988097
71286628961249256; 0.800583917386213945122127974713613725871538564058745448
9875530770745609066080593 0.50721047569577046851612130484666671502416369802
96649521676657783022179692240274; 0.870531778220792334981714576422969126485
2403616565249166742072121906959959091394 0.59105356350090042390541268736732
23961745480343817729039253091287960618497534132; 0.608916461996469615887660
1528014966321504333878955047821579572918348250560415089 0.93371742564251247
79287687324282747361067050861860266829166018689300067523573028]
 [-1.0680872643895192608967856389361936062604201135032587867220503261749095
347613 -1.38039593755857488620517070856575417622384353701667932374377441009
2382911902784; 0.5348075723758538927457890802943566452525449916990790014087
586595405255736311245 -0.26939986162458016194411274341289010408660949392457
52620312634262392621413861391; 1.576241165640527184272407315741870238390481
660126733898823939256816391608961134 1.244599569950719716433272074881346260
195920151952047003230584661970565606602582; 0.66911824071279589133588366998
22886784312931633379837556255450050528638155852949 1.9005927864967098363749
76113939758665831393694491007318899795082148532284760549]
 [-1.9449642889737408085440161875151162663670360123696807024065831151229989
97939808 -4.959352545643675322430007166272162655308426121277807546219810638
922834439618677; 1.00464884854757900524477465188403065954142379524227781993
000529261878664381828 -1.39306225750866632584658379564153810442537041383970
4126563204925650526500281618; 3.5719008437274272041054028305422424858636556
60414206629667899786475765293923282 5.5372646955157232107775710958600400084
99602142404417641948211015074976877728034; -0.35482865341099980038565593057
16858651854584268932023285412931735173050499203127 1.1887329037036185416619
07992439578020425075975087545453335895888204195494174665]
 [-0.0493512820139255354776049321250292055572674112441557901169237526285244
4754681337 -5.8218981049343936550272056777861224685827924330076197523774729
93243374198867612; 4.752043366920882341256435540424732447914660509051083765
093277711594700856490336 3.340437882699244542132809431835009015994101675563
918994721621667092518370923899; 4.84893575176579252963001998579433600755572
0185162817296837157489101838531010464 12.9306939555687068072614748616322490
7933337001943082874057505141631451134545785; -2.700764647635840882401398036
458257028471889685398534772623582840182362322295917 -3.13530696393598335213
8979614318112248961695624540587439146489768420200941619209]
 [0.52844597052349079151538331621448153033690227281798084771312347755742086
34349984 -5.343561430923411916663579551608836948519276500484634365369391490
015126607191726; 5.48378616614306503933516814368036173716081298410630098351
0357479090418801709764 4.67027132471155859206030867458042408533511781742260
9095316965127635921892131967; 4.6983007492415210057724752807122958927871903
89956065317926436605055643810036407 13.680147600662823761073911195613376381
89135002982548515676927557538736649248298; -3.02210133688816024363585133831
2441466986283042393158656761269299324875964156094 -3.9267888772087156425233
18545199414674998124445370458175062526598216862645434237]
````



````julia

sol[1,3]
````


````
-1.068087264389519260896785638936193606260420113503258786722050326174909534
7613
````





To really make use of this, we would want to change `abstol` and `reltol` to be small! Notice that the type for "time" is different than the type for the dependent variables, and this can be used to optimize the algorithm via keeping multiple precisions. We can convert time to be arbitrary precision as well by defining our time span with `BigFloat` variables:

````julia

prob = ODEProblem(f,big_u0,big.(tspan))
sol = solve(prob)
````


````
retcode: Success
Interpolation: automatic order switching interpolation
t: 6-element Array{BigFloat,1}:
 0.0
 0.099663720287805519850364474986984514692436826154524436670624328377388155
82142945
 0.343029443917207968135337710647404410910375827827586442269037916172444070
4967693
 0.640442441766117732304786095599486692725838678779676482528375313036180311
5985089
 0.960267346232963927519654549117369480718664969067561693257288212725417817
6310427
 1.0
u: 6-element Array{Array{BigFloat,2},1}:
 [0.16644874626312411436401816899888217449188232421875 0.764855668938691390
934536684653721749782562255859375; 0.75742144643005526738477328763110563158
98895263671875 0.2494244566324297363024697915534488856792449951171875; 0.83
846774416520464257018829812295734882354736328125 0.790367048939648819683156
9439382292330265045166015625; 0.4184056814230097121054541275952942669391632
080078125 0.356767345771566102285987653885968029499053955078125]
 [-0.0875706938487612629901597654554148088713632176501230920399944261274521
0682965103 0.50737623740836116309651712520917255130187766904330917456199197
08026412621003387; 0.800583917386213943047435444256462011690426624406598353
9261863913151344228503384 0.50721047569577047256740463526479614604542272778
44457599312625898137192457267903; 0.870531778220792341703179477584125863289
6216800962729188430678038062905829889684 0.59105356350090041622682299038097
83600177142163866652818383164209277100818861217; 0.608916461996469626600470
3159331512709421363340561432053399503166390898478902564 0.93371742564251251
65093887647374719986139780400888059061692310508387643087116769]
 [-1.0680872643895193942100393698834892771560388228407548640642222032658690
10698407 -1.380395937558575214932360739446347002782972858037762019902301505
881349434650393; 0.53480757237585386120152029645978933788473151516491020946
25204127290467036816082 -0.269399861624580334296335872304840415837000736234
4196266012314573647436225309273; 1.5762411656405273335272954349275933698445
58391460520683145873518806485977421887 1.2445995699507199406169138772897780
79892344290281174906681777961031096125389886; 0.669118240712795853573696302
7622884983529817491334111687903520161053232590034331 1.90059278649670989158
2811701475677926757636807923689741730731427730892406605814]
 [-1.9449642889737408191018037105216714501507621389728494928596322163410319
36290607 -4.959352545643675996304319820668920931598665819050693483167232615
481668420703067; 1.00464884854757934906877324312016940967420153219271803655
4283261380679272051385 -1.3930622575086662311857123353514234679785547261445
08890428345343339358418102478; 3.571900843727427663018376959440160342294320
804453563049544296181796343050931521 5.537264695515724510324115136712382977
931321412358295210407785314388526582700816; -0.3548286534110001498877160420
380081892482367916622206724446400098442747903885741 1.188732903703618086155
566625448380690600911804680071845593908721045994219534006]
 [-0.0493512820139246846638157856131727246506121083259107073194114566567024
3669480624 -5.8218981049343930318543770785883208878601897805242216779594758
18336434657345035; 4.752043366920883466647725033374185579123804324469736236
424627627120721104306346 3.340437882699246512596758027205413491663300892933
56139349717882752238286250962; 4.848935751765792371326274432980422446003022
630837626938087642310854440309646044 12.93069395556870808162148960093429257
192974243829752960884337064506848753878475; -2.7007646476358413981118748735
29572876238745489407714572580800445593801784255244 -3.135306963935984574810
417043610997317816493537815287203544131560716906119591985]
 [0.52844597052349035744725627049129930609389714689979840719744333009111819
88123024 -5.343561430923412313301338750243782503332517719001405244013649171
972710521801053; 5.48378616614306451182048494803244509884868586167632095592
4482380876896096954667 4.67027132471155759879459776895328836102555904113257
7092760104250605613228887431; 4.6983007492415211483218860026561758773029009
18729143561226682901460937202392995 13.680147600662823276809615386835137281
38683357233436296204227727820226742384662; -3.02210133688816002172875924079
1134290178607739838387008686883323075755601029812 -3.9267888772087150741329
32363470698332503176668453472690492731271153530920982005]
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
Interpolation: automatic order switching interpolation
t: 11-element Array{Float64,1}:
 0.0
 0.03853083417381429
 0.09826500352088789
 0.1684337899446925
 0.24569284333590585
 0.33515674170740595
 0.4583576152442672
 0.6004143088434256
 0.7563404507387839
 0.9228596575369267
 1.0
u: 11-element Array{StaticArrays.SArray{Tuple{4,2},Float64,2,8},1}:
 [0.8465039773841412 0.12487166414587691; 0.8772247094064545 0.497955640492
99057; 0.22884474240226926 0.06307747595861102; 0.4567288576756148 0.866106
5467944384]
 [0.7739965120707026 -0.04880693327064392; 0.8985029791127979 0.37669685515
57989; 0.12429267329343546 0.09196362746020416; 0.6204474296115354 0.950867
7100254495]
 [0.5934233204670447 -0.3613596209850293; 0.8396541330153405 0.149785203565
87963; 0.0036427491064897066 0.20007822970610006; 0.8619760077530711 1.0559
784130643715]
 [0.2763890357067625 -0.7865759692215154; 0.6408039555377326 -0.15464661641
391972; -0.05139672572401409 0.43691787963817835; 1.1167421936218247 1.1306
849432446628]
 [-0.19733658863852283 -1.3084591286568084; 0.28859649361354983 -0.49708437
63065824; 0.028684558830358253 0.8477525209551958; 1.3453636331141348 1.139
9195600814844]
 [-0.8895079472812393 -1.9458323148873498; -0.2350224242915035 -0.834894818
4935793; 0.34791320695438954 1.5280707377228873; 1.5179155740116161 1.03760
92128248002]
 [-2.023887218176504 -2.777167411599006; -0.9965987142591297 -1.03584153684
99858; 1.2507138891278475 2.804832347054119; 1.5408833664301784 0.665033234
0910021]
 [-3.400715580102548 -3.435426615630712; -1.585000719753579 -0.596280942808
5849; 2.993189947039186 4.634756275083157; 1.1697794469961615 -0.1393145658
902617]
 [-4.595818658536972 -3.374105377261275; -1.3055210396451593 1.093261144524
8255; 5.630590264344536 6.692498195370133; 0.1656488252670658 -1.4943400903
90621]
 [-4.834281495887991 -1.826952431254577; 0.756470293843353 4.57708439097409
; 8.739613490549637 8.125131590755691; -1.6552043912176375 -3.3802343400003
98]
 [-4.3478033134501235 -0.4307208255563493; 2.48777093171584 6.7673715342326
3; 9.99926127415837 8.188358126113654; -2.7418017609232024 -4.3286047721021
84]
````



````julia

sol[3]
````


````
4×2 StaticArrays.SArray{Tuple{4,2},Float64,2,8} with indices SOneTo(4)×SOne
To(2):
 0.593423    -0.36136
 0.839654     0.149785
 0.00364275   0.200078
 0.861976     1.05598
````





## Conclusion

These are the basic controls in DifferentialEquations.jl. All equations are defined via a problem type, and the `solve` command is used with an algorithm choice (or the default) to get a solution. Every solution acts the same, like an array `sol[i]` with `sol.t[i]`, and also like a continuous function `sol(t)` with a nice plot command `plot(sol)`. The Common Solver Options can be used to control the solver for any equation type. Lastly, the types used in the numerical solving are determined by the input types, and this can be used to solve with arbitrary precision and add additional optimizations (this can be used to solve via GPUs for example!). While this was shown on ODEs, these techniques generalize to other types of equations as well.


## Appendix

 This tutorial is part of the SciMLTutorials.jl repository, found at: <https://github.com/SciML/SciMLTutorials.jl>.
 For more information on doing scientific machine learning (SciML) with open source software, check out <https://sciml.ai/>.

To locally run this tutorial, do the following commands:
```
using SciMLTutorials
SciMLTutorials.weave_file("introduction","01-ode_introduction.jmd")
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
[91a5bcdd-55d7-5caf-9e0b-520d859cae80] Plots 1.6.3
[90137ffa-7385-5640-81b9-e52037218182] StaticArrays 0.12.4
[c3572dad-4567-51f8-b174-8c6c989267f4] Sundials 4.3.0
[37e2e46d-f89d-539d-b4ee-838fcccc9c8e] LinearAlgebra
```
