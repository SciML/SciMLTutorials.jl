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
 0.051702071063330475
 0.13562999422444041
 0.23037554474810376
 0.3426983013613073
 0.47176712202063925
 0.6149940094645234
 0.7743407691370402
 0.9475780987555275
 1.0
u: 10-element Array{Array{Float64,2},1}:
 [0.7224595586568736 0.24874577172817847; 0.5942403386924497 0.898827791378
2818; 0.649495010203508 0.5527041618294715; 0.11316625617514475 0.392222860
6003171]
 [0.7032596978005325 0.14675221643974926; 0.7581025440265211 0.891797644539
154; 0.5123816816024496 0.5338172965304608; 0.32172852001866103 0.474372304
18566056]
 [0.5532545184843075 -0.07090838325944357; 0.8533538446961457 0.81139425452
30367; 0.33807143174410337 0.563677012549482; 0.6449789878234523 0.57210509
40797825]
 [0.2072110953866354 -0.3764846100934285; 0.7288263285319521 0.654145579854
575; 0.26486357638022745 0.7038101630475819; 0.9635794826135089 0.617612177
7569383]
 [-0.4270525821375436 -0.7790021312307465; 0.33440960891894883 0.4518175021
615276; 0.43064528231091403 1.0303547949712062; 1.2355358026772887 0.564357
9967875104]
 [-1.3871904507622688 -1.2094926918989968; -0.2783593858431989 0.3298704055
553676; 1.059698822660427 1.6077684404603216; 1.3405450444008722 0.33823263
20677363]
 [-2.5718602674645084 -1.4938959714599345; -0.8332194726232629 0.5245758837
765426; 2.373809264219824 2.4137716471680077; 1.1073669827138528 -0.1331476
6741844347]
 [-3.6823033241799488 -1.327218869505374; -0.7685836050181956 1.37163172235
86019; 4.512271651980489 3.267173162648726; 0.3133804601951733 -0.910751823
3601639]
 [-4.0142028619691725 -0.24773054065126598; 0.8044462675615107 3.1825209520
952367; 7.183017378159402 3.6295324180704194; -1.2500570287843304 -1.951519
2920030515]
 [-3.813885057351085 0.30721078510841193; 1.681269283693354 3.8997616762895
784; 7.925469374380137 3.514504534721358; -1.8593116515826504 -2.2735375462
877445]
````





There is no real difference from what we did before, but now in this case `u0` is a `4x2` matrix. Because of that, the solution at each time point is matrix:

````julia

sol[3]
````


````
4×2 Array{Float64,2}:
 0.553255  -0.0709084
 0.853354   0.811394
 0.338071   0.563677
 0.644979   0.572105
````





In DifferentialEquations.jl, you can use any type that defines `+`, `-`, `*`, `/`, and has an appropriate `norm`. For example, if we want arbitrary precision floating point numbers, we can change the input to be a matrix of `BigFloat`:

````julia

big_u0 = big.(u0)
````


````
4×2 Array{BigFloat,2}:
 0.72246   0.248746
 0.59424   0.898828
 0.649495  0.552704
 0.113166  0.392223
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
 0.07547982651819603
 0.3051645310535731
 0.6071222946756546
 0.9253561445439225
 1.0
u: 6-element Array{Array{BigFloat,2},1}:
 [0.7224595586568736127475176544976420700550079345703125 0.2487457717281784
663754251596401445567607879638671875; 0.59424033869244974326306873990688472
9862213134765625 0.89882779137828183735337006510235369205474853515625; 0.64
94950102035079719797749930876307189464569091796875 0.5527041618294714631076
65840652771294116973876953125; 0.113166256175144752660344238393008708953857
421875 0.392222860600317080326249197241850197315216064453125]
 [0.67581734147909378580896665477205728428628440533112412349928154999418656
38221693 0.0912267369714313614816760407409668048679183191268161199057856107
1734852406713231; 0.8061880528312891426061770976393405332150407275244938781
412269184548194660623026 0.876795458241779563829083764218313330814859967399
0162774226459217474824749523651; 0.4554773794280405331550788496293619618164
34223166858753020147768198155249910871 0.5341319978064587459689221614144620
248199282212197643810597337220631379900407861; 0.41587268661909300224805822
44376659888099690021832219181827760769716794495331162 0.5068783200262910615
069327117818815591992177917131834405810542525321495998382656]
 [-0.1904393428781420394984313564322328726867419542806385284800762757078884
408459078 -0.64317965214501450478187911394991422757452750260578674694632923
80223753241108066; 0.489227455705703264412448926535167933132391428383291683
5658692247756674501252196 0.51553016383526395427397919075294791021777091646
58078236457092695445210909231248; 0.339560106669062212635131737376296938160
6236862441834839988873062182718562195938 0.90175834096187655273794561548455
9087342172498291652554150578280797501534093711; 1.1604226515246721511293018
83929390766021989663855287982832174237770476464448024 0.5962066464031061455
67431437436465089847277221978415922355069510356458085035063]
 [-2.5078098966852356758518696434372920075101643137476733866069249660135527
73927736 -1.486775273523266590568102314652692117597174642307310213005611752
278343112166893; -0.8131136268391734730313825745740209981392775140640277771
364152539440582868996759 0.501537012172555023537087933737720227718398416858
1626131563652322909610424209927; 2.2850279611205931562498663978728806331028
43134192675420504451443525281886569945 2.3677664486312874073586597244041110
90904312648346189342174186888921795900568758; 1.131264759357704092303731085
639043542032136550356752388750019602077892472486795 -0.10128404443417442684
84001263087455777964104462768471428568900047561020369016816]
 [-4.0499944397859277808980296608761644684741821731337253430910453810936337
56764685 -0.449118939912150863334884047987430574540131956884402680246529070
2871798316396134; 0.4927788820906759288226765489518973965131244777202811793
800915989199835806731487 2.900011625494969984624674214709393481956523085643
859277692177382526852617530636; 6.84951966736828958305660857574540721917575
9501599414869200824807354772824641791 3.64025490678090610900958342070961421
1927285668077262798025176774209344524079932; -1.009708088419327135746545424
918262242320766465703765858633376182464612724564408 -1.81363693019522104143
1114868541210879824719621639633310366797253482026761746354]
 [-3.8138870367582791223318544666819391289939337699419814332851718293990503
19351548 0.3072203604937003709979962049538138835143880357240232385194604170
981873859859486; 1.68128831901061408181785017254212287253718445367179820451
9553703012000535950065 3.89977956564943673985282958711861300093959042715347
6790334461363069106913763622; 7.9254923680345958950438870819418045146863729
48152672009201293938384630806738748 3.5145069769139993062715999559101538549
11376137416042014941312951514498544409901; -1.85932367123712800927823716737
6847355841945396924778892957140187495021832524299 -2.2735459674734173979031
03422847233763521292904689382281851304603941004153488952]
````



````julia

sol[1,3]
````


````
-0.190439342878142039498431356432232872686741954280638528480076275707888440
8459078
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
 0.075479826518196020944626426131893343501211414735435567331157261844111936
57668153
 0.305164531053573089087551485038087012529163394675542355930797213294530696
7451567
 0.607122294675654530192813149991753325572084922706560607086345250167758287
6733069
 0.925356144543922327688409180801795514694493132258846925489902051566095127
3237512
 1.0
u: 6-element Array{Array{BigFloat,2},1}:
 [0.7224595586568736127475176544976420700550079345703125 0.2487457717281784
663754251596401445567607879638671875; 0.59424033869244974326306873990688472
9862213134765625 0.89882779137828183735337006510235369205474853515625; 0.64
94950102035079719797749930876307189464569091796875 0.5527041618294714631076
65840652771294116973876953125; 0.113166256175144752660344238393008708953857
421875 0.392222860600317080326249197241850197315216064453125]
 [0.67581734147909379517566817849591675708849899686023694323016403401481494
67898453 0.0912267369714313777863746774036237766099483744065580220057955477
9344183819622455; 0.8061880528312891314934047208947664830192018099718395702
275105702739711729973499 0.876795458241779568986379341879191159274790715566
4269520573387026569201110206426; 0.4554773794280405484202433421330690372343
942824332329777554907834579999537230324 0.534131997806458745021462383857467
9392106863761554457282342333523612838637923921; 0.4158726866190929760522970
627896787503257265813695400856208164782201514974575288 0.506878320026291052
888370601951529086406451505008523813791310055797154902868144]
 [-0.1904393428781418487271876148148369972489231168488428958167380468058150
436036895 -0.64317965214501438940607108800743787277278726965214009210458808
82528898241286688; 0.489227455705703387397612012679355077537878978988275605
8272465597537170987452489 0.51553016383526401111059273903025584979213368553
92841326120498165308590404360151; 0.339560106669062151443065443103792567114
2631234296319318532924481451323333554582 0.90175834096187645185609236999872
61364255111198258652076780721531982470771112375; 1.160422651524672080146259
904436609386971076489697885320085385803534952151837583 0.596206646403106166
4135065883188447185176676097954622647131987897154552282804992]
 [-2.5078098966852350365589256084152273521059980766447118095935877272502794
04508604 -1.486775273523266513801359990915965798719888226043558002705334430
460515301981105; -0.8131136268391732648421695412953179785738515383652513849
941921000276733493059297 0.501537012172554802345423424363577507889093215184
4102795197559025125267343341279; 2.2850279611205922821693173466388290329515
9895832797625869382807842767376076557 2.36776644863128694960163329633983008
9745415286042103797729202194864034620285799; 1.1312647593577043232231673780
60390881609420057778161499303142477505139785662504 -0.101284044434174113215
526705482222371929142680633781246482347675893417719505296]
 [-4.0499944397859279210876351253901314377141613676462771805834072011727114
29777592 -0.449118939912152073396588204380885439312311768240657004546720852
6530939368798249; 0.4927788820906740697772502631404411633228591002210530604
827817496381959115687178 2.900011625494968242981547907251070265099212386792
260537646902140220754052093768; 6.84951966736828745044015514671227433487055
9591705856085865830681618375121650383 3.64025490678090611141903638256522238
4409663153975582587019708069955592003644037; -1.009708088419325652401848778
200964463535022435445647237338512045158541976603757 -1.81363693019522017014
5818208692052521929629447482832536664156374894009269355727]
 [-3.8138870367582791984201767991166345158945860171600853396858642608875190
23587159 0.3072203604937002089754871662680334338468748929006840425221761127
825843109482911; 1.68128831901061382283309336057126703947441886176033469891
4599151028845911488442 3.89977956564943654128968417784382956479809727909174
1459043764816880185655595219; 7.9254923680345957091341536079467119067347269
73331483682883264079348020135136683 3.5145069769139993548776858875459687951
03389022710718816538727190620729367375269; -1.85932367123712784053827101424
0190859332342319173112883202426051757648495325054 -2.2735459674734173138720
23276730215516186821919658606473999319298826086187264111]
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
t: 10-element Array{Float64,1}:
 0.0
 0.04988475413263506
 0.12496075270086432
 0.20952999074904327
 0.31143311490330344
 0.4330684020433413
 0.5731646337595461
 0.7299987573564539
 0.8905974335931316
 1.0
u: 10-element Array{StaticArrays.SArray{Tuple{4,2},Float64,2,8},1}:
 [0.21564908876620326 0.769633786771917; 0.3072951153541328 0.9407682161751
918; 0.16909664737853913 0.45679515841168117; 0.32674323699202956 0.3962499
0947154815]
 [0.13214944516816202 0.6822686221890325; 0.28878773683617914 0.99372235681
14463; 0.1524304317728905 0.3358035631550183; 0.41188724757969186 0.5948947
352267456]
 [-0.040136253366522184 0.44941630707686775; 0.20760001722174232 0.93512147
11666537; 0.17246951638371333 0.21806882289440882; 0.522148889295442 0.8703
924267576695]
 [-0.2950666935501699 0.04629270231216687; 0.05789476418731393 0.6964941729
450265; 0.27572443341430203 0.2151908032155574; 0.6115191367455922 1.129449
265367407]
 [-0.6695454507757463 -0.6154127888110585; -0.16184791894496806 0.234683996
7734522; 0.5354445858150272 0.45210245360121715; 0.654157557165632 1.338094
3960289737]
 [-1.1653418834839968 -1.5890248661190043; -0.39082003989237646 -0.40856968
013802386; 1.0579145455477859 1.1505580020541686; 0.5859953665167703 1.3863
469857905024]
 [-1.6944927452613943 -2.7884334503447237; -0.4444175963405534 -0.961730734
327631; 1.9302031654750036 2.557252859150692; 0.31108058827598384 1.0923040
435409612]
 [-2.0239775952927754 -3.874253812408612; 0.013005121850123813 -0.828797323
8695163; 3.1245254846749577 4.784921984660782; -0.27909967100606137 0.22715
70656934636]
 [-1.7695189917591365 -4.14920549747011; 1.300400825719739 0.72486962614538
09; 4.270278119370374 7.347883018362731; -1.1856197560186312 -1.29133053053
72706]
 [-1.0794377722139354 -3.505614680035695; 2.723187104220299 2.8455236559221
264; 4.7413164070116105 8.838706185315274; -1.9364260835270608 -2.659306294
7413317]
````



````julia

sol[3]
````


````
4×2 StaticArrays.SArray{Tuple{4,2},Float64,2,8} with indices SOneTo(4)×SOne
To(2):
 -0.0401363  0.449416
  0.2076     0.935121
  0.17247    0.218069
  0.522149   0.870392
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
[91a5bcdd-55d7-5caf-9e0b-520d859cae80] Plots 1.6.6
[90137ffa-7385-5640-81b9-e52037218182] StaticArrays 0.12.4
[c3572dad-4567-51f8-b174-8c6c989267f4] Sundials 4.3.0
[37e2e46d-f89d-539d-b4ee-838fcccc9c8e] LinearAlgebra
```
