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
 0.04018510695223035
 0.11203472865898827
 0.19549824155043496
 0.29790008667570517
 0.4228354632291847
 0.5687866645125571
 0.729379105150639
 0.9002384713134701
 1.0
u: 10-element Array{Array{Float64,2},1}:
 [0.7637548238660845 0.015803866809346978; 0.3474524613246084 0.37343280589
831185; 0.5249741961172385 0.5965552620587542; 0.6351532866089222 0.0636681
1036365938]
 [0.6402175098666107 0.0006037801944166459; 0.41601993251156627 0.429164140
30415967; 0.4418074384826368 0.5982663265278898; 0.8757665695906751 0.09041
848621022988]
 [0.2862594757753755 -0.04020129399261979; 0.37623997379820256 0.5061412983
900834; 0.38224751510469845 0.6114691667062325; 1.2785601486874816 0.126695
7048707651]
 [-0.333936569561097 -0.10405481326601594; 0.10151992274017108 0.5671500340
430028; 0.5080872193917743 0.6467545787533748; 1.6760665824784855 0.1465985
5007833558]
 [-1.3705873843789984 -0.1926258390944256; -0.4733912111137165 0.6198228110
796913; 1.0380809960525847 0.7221276332374202; 2.0078384004987537 0.1325966
1414373825]
 [-2.9368209648091046 -0.28410534810870725; -1.2819041969257159 0.691741337
731596; 2.3641153311508862 0.855133804001645; 2.0839418186986336 0.05157366
857847365]
 [-4.884070261324253 -0.3045315821676079; -1.8438324175551326 0.85566752670
63094; 4.926651770014615 1.0285480661180786; 1.562457201740442 -0.133292817
8720064]
 [-6.532010306175154 -0.12157235068007999; -1.1066229366117604 1.2141676174
186506; 8.787568370885724 1.1355810046432984; 0.05248097613732483 -0.429429
00111682425]
 [-6.6353434696133835 0.41836916987764117; 2.3167982303665955 1.82938329738
57505; 13.20978405135482 0.9528893310533812; -2.717831553464249 -0.78796784
45840644]
 [-5.374418625048119 0.9249705080340971; 5.909527417590797 2.26648785453837
04; 15.271296611398837 0.6014601259756733; -4.83088087362685 -0.97271848081
71725]
````





There is no real difference from what we did before, but now in this case `u0` is a `4x2` matrix. Because of that, the solution at each time point is matrix:

````julia
sol[3]
````


````
4×2 Array{Float64,2}:
 0.286259  -0.0402013
 0.37624    0.506141
 0.382248   0.611469
 1.27856    0.126696
````





In DifferentialEquations.jl, you can use any type that defines `+`, `-`, `*`, `/`, and has an appropriate `norm`. For example, if we want arbitrary precision floating point numbers, we can change the input to be a matrix of `BigFloat`:

````julia
big_u0 = big.(u0)
````


````
4×2 Array{BigFloat,2}:
 0.763755  0.0158039
 0.347452  0.373433
 0.524974  0.596555
 0.635153  0.0636681
````





and we can solve the `ODEProblem` with arbitrary precision numbers by using that initial condition:

````julia
prob = ODEProblem(f,big_u0,tspan)
sol = solve(prob)
````


````
retcode: Success
Interpolation: Automatic order switching interpolation
t: 6-element Array{Float64,1}:
 0.0
 0.11066848054444339
 0.34969130745307003
 0.6507768497178879
 0.971271619312092
 1.0
u: 6-element Array{Array{BigFloat,2},1}:
 [0.7637548238660845090208795227226801216602325439453125 0.0158038668093469
7807288102922029793262481689453125; 0.3474524613246083859507962188217788934
7076416015625 0.37343280589831184812510400661267340183258056640625; 0.52497
41961172385007472485085600055754184722900390625 0.5965552620587541987617896
66605181992053985595703125; 0.635153286608922229561358108185231685638427734
375 0.0636681103636593803685173043049871921539306640625]
 [0.29457217844753972330140792065331594335045078250348221084840344006681996
36013995 -0.039283332922574536292861102069405979455975852780635151012477307
91313533154599686; 0.378821620917826329630878685949135228948918468625509283
9985488125281907251352488 0.50491744373914835297948695844662139060070202938
21290574153081912499904303620086; 0.382092777861873759716395691605006861620
2623255318835747173844714976113493555808 0.61107924689751476454293377498209
62238679820450917009856990346872227630909770667; 1.271337872187865771412738
984550766760302951988345904091214481277974975569039645 0.126160049694481879
734966175540987662335114857270324644603094015379590786959125]
 [-1.9892926282818621900197912256130028748773076217535806825860776340085075
7369819 -0.2351111207844671210211505224171749474368397131790208752209810956
846134655607402; -0.8133712259133945584805087295674760513317379625985085575
115213292710851175471148 0.645802717579950494643509518655425191537427126224
3825181914192099180845988940512; 1.4915723320317631671842021633356929151593
64968638731771257016633304947406788783 0.7727846435558810491034147158424620
989153631309331168539731719540580874912342454; 2.08952603477208884626517129
3161910351816051512848895274730321996406852790161419 0.10779206885146547658
06990618400927263515616631804273476625316639950905334914731]
 [-5.8436678475230153350057197994238684200662924413580915017424480179935398
56803902 -0.244621598701266575261780725884407600658135059943225268057142727
0232645890027573; -1.707192889125725729748231339406664265479135727660948471
152989385581618194225677 1.011739730724958513208038280166042812406417717434
136780577270347345730775765552; 6.79451521497870498424179888848328465440692
560048245549399469326702090439283259 1.103749039986005461291303388277713990
344595504132308186701372766824886185862815; 0.92290019864425632829091723206
96043072999834114360721198367816014686019517105582 -0.274841299713644495789
3351213289128247736164627121582177047190927484441042263108]
 [-5.8605010633542192002610608231047198895904290914052445384610152336878782
6628228 0.76447118867664096396705092409965622742549855429590243177660179070
98462275396826; 4.747241385252438316383061628863536296997881452742151793823
443462735151882621724 2.137859881705940420934978430112733223191928575611197
235184582675582087470871472; 14.7543422896124773182365625932780330909063019
8563755140357872728085152415882592 0.72566897389483529103747384956695807792
81212601662062972286121331628784140912331; -4.19139067862137130664772459540
5946501433321045338783885186878924171298648796172 -0.9238562474921242976030
09595745951179317391997216675117854015535477431442433709]
 [-5.3744128793052812184455230021782629546357551704322158854726675104467278
66417488 0.9249765388304762539468854311212646297128596498894704290365151779
679263913522488; 5.90957132286666258627380774829416831594863651897120967107
0572522504336108686603 2.26649528024961185300600861691960893234475462925129
7955652763649051904326019956; 15.271333102751646217641437324766575399053905
72456413350361989120100755939894328 0.6014570411302269580019589290883742856
703714351035806483402069177125374063041795; -4.8309053298570960090548614911
05422454888996112667949054566435446231737319341218 -0.972721147184919753984
0645381251440841349220377394155925147718825670339577242146]
````



````julia
sol[1,3]
````


````
-1.989292628281862190019791225613002874877307621753580682586077634008507573
69819
````





To really make use of this, we would want to change `abstol` and `reltol` to be small! Notice that the type for "time" is different than the type for the dependent variables, and this can be used to optimize the algorithm via keeping multiple precisions. We can convert time to be arbitrary precision as well by defining our time span with `BigFloat` variables:

````julia
prob = ODEProblem(f,big_u0,big.(tspan))
sol = solve(prob)
````


````
retcode: Success
Interpolation: Automatic order switching interpolation
t: 6-element Array{BigFloat,1}:
 0.0
 0.110668480544443385144722681623580989434451154580994219717861493916260649
8027799
 0.349691307453069995032215488856673647275777284160775862405509510342263101
4312538
 0.650776849717887776346263402349261397623796534801642043607649870053772486
2142495
 0.971271619312091948957930567740879115714675683770279575094010499569202415
9077896
 1.0
u: 6-element Array{Array{BigFloat,2},1}:
 [0.7637548238660845090208795227226801216602325439453125 0.0158038668093469
7807288102922029793262481689453125; 0.3474524613246083859507962188217788934
7076416015625 0.37343280589831184812510400661267340183258056640625; 0.52497
41961172385007472485085600055754184722900390625 0.5965552620587541987617896
66605181992053985595703125; 0.635153286608922229561358108185231685638427734
375 0.0636681103636593803685173043049871921539306640625]
 [0.29457217844753976956794814101607456694510434469852826400545590074284410
18512334 -0.039283332922574531178732154026607272871273293195963562158898804
15259782944801773; 0.378821620917826343864678140096278331421648564292341421
2352413183448967974505587 0.50491744373914834611925316914238591319220641843
32191765745717723101353684876187; 0.382092777861873759006237424983547623700
7127743430653265811901955884564812842144 0.61107924689751476238081710364253
44713661746637196510174346769365511584466954883; 1.271337872187865731012978
478823867003726651653645995407546118684087277042506695 0.126160049694481876
7249738495845282735927546111132028070493637791268222071273071]
 [-1.9892926282818617213240352055427900600018319298058145923659045998012856
54028714 -0.235111120784467091849625313943648820159085439173821665548482789
9847587802223395; -0.813371225913394308521073340158755896264088157703000247
8340774692927415972822211 0.64580271757995047445474083578184225196994747568
30282919468283988048549068462689; 1.491572332031762788566104404405841023565
501199691316489350905513981020706849705 0.772784643555881009599750648971518
4248589311334653244853686473727688680587681075; 2.0895260347720888111419002
93447729680719551911839150040105633776524782622527076 0.1077920688514654991
249483638957735418886356257930645139826182621601987936169665]
 [-5.8436678475230145559003854599723534008955592214399688473743480348559670
40730172 -0.244621598701266659412815405384677412649964563254688516098025581
6090257332568621; -1.707192889125726061192346774799607909714152473589958414
847727749571644452094953 1.011739730724958346516282298631557561506983002999
382983968892012679679665624513; 6.79451521497870317413826390778699003480897
7099874321061793239181996029918146766 1.10374903998600540887165121048819952
8122411881177678714299407608197950715868868; 0.9229001986442570320090273422
002826170810533780410310038624228847840313695021898 -0.27484129971364435695
51732426008863968723967789200349515204777196169994997521221]
 [-5.8605010633542207161236924181024214441631109165756855341676716742021308
68055959 0.7644711886766404233745864222366199877188140868129878084792735637
017321907148534; 4.74724138525243443495808168532036148717147951382437163082
8072117068754136299995 2.13785988170593997345873110429821526170247409892348
074551181062742286652617896; 14.7543422896124753852473461054474867628088907
9392676665686535910762052379856644 0.72566897389483569085142652331044486534
06569027086927786455719640515335201107109; -4.19139067862136911139645267600
5122693145159813989552933040222072328265793321175 -0.9238562474921241195131
32605179671938764041694184013554431872149434086221214728]
 [-5.3744128793052812184456535476578691686998539935544727307983989440580883
99021822 0.9249765388304762539468899210499841827662511126108836381305616581
395393901378754; 5.90957132286666258627381236643856433428560567417523805948
4723055014279238725487 2.26649528024961185300603818488951248912132517822791
3594899289512245332760478566; 15.271333102751646217641691488776394852140018
60815484588117891702797151001202751 0.6014570411302269580019855254582193867
445483775727938036664581254761801864414245; -4.8309053298570960090549109490
72624547520463158849500811539963393188081307745328 -0.972721147184919753984
080364927216402798425139699116221968551185489280069931983]
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
 0.03205239800732855
 0.09895426141282566
 0.1823778419174073
 0.2828280692410671
 0.40449814208297546
 0.5457855570638858
 0.7045526160891851
 0.856138015131654
 1.0
u: 10-element Array{StaticArrays.SArray{Tuple{4,2},Float64,2,8},1}:
 [0.5804069701197787 0.6081566479542033; 0.9418641877739737 0.0481315644732
041; 0.7092190584215454 0.6294656779019507; 0.07430230080546174 0.595508208
7592192]
 [0.58035488061506 0.5155469883427757; 1.0282947861540435 0.125598587430595
; 0.6384024027164247 0.5792156890453118; 0.15862052896702722 0.785276933962
4651]
 [0.5365781476545716 0.21440118695163585; 1.1366111386896796 0.157719362397
953; 0.5039298469630689 0.5435714595661453; 0.3262219863912431 1.1603404715
170809]
 [0.4005489066481007 -0.36146233252069543; 1.1437208813267146 -0.0207585184
51784086; 0.38025813298918343 0.676328262279411; 0.5110010379341184 1.56409
24369417673]
 [0.12572279451909574 -1.3230196745037563; 0.9933631123650289 -0.4723030050
9434076; 0.33117516765677096 1.181581124881154; 0.6813476474737132 1.905664
8972257086]
 [-0.3380972505296776 -2.7828271755195564; 0.6570132626483132 -1.1464402771
714133; 0.4684506923757258 2.413198019450224; 0.7833627573684587 2.01662778
98309652]
 [-0.9731465709232148 -4.618705973388045; 0.2229170315922806 -1.63644774445
98433; 0.945498453166232 4.769138091097415; 0.7171657273661947 1.5866294378
86209]
 [-1.6302752380850851 -6.276077662813274; -0.0268886367410181 -0.9984804244
18369; 1.8712796590107772 8.428461738252519; 0.3584308449045448 0.213806569
0550377]
 [-1.952912467467896 -6.634706900240984; 0.28394044297409693 1.640077555622
5646; 2.9889379153279836 12.29874458840961; -0.2874779700449619 -2.06136048
40772743]
 [-1.7201214589663556 -4.99596692502776; 1.3051093798911246 6.5574719311210
3; 3.976240515854327 15.257533807546377; -1.1523744802648688 -5.01009305632
1093]
````



````julia
sol[3]
````


````
4×2 StaticArrays.SArray{Tuple{4,2},Float64,2,8} with indices SOneTo(4)×SOne
To(2):
 0.536578  0.214401
 1.13661   0.157719
 0.50393   0.543571
 0.326222  1.16034
````





## Conclusion

These are the basic controls in DifferentialEquations.jl. All equations are defined via a problem type, and the `solve` command is used with an algorithm choice (or the default) to get a solution. Every solution acts the same, like an array `sol[i]` with `sol.t[i]`, and also like a continuous function `sol(t)` with a nice plot command `plot(sol)`. The Common Solver Options can be used to control the solver for any equation type. Lastly, the types used in the numerical solving are determined by the input types, and this can be used to solve with arbitrary precision and add additional optimizations (this can be used to solve via GPUs for example!). While this was shown on ODEs, these techniques generalize to other types of equations as well.


## Appendix

 This tutorial is part of the DiffEqTutorials.jl repository, found at: <https://github.com/JuliaDiffEq/DiffEqTutorials.jl>

To locally run this tutorial, do the following commands:
```
using DiffEqTutorials
DiffEqTutorials.weave_file("introduction","01-ode_introduction.jmd")
```

Computer Information:
```
Julia Version 1.4.2
Commit 44fa15b150* (2020-05-23 18:35 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Xeon(R) CPU @ 2.30GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-8.0.1 (ORCJIT, haswell)
Environment:
  JULIA_CUDA_MEMORY_LIMIT = 536870912
  JULIA_DEPOT_PATH = /builds/JuliaGPU/DiffEqTutorials.jl/.julia
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
