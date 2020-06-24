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
t: 11-element Array{Float64,1}:
 0.0
 0.04607057407067815
 0.1124668145681987
 0.1821322854968907
 0.27290665707088113
 0.38822018198766434
 0.5208308076580936
 0.643747804384864
 0.785728049362553
 0.9551021713322809
 1.0
u: 11-element Array{Array{Float64,2},1}:
 [0.7941532891032377 0.4907734093293066; 0.5862181620470468 0.9455968500164
746; 0.17738770304358664 0.9271748093737235; 0.7438330652698997 0.572531308
0506806]
 [0.6276311358533029 0.35722914295528574; 0.5682531316110733 1.007931562807
279; 0.08530979883863311 0.8789494696303977; 0.9865032644617893 0.755546676
9971214]
 [0.2761469047001174 0.08190797533615646; 0.40922861319528 0.99122315447337
52; 0.03884261294722911 0.876614598527304; 1.3098877852101074 0.98687156123
89992]
 [-0.23000099518297468 -0.30374667127066973; 0.09681577222203658 0.86120757
85832805; 0.13095282834491492 0.9810622182569285; 1.6001764982037967 1.1758
123515847376]
 [-1.0812281457148751 -0.9292729041735204; -0.4708810435794842 0.5764610956
11605; 0.5215662579937026 1.315360996510325; 1.872193503717329 1.3157626303
991383]
 [-2.4056927920844453 -1.8487150076708467; -1.2871104494811305 0.1703590781
570687; 1.5425064471823777 2.104828641997566; 1.981603611173861 1.276084917
6482931]
 [-4.0765400335903 -2.898116140145526; -2.0014967755905984 -0.0699032475790
351; 3.5090556510040787 3.51626244563223; 1.6760471818246647 0.865211473707
7802]
 [-5.480917797680895 -3.6159027799712002; -2.01896531971339 0.2480010846624
4825; 6.034180566884379 5.196215933802348; 0.8837186347882154 0.08938142426
766538]
 [-6.441027631393474 -3.753997249725007; -0.6945994516424794 1.649670784879
0688; 9.487396518871714 7.255076106962838; -0.7239426815081805 -1.289723982
4006614]
 [-5.744628878475986 -2.303742281673544; 3.512929443628275 5.17072862725652
6; 13.406839080073667 9.028167314944797; -3.606276777438641 -3.499456119340
925]
 [-5.07730655915073 -1.5322359772549357; 5.177298567770851 6.45719382049693
5; 14.199891218877683 9.203665415216822; -4.51978430231523 -4.1495596973658
03]
````





There is no real difference from what we did before, but now in this case `u0` is a `4x2` matrix. Because of that, the solution at each time point is matrix:

````julia
sol[3]
````


````
4×2 Array{Float64,2}:
 0.276147   0.081908
 0.409229   0.991223
 0.0388426  0.876615
 1.30989    0.986872
````





In DifferentialEquations.jl, you can use any type that defines `+`, `-`, `*`, `/`, and has an appropriate `norm`. For example, if we want arbitrary precision floating point numbers, we can change the input to be a matrix of `BigFloat`:

````julia
big_u0 = big.(u0)
````


````
4×2 Array{BigFloat,2}:
 0.794153  0.490773
 0.586218  0.945597
 0.177388  0.927175
 0.743833  0.572531
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
 0.15430179028206634
 0.40356148993249996
 0.7103297612192179
 1.0
u: 5-element Array{Array{BigFloat,2},1}:
 [0.7941532891032376806350612241658382117748260498046875 0.4907734093293065
757990234487806446850299835205078125; 0.58621816204704679797998778667533770
20359039306640625 0.945596850016474643751962503301911056041717529296875; 0.
1773877030435866419821877570939250290393829345703125 0.92717480937372354077
297131880186498165130615234375; 0.74383306526989967366603195841889828443527
2216796875 0.57253130805068064290708207408897578716278076171875]
 [-0.0114699974973189404115027131132974520883943206398908031440357944419537
9401829513 -0.1386119005901841679092302633844412972279876464353025057116785
199745833576137993; 0.23739114277680506116225653957728100060022493446745232
39898297943918086684181474 0.9249909911322179663051824598707513139216518264
892461971535221512284899907892996; 0.07466519009435353220757384858893622319
777205352857632414172096137258516275933286 0.924768138379166422364606849233
0431936571379648608566458566082803871324662584682; 1.4913942085606147935131
20552829935066038813934571403996289918041989598376916776 1.1078506154991132
98679535787764081165798065420539115910735238484301945319299405]
 [-2.5957279825725918222299502762112005986440620812196828760552540863465457
57274609 -1.974952528666520776484875852638844198030883252863049742192862880
998006875313771; -1.3902413248670703460890626896633956115533723030245721505
95939241642850427077915 0.1237116208267421282568731755461698329295358800497
937210601265909390264967392005; 1.72628725519265558888063584601056218802749
346406486777691590417364395558733768 2.241531078724695902890809029216693371
239253763373023224474554323315437143075917; 1.97200122279536814430473855086
6236117932300488118852869970636404718683248177782 1.24973554951010001116382
7831037529549679985722094657164796292646891618681980606]
 [-6.0551444604027631684369266943629870212559786130419001258559648678843305
49844528 -3.799442849609254652666324218884501010758512693454455388951367233
767272059735532; -1.6122875206229236362535901738675168870619207785744820622
75089864374475917153173 0.7458325815295959274783338135555207934795550806941
712756482689851752257632960027; 7.61222369511939228249820288103933617754216
9869777554489832815235590189773696319 6.17485295405117932642392601265224875
7457320789576977971877185563340225172764634; 0.2249604362032050238721503025
840102097068480249371381190631600554119222506386768 -0.49486380465149302988
95771169162354672198542440287748393485444609508801087134184]
 [-5.0773019391061987340321955245061097284414648042022449291586130700336070
55886608 -1.532226398591335270104745568561844772038909552994732214250212084
306903603620918; 5.17732832764655302663957321534779536572105105079983397335
1106429128741692020324 6.45721928527921474305467381274437329682819695546661
3336582943192556835271717874; 14.199915178241088065216125871300357568127730
5879705047808358162998383751844336 9.20367562011533386701587378354910357351
400169798491011215458616566164452286847; -4.5198010150846592178183998833644
83844634449056442769218436899545068761532776986 -4.149572356245713365596243
089869062248960229612680947536189110210056709031356784]
````



````julia
sol[1,3]
````


````
-2.595727982572591822229950276211200598644062081219682876055254086346545757
274609
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
 0.154301790282066334282984754228530361635045244524561270578218607216071593
9716733
 0.403561489932499970784037349339584588677311990157368071715104766286773833
8601471
 0.710329761219217814618631470305626707466802557043095448071604716145630073
9301024
 1.0
u: 5-element Array{Array{BigFloat,2},1}:
 [0.7941532891032376806350612241658382117748260498046875 0.4907734093293065
757990234487806446850299835205078125; 0.58621816204704679797998778667533770
20359039306640625 0.945596850016474643751962503301911056041717529296875; 0.
1773877030435866419821877570939250290393829345703125 0.92717480937372354077
297131880186498165130615234375; 0.74383306526989967366603195841889828443527
2216796875 0.57253130805068064290708207408897578716278076171875]
 [-0.0114699974973189221010049351458049857778324213847688158268458792490327
0149765141 -0.1386119005901841539887167986676404398166819836273381365140438
645413717946849353; 0.23739114277680507267596720478823601599260571011517704
03202321861350000978139759 0.9249909911322179712794911914103237402728506476
920794756379635139325956421091698; 0.07466519009435352843861381355625982961
752340144015590269922148604603100786295877 0.924768138379166418289123403213
9396201421648784489760910010303717045670745596987; 1.4913942085606147834822
1618979005963792981780643503615452297226788841437047003 1.10785061549911329
2231407877086494252351670537764253491018439414910079205520004]
 [-2.5957279825725915909446099763768008762078983941320612235986808469388631
32679397 -1.974952528666520623783709668377033703475787361246021345492877098
68880171100826; -1.39024132486707022328974877972601295484794584033552018036
7718138804776645343298 0.12371162082674218266868353361246720772165454958618
77274021360261235976962965772; 1.726287255192655359467396374002807703537410
086125153140801875020612991234250752 2.241531078724695732996704721446078975
921253410457079896062149747322173277122928; 1.97200122279536815970884791722
3662208197893479178156342597361726751554535952994 1.24973554951010004625691
1744932293415798880004808339291829088393878831004938507]
 [-6.0551444604027625711617372341991225588112101077812847922900464597849172
97621583 -3.799442849609254542433860367518045289866626213246846073888122376
700283572782045; -1.6122875206229243664664589511759701347831908657343467587
74940226085506800453871 0.7458325815295951376576535175114240564412424515755
570349392163052058087127018548; 7.61222369511939024895729832222222122377103
6492834837082444488508232937535434912 6.17485295405117810333925248235201345
6596001971744950637271183156728433448636223; 0.2249604362032059515539301180
499457880320481736196537692776401694778083857656278 -0.49486380465149222932
13844941308232582604388887318832146565658145136046637476256]
 [-5.0773019391061992203572112710974711247302234415975367733011400953755009
78809992 -1.532226398591335803445798316861077058961701502105732747844894261
299443606208807; 5.17732832764655192487646495690854381125113634849360454655
8097523510739697384979 6.45721928527921390427910708022287556044422431787907
9113369123777947865839154976; 14.199915178241087626972164378537634087794035
17812980015919667140479770656026898 9.2036756201153338120784252738517117424
91454454833675939773488836431710363456092; -4.51980101508465863770637702046
6406843044748096368229914628132580130731711990742 -4.1495723562457129598947
82919608400456096619949242972042489518282324423328927483]
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
 0.05173733926519386
 0.12968216086995493
 0.2181507789071571
 0.33377183039595937
 0.45903347705795666
 0.599523478893951
 0.7669966844409005
 0.9491428902520377
 1.0
u: 10-element Array{StaticArrays.SArray{Tuple{4,2},Float64,2,8},1}:
 [0.24174231514370637 0.9588713651058418; 0.6714930389257048 0.792093961019
5623; 0.9934337836950857 0.22231721705021368; 0.1161824600446737 0.11835856
725928973]
 [0.20916618198865156 0.9510037416277852; 0.8160929902766707 0.894604239742
1133; 0.9551347559877428 0.035260011767177746; 0.22511190179953078 0.326297
65025773394]
 [0.10570284184988046 0.8337484301856629; 0.9505738150284061 0.893221579003
4673; 0.9276458416039985 -0.2087975574517865; 0.3654522141473526 0.63480830
76894289]
 [-0.08140127632463262 0.5415480268328536; 1.0011668688850042 0.67766913424
20539; 0.9591168772975129 -0.38630410321240194; 0.47822088695366133 0.96109
5698334373]
 [-0.40583827597685007 -0.09195505835799056; 0.9612263986584619 0.111067883
27375428; 1.128585574534724 -0.3684896675387904; 0.5286896754699721 1.31144
07836581035]
 [-0.7906681231386433 -1.0626593937565882; 0.8924317438517294 -0.7250421422
974695; 1.4907825474210916 0.09008629759693254; 0.4322439246974083 1.532707
2838771854]
 [-1.1327446348328685 -2.3864391270661436; 0.9626133966359122 -1.6512929639
783387; 2.0804359005017963 1.2684130246537533; 0.11291985713218994 1.492569
8001950711]
 [-1.1674822889530043 -3.9716177631281235; 1.518745680102495 -2.17208675706
43723; 2.8487978126797584 3.620423382380854; -0.557978697553646 0.904185110
3614993]
 [-0.37579837691369056 -5.036291633213437; 2.9550601185846865 -1.1464850589
659084; 3.285465781672628 7.016116774022278; -1.5564052588686832 -0.5567971
782024477]
 [0.04941428639581352 -5.081540770528331; 3.516767773176111 -0.439809712370
85775; 3.2335390545194254 8.005136666685008; -1.8565253332994112 -1.1272786
856389185]
````



````julia
sol[3]
````


````
4×2 StaticArrays.SArray{Tuple{4,2},Float64,2,8} with indices SOneTo(4)×SOne
To(2):
 0.105703   0.833748
 0.950574   0.893222
 0.927646  -0.208798
 0.365452   0.634808
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
