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

```julia
f(u,p,t) = 0.98u
```

```
f (generic function with 1 method)
```





with $ u_0 = 1.0 $. If we want to solve this model on a time span from `t=0.0` to `t=1.0`, then we define an `ODEProblem` by specifying this function `f`, this initial condition `u0`, and this time span as follows:

```julia
using DifferentialEquations
f(u,p,t) = 0.98u
u0 = 1.0
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
```

```
ODEProblem with uType Float64 and tType Float64. In-place: false
timespan: (0.0, 1.0)
u0: 1.0
```





To solve our `ODEProblem` we use the command `solve`.

```julia
sol = solve(prob)
```

```
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
```





and that's it: we have succesfully solved our first ODE!

#### Analyzing the Solution

Of course, the solution type is not interesting in and of itself. We want to understand the solution! The documentation page which explains in detail the functions for analyzing the solution is the [Solution Handling](https://docs.sciml.ai/dev/basics/solution/) page. Here we will describe some of the basics. You can plot the solution using the plot recipe provided by [Plots.jl](http://docs.juliaplots.org/dev/):

```julia
using Plots; gr()
plot(sol)
```

![](figures/01-ode_introduction_4_1.png)



From the picture we see that the solution is an exponential curve, which matches our intuition. As a plot recipe, we can annotate the result using any of the [Plots.jl attributes](http://docs.juliaplots.org/dev/attributes/). For example:

```julia
plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
```

![](figures/01-ode_introduction_5_1.png)



Using the mutating `plot!` command we can add other pieces to our plot. For this ODE we know that the true solution is $u(t) = u_0 exp(at)$, so let's add some of the true solution to our plot:

```julia
plot!(sol.t, t->1.0*exp(0.98t),lw=3,ls=:dash,label="True Solution!")
```

![](figures/01-ode_introduction_6_1.png)



In the previous command I demonstrated `sol.t`, which grabs the array of time points that the solution was saved at:

```julia
sol.t
```

```
5-element Array{Float64,1}:
 0.0
 0.10042494449239292
 0.35218603951893646
 0.6934436028208104
 1.0
```





We can get the array of solution values using `sol.u`:

```julia
sol.u
```

```
5-element Array{Float64,1}:
 1.0
 1.1034222047865465
 1.4121908848175448
 1.9730384275622996
 2.664456142481451
```





`sol.u[i]` is the value of the solution at time `sol.t[i]`. We can compute arrays of functions of the solution values using standard comprehensions, like:

```julia
[t+u for (u,t) in tuples(sol)]
```

```
5-element Array{Float64,1}:
 1.0
 1.2038471492789395
 1.7643769243364813
 2.66648203038311
 3.664456142481451
```





However, one interesting feature is that, by default, the solution is a continuous function. If we check the print out again:

```julia
sol
```

```
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
```





you see that it says that the solution has a order changing interpolation. The default algorithm automatically switches between methods in order to handle all types of problems. For non-stiff equations (like the one we are solving), it is a continuous function of 4th order accuracy. We can call the solution as a function of time `sol(t)`. For example, to get the value at `t=0.45`, we can use the command:

```julia
sol(0.45)
```

```
1.554261048055312
```





#### Controlling the Solver

DifferentialEquations.jl has a common set of solver controls among its algorithms which can be found [at the Common Solver Options](https://docs.sciml.ai/dev/basics/common_solver_opts/) page. We will detail some of the most widely used options.

The most useful options are the tolerances `abstol` and `reltol`. These tell the internal adaptive time stepping engine how precise of a solution you want. Generally, `reltol` is the relative accuracy while `abstol` is the accuracy when `u` is near zero. These tolerances are local tolerances and thus are not global guarantees. However, a good rule of thumb is that the total solution accuracy is 1-2 digits less than the relative tolerances. Thus for the defaults `abstol=1e-6` and `reltol=1e-3`, you can expect a global accuracy of about 1-2 digits. If we want to get around 6 digits of accuracy, we can use the commands:

```julia
sol = solve(prob,abstol=1e-8,reltol=1e-8)
```

```
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
```





Now we can see no visible difference against the true solution:


```julia
plot(sol)
plot!(sol.t, t->1.0*exp(0.98t),lw=3,ls=:dash,label="True Solution!")
```

![](figures/01-ode_introduction_13_1.png)



Notice that by decreasing the tolerance, the number of steps the solver had to take was `9` instead of the previous `5`. There is a trade off between accuracy and speed, and it is up to you to determine what is the right balance for your problem.

Another common option is to use `saveat` to make the solver save at specific time points. For example, if we want the solution at an even grid of `t=0.1k` for integers `k`, we would use the command:

```julia
sol = solve(prob,saveat=0.1)
```

```
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
```





Notice that when `saveat` is used the continuous output variables are no longer saved and thus `sol(t)`, the interpolation, is only first order. We can save at an uneven grid of points by passing a collection of values to `saveat`. For example:

```julia
sol = solve(prob,saveat=[0.2,0.7,0.9])
```

```
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
```





If we need to reduce the amount of saving, we can also turn off the continuous output directly via `dense=false`:

```julia
sol = solve(prob,dense=false)
```

```
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
```





and to turn off all intermediate saving we can use `save_everystep=false`:

```julia
sol = solve(prob,save_everystep=false)
```

```
retcode: Success
Interpolation: 1st order linear
t: 2-element Array{Float64,1}:
 0.0
 1.0
u: 2-element Array{Float64,1}:
 1.0
 2.664456142481451
```





If we want to solve and only save the final value, we can even set `save_start=false`.

```julia
sol = solve(prob,save_everystep=false,save_start = false)
```

```
retcode: Success
Interpolation: 1st order linear
t: 1-element Array{Float64,1}:
 1.0
u: 1-element Array{Float64,1}:
 2.664456142481451
```





Note that similarly on the other side there is `save_end=false`.

More advanced saving behaviors, such as saving functionals of the solution, are handled via the `SavingCallback` in the [Callback Library](https://docs.sciml.ai/dev/features/callback_library/#saving_callback-1) which will be addressed later in the tutorial.

#### Choosing Solver Algorithms

There is no best algorithm for numerically solving a differential equation. When you call `solve(prob)`, DifferentialEquations.jl makes a guess at a good algorithm for your problem, given the properties that you ask for (the tolerances, the saving information, etc.). However, in many cases you may want more direct control. A later notebook will help introduce the various *algorithms* in DifferentialEquations.jl, but for now let's introduce the *syntax*.

The most crucial determining factor in choosing a numerical method is the stiffness of the model. Stiffness is roughly characterized by a Jacobian `f` with large eigenvalues. That's quite mathematical, and we can think of it more intuitively: if you have big numbers in `f` (like parameters of order `1e5`), then it's probably stiff. Or, as the creator of the MATLAB ODE Suite, Lawrence Shampine, likes to define it, if the standard algorithms are slow, then it's stiff. We will go into more depth about diagnosing stiffness in a later tutorial, but for now note that if you believe your model may be stiff, you can hint this to the algorithm chooser via `alg_hints = [:stiff]`.

```julia
sol = solve(prob,alg_hints=[:stiff])
```

```
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
```





Stiff algorithms have to solve implicit equations and linear systems at each step so they should only be used when required.

If we want to choose an algorithm directly, you can pass the algorithm type after the problem as `solve(prob,alg)`. For example, let's solve this problem using the `Tsit5()` algorithm, and just for show let's change the relative tolerance to `1e-6` at the same time:

```julia
sol = solve(prob,Tsit5(),reltol=1e-6)
```

```
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
```





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

```julia
function lorenz!(du,u,p,t)
    σ,ρ,β = p
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end
```

```
lorenz! (generic function with 1 method)
```





Notice here we used the in-place format which writes the output to the preallocated vector `du`. For systems of equations the in-place format is faster. We use the initial condition $u_0 = [1.0,0.0,0.0]$ as follows:

```julia
u0 = [1.0,0.0,0.0]
```

```
3-element Array{Float64,1}:
 1.0
 0.0
 0.0
```





Lastly, for this model we made use of the parameters `p`. We need to set this value in the `ODEProblem` as well. For our model we want to solve using the parameters $\sigma = 10$, $\rho = 28$, and $\beta = 8/3$, and thus we build the parameter collection:

```julia
p = (10,28,8/3) # we could also make this an array, or any other type!
```

```
(10, 28, 2.6666666666666665)
```





Now we generate the `ODEProblem` type. In this case, since we have parameters, we add the parameter values to the end of the constructor call. Let's solve this on a time span of `t=0` to `t=100`:

```julia
tspan = (0.0,100.0)
prob = ODEProblem(lorenz!,u0,tspan,p)
```

```
ODEProblem with uType Array{Float64,1} and tType Float64. In-place: true
timespan: (0.0, 100.0)
u0: [1.0, 0.0, 0.0]
```





Now, just as before, we solve the problem:

```julia
sol = solve(prob)
```

```
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
```





The same solution handling features apply to this case. Thus `sol.t` stores the time points and `sol.u` is an array storing the solution at the corresponding time points.

However, there are a few extra features which are good to know when dealing with systems of equations. First of all, `sol` also acts like an array. `sol[i]` returns the solution at the `i`th time point.

```julia
sol.t[10],sol[10]
```

```
(0.08368539694547242, [1.0888636764765296, 2.052326153029042, 0.07402570506
414284])
```





Additionally, the solution acts like a matrix where `sol[j,i]` is the value of the `j`th variable at time `i`:

```julia
sol[2,10]
```

```
2.052326153029042
```





We can get a real matrix by performing a conversion:

```julia
A = Array(sol)
```

```
3×1294 Array{Float64,2}:
 1.0  0.999643     0.996105    0.969359     …   4.52712   8.04367   9.97538
 0.0  0.000998805  0.0109654   0.0897706        6.89588  12.7116   15.1439
 0.0  1.78143e-8   2.14696e-6  0.000143802     16.5854   18.1254   21.0064
```





This is the same as sol, i.e. `sol[i,j] = A[i,j]`, but now it's a true matrix. Plotting will by default show the time series for each variable:

```julia
plot(sol)
```

![](figures/01-ode_introduction_29_1.png)



If we instead want to plot values against each other, we can use the `vars` command. Let's plot variable `1` against variable `2` against variable `3`:

```julia
plot(sol,vars=(1,2,3))
```

![](figures/01-ode_introduction_30_1.png)



This is the classic Lorenz attractor plot, where the `x` axis is `u[1]`, the `y` axis is `u[2]`, and the `z` axis is `u[3]`. Note that the plot recipe by default uses the interpolation, but we can turn this off:

```julia
plot(sol,vars=(1,2,3),denseplot=false)
```

![](figures/01-ode_introduction_31_1.png)



Yikes! This shows how calculating the continuous solution has saved a lot of computational effort by computing only a sparse solution and filling in the values! Note that in vars, `0=time`, and thus we can plot the time series of a single component like:

```julia
plot(sol,vars=(0,2))
```

![](figures/01-ode_introduction_32_1.png)



## Internal Types

The last basic user-interface feature to explore is the choice of types. DifferentialEquations.jl respects your input types to determine the internal types that are used. Thus since in the previous cases, when we used `Float64` values for the initial condition, this meant that the internal values would be solved using `Float64`. We made sure that time was specified via `Float64` values, meaning that time steps would utilize 64-bit floats as well. But, by simply changing these types we can change what is used internally.

As a quick example, let's say we want to solve an ODE defined by a matrix. To do this, we can simply use a matrix as input.

```julia
A  = [1. 0  0 -5
      4 -2  4 -3
     -4  0  0  1
      5 -2  2  3]
u0 = rand(4,2)
tspan = (0.0,1.0)
f(u,p,t) = A*u
prob = ODEProblem(f,u0,tspan)
sol = solve(prob)
```

```
retcode: Success
Interpolation: automatic order switching interpolation
t: 10-element Array{Float64,1}:
 0.0
 0.0435034500129451
 0.11612476188511982
 0.19932266299685725
 0.30384642691933383
 0.42320413461981526
 0.5657506895018174
 0.7194383153589604
 0.8790351034023395
 1.0
u: 10-element Array{Array{Float64,2},1}:
 [0.5390278496717593 0.38648009835056407; 0.22796668993565983 0.40032884348
80601; 0.7145573813827373 0.35412177420434787; 0.41491965964446 0.536517465
7342933]
 [0.4474126803601463 0.2681488503711519; 0.3398165859984035 0.4015258530774
339; 0.6506597235042132 0.3231714471247; 0.6250994565793442 0.6820627375603
341]
 [0.18371868540713732 -0.010055569454728464; 0.38857776316383974 0.30964621
41345637; 0.6131419372712692 0.3409230026858953; 0.9519832436304088 0.89890
62566771611]
 [-0.28357643411819805 -0.44450125242783955; 0.2582447012630125 0.087821281
20662219; 0.7178960524073217 0.4966028924668838; 1.267647987546922 1.091138
0037574459]
 [-1.0924023592064367 -1.135090040836276; -0.10688360296778188 -0.292247311
1549643; 1.1455526979343646 0.9438633187392842; 1.5307792703145302 1.215501
0634304557]
 [-2.2323244229277894 -2.0416162857068367; -0.6152605466993605 -0.708305430
5909887; 2.1212548235055086 1.8429916281780248; 1.5784701947300062 1.147533
3438431627]
 [-3.667345851739975 -3.0932010554902236; -0.9604615818432061 -0.8685646567
375258; 4.009844716643548 3.4516447838331867; 1.172257860480183 0.697683139
9604815]
 [-4.8490200169186695 -3.8180648311572583; -0.39811662918447255 -0.16764239
522095592; 6.7659417679386635 5.646909481507276; 0.05208029271654713 -0.309
31063817053017]
 [-4.963982968059209 -3.566162972993426; 1.9260670303722893 2.0165255242804
254; 9.842743905156954 7.903325247917236; -1.9106005716118946 -1.9365153383
680427]
 [-3.7734670371079084 -2.3022040864455926; 5.199605930692368 4.852333129223
432; 11.66387159516882 9.041121113985; -3.8757631702240865 -3.4882953941880
825]
```





There is no real difference from what we did before, but now in this case `u0` is a `4x2` matrix. Because of that, the solution at each time point is matrix:

```julia
sol[3]
```

```
4×2 Array{Float64,2}:
 0.183719  -0.0100556
 0.388578   0.309646
 0.613142   0.340923
 0.951983   0.898906
```





In DifferentialEquations.jl, you can use any type that defines `+`, `-`, `*`, `/`, and has an appropriate `norm`. For example, if we want arbitrary precision floating point numbers, we can change the input to be a matrix of `BigFloat`:

```julia
big_u0 = big.(u0)
```

```
4×2 Array{BigFloat,2}:
 0.539028  0.38648
 0.227967  0.400329
 0.714557  0.354122
 0.41492   0.536517
```





and we can solve the `ODEProblem` with arbitrary precision numbers by using that initial condition:

```julia
prob = ODEProblem(f,big_u0,tspan)
sol = solve(prob)
```

```
retcode: Success
Interpolation: automatic order switching interpolation
t: 6-element Array{Float64,1}:
 0.0
 0.1331134831839388
 0.38839868455954096
 0.684830725967047
 0.9961979315478988
 1.0
u: 6-element Array{Array{BigFloat,2},1}:
 [0.5390278496717593181841721161617897450923919677734375 0.3864800983505640
719073426225804723799228668212890625; 0.22796668993565982752613763295812532
30571746826171875 0.4003288434880600821230700603337027132511138916015625; 0
.7145573813827372777751634203013963997364044189453125 0.3541217742043478722
28026972152292728424072265625; 0.414919659644459981606701148848515003919601
4404296875 0.536517465734293264034704407094977796077728271484375]
 [0.10228326861112745485553070389744769720258290833394829337563414461377781
22746529 -0.089166590705624325657438436018173779742740862463913052740998287
37890751849541469; 0.376851195814223596169586773238631524346360054505195565
0306306710298031505545346 0.27320053644519862228988584989831737568034411513
62566519051421533599100790688084; 0.620157668366394563788975612128633370179
0849981543443641987813209748161117361396 0.35992027880868300734384081162878
09039471511562779190087920874281021801640063118; 1.022413667998831164395769
647040999408428736305118541065709685498635716712943413 0.943666542957290381
4398724413747897558799973287436808595719286454029545305478708]
 [-1.8839812288676790355725187046092380052280590624660548977089249132965982
88038439 -1.771208599691473175997409551169662126584051668109024142906915702
402636103373277; -0.4716553017672795977450383453958976939523623691008985515
965742485734428143912133 -0.60099373056960719642451983107354824117324324397
69535072789261770009730751024986; 1.779472160437749583016464530100968209388
262302311492114076753329247123614741166 1.536782368484084819690162400620559
821609758732753098733295288831970832858040864; 1.59709179809514948933811473
5308495189304427635129514070472783090573561480565725 1.19374958459272908321
8173761600781353985035810690528056788655793531447565159723]
 [-4.6474472341532992224330809788895847878331677667505459160630494058814275
26663023 -3.716785486943640584327587991864272774956912138173899478320947370
945857449974818; -0.6449191641221221488151679216750180518852701569001534999
584414140833185767923404 -0.43014216877018618342418345242798698899542042779
99247083410765902243513782333473; 6.100674710316737128385409345785800655353
213116276089963118164069990979623434925 5.130726279436337797317609966798346
339736352032553318513006995517299660604073792; 0.37038356674040398684065044
69581350165906950460269884966549252588541356164275061 -0.033273083991408895
42360742743066086666557616806571483113446921972452476812182355]
 [-3.8320480648080671470902721770187978182724725852333504734010020266239600
91538642 -2.359158174753792756300436448783403954343384607772483383136126741
187852103541448; 5.07564699565013672195908107367682464617758930937625276315
6577277622250898009879 4.74747606561705636253511320664689141800111216587674
4880720107910214115173691657; 11.620663656311508131336986193349313566432637
0725099528279746723244973441094635 9.01885082711463107443987508282085865916
4740840381473750223026342370974192905036; -3.809119745643721416777274449585
955213624872457276662109101872107257080553796311 -3.43667192333413546373982
1617791420977100934269151113947500396666958629475634869]
 [-3.7734610190424858359257484186834214904398293505078216905539343146977276
52831761 -2.302197051560771483301000854476534725398783202307495885501703649
642113926645434; 5.19963321437864573435655987430447937195721868779458689841
8688292870492791062273 4.85235605452318276326227033407245636547862845817434
3470917072545602047901715924; 11.663889706714553359743178231044178643496900
87525785496583614719647300639757844 9.0411333491996389024070766023449508487
94204248326097007266401013232326325505722; -3.87577754412719611270852468707
9776513196181158731702128381231351715781167364548 -3.4883071576857996891045
76041885208023578865223984296123033993130073340829104117]
```



```julia
sol[1,3]
```

```
-1.883981228867679035572518704609238005228059062466054897708924913296598288
038439
```





To really make use of this, we would want to change `abstol` and `reltol` to be small! Notice that the type for "time" is different than the type for the dependent variables, and this can be used to optimize the algorithm via keeping multiple precisions. We can convert time to be arbitrary precision as well by defining our time span with `BigFloat` variables:

```julia
prob = ODEProblem(f,big_u0,big.(tspan))
sol = solve(prob)
```

```
retcode: Success
Interpolation: automatic order switching interpolation
t: 6-element Array{BigFloat,1}:
 0.0
 0.133113483183938807201330765667301604334443571094789747354191779420180035
3616041
 0.388398684559540902866767339950787752520659977425419672717096513630516871
712286
 0.684830725967046954864338456867805990885293651322234627006957635860209903
3454134
 0.996197931547898744610352459951323575574720538736483528564123693067586360
9303861
 1.0
u: 6-element Array{Array{BigFloat,2},1}:
 [0.5390278496717593181841721161617897450923919677734375 0.3864800983505640
719073426225804723799228668212890625; 0.22796668993565982752613763295812532
30571746826171875 0.4003288434880600821230700603337027132511138916015625; 0
.7145573813827372777751634203013963997364044189453125 0.3541217742043478722
28026972152292728424072265625; 0.414919659644459981606701148848515003919601
4404296875 0.536517465734293264034704407094977796077728271484375]
 [0.10228326861112748537026497932679500967546998443272748814816219667489799
44001627 -0.089166590705624296374832147232657431211996205236251705289846984
54379875802123162; 0.376851195814223601841426802262813000110048903575605085
4118043926331286323326603 0.27320053644519863626505176360192262882333428542
76265720659083572206671436434579; 0.620157668366394560053467185193725176069
6184981962264797071943000589776596951212 0.35992027880868299942347845554939
07978217028909016988754005442559962038245269339; 1.022413667998831139634100
043224767937705135684398362034476541883183551896249504 0.943666542957290365
8553341871824892498326896404089079001579589117149834263618407]
 [-1.8839812288676787144729187857890559713594267239684648566954882255361863
68236399 -1.771208599691472924179992729143751233612021849948534134823846298
42348070881541; -0.47165530176727945895183742530580274844378214249624822138
79676452579872670914269 -0.600993730569607088507913097098398138119569404322
75248770857031958140750556178; 1.779472160437749285876204327559594063607802
845067346632430750991249893098400507 1.536782368484084550348643315292564474
749122544318135423788041778979293045128218; 1.59709179809514949344971957789
6978162310749778232152148519045689857730235493003 1.19374958459272911572824
7779488561717239618738508081097643768311537999881100001]
 [-4.6474472341532989336105130820203153383392594574003198327369767276132520
58266769 -3.716785486943640426551958124908667716501493800647970763292517397
332082981304654; -0.6449191641221224150732943576585409062016659456141711219
209285812261121616428136 -0.43014216877018647742394614911327951221597678427
53096800553755856185663447960672; 6.100674710316736285822197845240397028960
15327623460938253218054036982601688662 5.1307262794363371381214893327790741
65805525655123060215184704354826770909622859; 0.370383566740404370563282997
578022453734596026749787205915033636834128528463088 -0.03327308399140855937
859661079938043076702517362443449895945916284910505110307044]
 [-3.8320480648080678473416085765103056905383226477142252735930277890313802
62491512 -2.359158174753793438630749664420139536550050657670710111447331252
066885735951842; 5.07564699565013522924248052942535176640185580888510242095
8034975663646927368104 4.74747606561705509888641418833584705486026142295044
0446811339118445525371694983; 11.620663656311507601135637016710426568013697
80164301875582692097371262288000702 9.0188508271146307982729209255820260907
63670971635791537813425548494653354300458; -3.80911974564372061139699564025
5939920935825588526964693238153933022560864696253 -3.4366719233341348394576
46964696843387848702413409828308673474310431203162940894]
 [-3.7734610190424854027881658014737642921491021455253931509626038186430390
48186433 -2.302197051560771063099952156971927609417739717010701495821830481
592658233595505; 5.19963321437864664445581006860940005574749646173591360698
4178033999825162526684 4.85235605452318353253455064629250786232437302135832
9529452417114292685839026707; 11.663889706714553671107113318426696843595518
99018576577041821456262898147695064 9.0411333491996390611823534582499787753
10603450777022214787836378159833747421812; -3.87577754412719660026643813530
7900451557871489300757303936613640771976522920298 -3.4883071576858000665346
76080206391202238336970609203393254159447843735074770321]
```





Let's end by showing a more complicated use of types. For small arrays, it's usually faster to do operations on static arrays via the package [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl). The syntax is similar to that of normal arrays, but for these special arrays we utilize the `@SMatrix` macro to indicate we want to create a static array.

```julia
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
```

```
retcode: Success
Interpolation: automatic order switching interpolation
t: 11-element Array{Float64,1}:
 0.0
 0.0391413918763094
 0.10266328228254122
 0.1724728966332868
 0.2641135700213694
 0.3761524540788478
 0.5120139719843569
 0.664106427674025
 0.8238542822757119
 0.9944937932372135
 1.0
u: 11-element Array{StaticArrays.SArray{Tuple{4,2},Float64,2,8},1}:
 [0.8665948168269284 0.527373635869063; 0.3344267285851159 0.17770988311023
86; 0.26957803094920374 0.7946659774538571; 0.9764659482977645 0.8127403673
11275]
 [0.678932725950976 0.3620825958427153; 0.334255838446851 0.242750738082616
9; 0.19141525178194232 0.760910563771684; 1.2506887837999667 1.054155096447
2248]
 [0.24470156884274463 -0.019576877051132036; 0.1837499023451814 0.217273632
47258763; 0.16357274005557548 0.79311127981975; 1.666557613525322 1.4147742
248109585]
 [-0.4130637381673363 -0.5943127436113538; -0.16867106538729043 0.027932673
3834121; 0.3133865669676923 0.9860226131918394; 2.0620460179590543 1.747917
5793697928]
 [-1.5390006240196574 -1.5676392442575506; -0.8476107581543646 -0.401912811
51317576; 0.8704197215554264 1.5505627883747037; 2.4404649787975567 2.04534
42295279864]
 [-3.2342233576915804 -3.005319339865909; -1.805274322920641 -1.01800031091
0972; 2.2155324515346972 2.804332081893478; 2.6034181044695326 2.1193717419
922367]
 [-5.5034761430137085 -4.8604139914323845; -2.6742415425245603 -1.464784762
868247; 4.923008861671121 5.2084122586190045; 2.20979441744138 1.6545438508
113228]
 [-7.734833594594397 -6.526133715402578; -2.4023193817560164 -0.81299945326
09587; 9.229863496227278 8.870806823210254; 0.8009647556905282 0.2569231877
462672]
 [-8.710546214006884 -6.882358630889965; 0.4307483751364929 2.1368887478890
537; 14.51892302966374 13.111991637870183; -1.9326735730222897 -2.296058606
427409]
 [-6.700526915450677 -4.379210887757312; 7.429323372193585 8.68034237248645
9; 19.328741931072592 16.46112534131533; -6.1984265354751225 -6.09815243497
2047]
 [-6.564250253991363 -4.233165915485553; 7.730395112840296 8.95348502070396
5; 19.44026876023836 16.522025174598916; -6.35470808873006 -6.2340477071660
4]
```



```julia
sol[3]
```

```
4×2 StaticArrays.SArray{Tuple{4,2},Float64,2,8} with indices SOneTo(4)×SOne
To(2):
 0.244702  -0.0195769
 0.18375    0.217274
 0.163573   0.793111
 1.66656    1.41477
```





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
[91a5bcdd-55d7-5caf-9e0b-520d859cae80] Plots 1.6.9
[90137ffa-7385-5640-81b9-e52037218182] StaticArrays 0.12.4
[c3572dad-4567-51f8-b174-8c6c989267f4] Sundials 4.3.0
[37e2e46d-f89d-539d-b4ee-838fcccc9c8e] LinearAlgebra
```
