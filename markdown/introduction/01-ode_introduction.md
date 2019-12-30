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
 1.9730384275623003
 2.664456142481452
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
 1.9730384275623003
 2.664456142481452
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
 2.6664820303831105
 3.664456142481452
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
 1.9730384275623003
 2.664456142481452
````





you see that it says that the solution has a order changing interpolation. The default algorithm automatically switches between methods in order to handle all types of problems. For non-stiff equations (like the one we are solving), it is a continuous function of 4th order accuracy. We can call the solution as a function of time `sol(t)`. For example, to get the value at `t=0.45`, we can use the command:

````julia
sol(0.45)
````


````
1.5542610480553116
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
 1.1029627851292922
 1.2165269512238264
 1.341783821227542 
 1.4799379510586077
 1.6323162070541606
 1.8003833264983586
 1.9857565541588764
 2.1902158127997704
 2.4157257420844966
 2.664456142481452
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
 1.9857565541588764
 2.4157257420844966
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
 1.9730384275623003
 2.664456142481452
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
 2.664456142481452
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
 2.664456142481452
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
 2.043449143475479 
 2.6418256160577602
 2.6644526430553808
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
 2.0905717422346277
 2.4861021714470244
 2.6644562434913373
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
 0.04841602236644631
 0.12841436984565474
 0.21949854815526873
 0.3286354523115623 
 0.45363089733018575
 0.5884696836694016 
 0.717622136733705  
 0.8905885900507794 
 1.0                
u: 10-element Array{Array{Float64,2},1}:
 [0.815796610273396 0.3545669953455124; 0.4267258979958011 0.96717269152610
36; 0.7372938700681682 0.5956927989343674; 0.11307704703071297 0.5421638672
32089]     
 [0.7982516219133012 0.22247112299856342; 0.6265883479028581 0.954691017339
085; 0.5914099176816812 0.5684827098999932; 0.3562759911599949 0.6631347835
45533]     
 [0.6356696207300371 -0.0699302148681995; 0.7697844516711172 0.841288679626
7256; 0.40170368675372103 0.601560360496467; 0.74488658768305 0.82313885049
45347]     
 [0.24550896421170043 -0.49746826290108265; 0.6689314152817034 0.6127000948
258978; 0.3210023164355822 0.7827470326290281; 1.1410050499197564 0.9286284
03914769]  
 [-0.4934606559541962 -1.095635211018943; 0.25118935900975836 0.29153200894
323766; 0.5105991446378327 1.2306136014319375; 1.5016594168676047 0.9202440
087549442] 
 [-1.634366624658973 -1.7940823881816685; -0.44321650040428384 0.0297190129
9240043; 1.2339131101088978 2.0587100077032434; 1.6840174430711865 0.693580
4508152102]
 [-3.0480311116087115 -2.3825820922429375; -1.11438993032307 0.122475324965
65858; 2.713300622574013 3.2585413126275444; 1.5008487473767809 0.154619794
2334656]   
 [-4.315674439761319 -2.56010175573077; -1.2664137454639106 0.8176908537029
307; 4.783920249052397 4.5274176666857935; 0.8621017015588281 -0.6570117608
688686]    
 [-5.2551724564900235 -1.7760576760312996; 0.033198786240424605 2.974658736
900867; 8.187045156049942 5.873036379501265; -0.8050964994094623 -2.1324714
860548752] 
 [-4.979066381210218 -0.4478110029627467; 2.1075802094039346 5.105843203753
425; 10.289435170877807 6.094978361813493; -2.3379487634840435 -3.199628982
155013]
````





There is no real difference from what we did before, but now in this case `u0` is a `4x2` matrix. Because of that, the solution at each time point is matrix:

````julia
sol[3]
````


````
4×2 Array{Float64,2}:
 0.63567   -0.0699302
 0.769784   0.841289 
 0.401704   0.60156  
 0.744887   0.823139
````





In DifferentialEquations.jl, you can use any type that defines `+`, `-`, `*`, `/`, and has an appropriate `norm`. For example, if we want arbitrary precision floating point numbers, we can change the input to be a matrix of `BigFloat`:

````julia
big_u0 = big.(u0)
````


````
4×2 Array{BigFloat,2}:
 0.815797  0.354567
 0.426726  0.967173
 0.737294  0.595693
 0.113077  0.542164
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
 0.06045369748327955
 0.27340205261345235
 0.5673380394454617 
 0.8768404355382271 
 1.0                
u: 6-element Array{Array{BigFloat,2},1}:
 [0.81579661027339600565255750552751123905181884765625 0.354566995345512392
390219247317872941493988037109375; 0.42672589799580107339238566055428236722
9461669921875 0.9671726915261036339188649435527622699737548828125; 0.737293
8700681681556403646027320064604282379150390625 0.59569279893436744011125938
41389752924442291259765625; 0.113077047030712973807453636254649609327316284
1796875 0.5421638672320889806854893322451971471309661865234375]            
                                                                           
                                                                           
                                                             
 [0.78453389874846221553898335297811719141303603004981581706144464707825965
61517591 0.1841749122174734180241275479000261509497673090472138310534341926
311468789467257; 0.66289317810761029049487789214463061749720289709590151908
45177233998175043537675 0.9444993386789167249801366675163326962012750949260
582747200538849669245108075166; 0.55793828025090878222453612694682815867414
61681071456774240397081876119879033652 0.5668332492220226795504899048145930
519547314271210827451259106082885998826755636; 0.41617606457727363023482486
31760962171710002822562196060823948248819152461176165 0.6906256051268396980
888300031156988821968575668427301615023624170649120553585027]
 [-0.0846775097251102490215945734514424983135383115152048514415067727032320
6945534579 -0.7852325199444285523384769766280161716708042917964239993569081
78799827062017212; 0.496130049309447339918811571278818096585779437286381683
3489406005459553280927505 0.45325439581016648378211265332080028219163727635
68128176029943519868466353976442; 0.369356234395109586302510902031643892800
053092689477833356653491630320274407486 0.971330923520600407789640481087406
3769105022292776695526018883337531568230656202; 1.3381877231392475513574538
67187173753850073806010790721514266640156788723132819 0.9446542402378955135
193476780711685838892016498853817621130203359001492331365822]
 [-2.8241813672739062052389815803089568071499473863983306118006207255559859
89464914 -2.310996071903638597027703842205896752066037337118207438409392337
843017297399427; -1.0318783969991041346093522398653294127948891392330533496
9374588650243859238759 0.07130296797066886352178318621967223900932563271262
166328116249592655632542726982; 2.43275289347621211235514628123461612126340
2766246593002636537407914880892749721 3.05571370678162579858772963877386998
8519412172668191736862991565971065467966149; 1.5599244746510185125008915351
84119060285338170964566367823817403623361804486123 0.2600973398008936633816
279656602241628024537594658014574534542783947860589812622]   
 [-5.2326311290892955506248935756792903565141615842641474722135039612371034
78309721 -1.892946987985841812288179118782116730410850519171202679245195628
523040833054454; -0.1523480672761063451549827615993995121245504913287306312
969586214536440975871661 2.747965303818258576746976361604343958856338049609
278304862923300999851756972517; 7.90855012010486824089912758278313485427752
1662013064777338678839962611794638316 5.80053824479758964576528788986099173
9171698324319216259653383823817131795958317; -0.637795818656921713530613737
6078148917645669590214302110451275144383269341200675 -2.0028936957520679404
77145676753186727773855373641133956511610022527441255402704] 
 [-4.9790679285952016378395847231417721438369095769478853841145120442551566
63021545 -0.447802280501505450389481043189337484651582666328072672139801900
2489561879848852; 2.1075987426257136029395053914747849296483434819794208828
29858647730237880033385 5.1058625591508303548174482296915080606218190766208
25548663895664438708707017582; 10.28945688427075079299505940692616668301575
201331579571755989755070525734745872 6.094983889033111664924227182886231598
311701562453572216085154309795489378156499; -2.3379601823905082048378497613
79796740077996262451394343837415668057307422128537 -3.199638426555815837054
60485718405994613833573923009201594447063802470051293152]
````



````julia
sol[1,3]
````


````
-0.084677509725110249021594573451442498313538311515204851441506772703232069
45534579
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
        
 0.060453697483279549515813279712810712424266474860592303264226255412244400
07007405
 0.273402052613452341453010149914879461747934153856765188346236511356383092
9721304 
 0.567338039445461579355857095764525119958745780825868581004452504045682965
6861064 
 0.876840435538226928445764018292903152398656739455686996795825674545911093
5040272 
 1.0                                                                       
        
u: 6-element Array{Array{BigFloat,2},1}:
 [0.81579661027339600565255750552751123905181884765625 0.354566995345512392
390219247317872941493988037109375; 0.42672589799580107339238566055428236722
9461669921875 0.9671726915261036339188649435527622699737548828125; 0.737293
8700681681556403646027320064604282379150390625 0.59569279893436744011125938
41389752924442291259765625; 0.113077047030712973807453636254649609327316284
1796875 0.5421638672320889806854893322451971471309661865234375]            
                                                                           
                                                                           
                                                              
 [0.78453389874846221488012813627899349705194027142182293639557139524773841
19760084 0.1841749122174734163627145399135445035116620675061252727225516104
907176956141425; 0.66289317810761029191570074371895024799124307539108499015
42563886103111865522497 0.9444993386789167244938307777049449450868877912511
47837845386765960328698824054; 0.557938280250908780841127184039081944416435
0235096653844893295198402603628729691 0.56683324922202267952707322890967866
03170943929746262976935869700262615476846257; 0.416176064577273632756350242
6336902150897515310938304046681226325474143380512241 0.69062560512683969922
59761260987142662474028883074650402994297572003474349851551]  
 [-0.0846775097251102313966966587556028353570103311123795789512867236216042
2328982887 -0.7852325199444285380096227618831977773878942499632069115632664
278814205737348195; 0.49613004930944734998062934000171636884620609801520508
3726140460710928627355574 0.45325439581016649157558351156060397900369799198
34105134776941526161094934537715; 0.369356234395109581940523516094936942536
1151877849918004099630227436406875784989 0.97133092352060039716211890057306
87907685867685366142095718790909088753849422992; 1.338187723139247542675534
540998782479478912437469202962207115175693044127715264 0.944654240237895513
6651380518969051114538188284206022945382533235641678416325072]
 [-2.8241813672739057597314552705338300010168836593828443851901541155812251
47546352 -2.310996071903638445580755813611681745460210337645353868735673316
190890748487608; -1.0318783969991039592491676670567250038413868547098594239
36884244081174966882921 0.0713029679706687773048801055572617883177887197288
8991755295488333120091106649332; 2.4327528934762115732135510611748400375523
91605235430529475713669108592373126095 3.0557137067816254000355716144182102
18818067437947698948604592188864739194372311; 1.559924474651018617835807921
535147809849178216332401461119121225191597196399308 0.260097339800893864914
9029786988559491170338942197944734429738037344507983764936]   
 [-5.2326311290892953212132571644994826410856173075305687637833137112173867
42618145 -1.892946987985842723976051912437666859618380865789935787573044098
824453906302059; -0.1523480672761077956973888944822588569039874711281577028
028863110683665642919559 2.747965303818256764588819225017498226928383413594
280013476189003789037102595366; 7.90855012010486596292301651700762740566952
0385529814165533496068004049481879541 5.80053824479758902062481706965453241
5472330486811972109474424685828814606507908; -0.637795818656920371541088075
6844183677374026389856155422906369830058415644590413 -2.0028936957520668888
34620177730737792086840865352599650431225243985504729582682]  
 [-4.9790679285952012187541502109092060960922684884016566542141281981951370
58819989 -0.447802280501504479266980737110267283838110245554632687720931377
3867632507524396; 2.1075987426257151042534701997702667642095781331428343708
69666090586920555446921 5.1058625591508317272145030042532220995041342253579
32992592590725964808152707638; 10.28945688427075189076163985400129153757320
254794146611711502885318164179543556 6.094983889033111576967926294487602847
936135224788700249445201084563875910205805; -2.3379601823905091756553011358
48615086114228325654235362408126707833221420163013 -3.199638426555816452792
373928308093548958105305825776884830968576591645009052217]
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
 0.04324956523726192
 0.11321071740160932
 0.19766180313092657
 0.30258838818651473
 0.41227782936810903
 0.5472917650344715 
 0.6836815256385428 
 0.83571552907508   
 1.0                
u: 10-element Array{StaticArrays.SArray{Tuple{4,2},Float64,2,8},1}:
 [0.24599623348395272 0.1814108490131281; 0.7708739473919768 0.529358641479
1931; 0.2875377125306622 0.9317905900644217; 0.7375591147725589 0.951643897
7048698]   
 [0.08262806295193584 -0.04183421547170835; 0.6839400434743584 0.5236711713
173474; 0.29281353341537897 0.9642819812371652; 0.837663950789166 1.1392217
54034667]  
 [-0.23918483474994076 -0.5050271924926594; 0.48253716125752355 0.405843814
542289; 0.37649487916335544 1.1269321265179466; 0.9655816072157422 1.392108
0695872687]
 [-0.707356115818706 -1.2107082737819224; 0.17502667550873519 0.13554490669
329494; 0.6199786613335732 1.5395352836258278; 1.051298907574082 1.58986983
20170663]  
 [-1.368199607680919 -2.2454865535343576; -0.21983680350695328 -0.266165376
3342805; 1.164008997882098 2.4312843658344923; 1.0286328220119376 1.6244650
963768992] 
 [-2.0739169564939326 -3.386352082229682; -0.513349635226205 -0.54999398015
93023; 2.0237638377535716 3.8340200867094527; 0.821539450526163 1.351058091
9520478]   
 [-2.793480456055117 -4.588498937673334; -0.46913795883528414 -0.3127794049
2284995; 3.4280007077225534 6.13769861923577; 0.2706288031634154 0.50214236
43415781]  
 [-3.0963872009705744 -5.135754407467346; 0.3105135016542419 1.068674934858
6093; 5.0398497142145455 8.807426511676754; -0.6382575575385778 -0.97801097
35271589]  
 [-2.5465669507537756 -4.283558815399139; 2.3351470437134223 4.478239144417
056; 6.615096625572361 11.449932032919943; -2.0347484151890463 -3.321401138
795485]    
 [-0.4173836402111992 -0.7717855393236781; 6.005015507560791 10.63324567229
428; 7.209465218424538 12.487348154949128; -3.827977037104813 -6.3971076701
67191]
````



````julia
sol[3]
````


````
4×2 StaticArrays.SArray{Tuple{4,2},Float64,2,8} with indices SOneTo(4)×SOne
To(2):
 -0.239185  -0.505027
  0.482537   0.405844
  0.376495   1.12693 
  0.965582   1.39211
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
