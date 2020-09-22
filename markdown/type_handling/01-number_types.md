---
author: "Chris Rackauckas"
title: "Solving Equations in With Julia-Defined Types"
---


One of the nice things about DifferentialEquations.jl is that it is designed with Julia's type system in mind. What this means is, if you have properly defined a Number type, you can use this number type in DifferentialEquations.jl's algorithms! [Note that this is restricted to the native algorithms of OrdinaryDiffEq.jl. The other solvers such as ODE.jl, Sundials.jl, and ODEInterface.jl are not compatible with some number systems.]

DifferentialEquations.jl determines the numbers to use in its solvers via the types that are designated by `tspan` and the initial condition of the problem. It will keep the time values in the same type as tspan, and the solution values in the same type as the initial condition. [Note that adaptive timestepping requires that the time type is compaible with `sqrt` and `^` functions. Thus dt cannot be Integer or numbers like that if adaptive timestepping is chosen].

Let's solve the linear ODE first define an easy way to get ODEProblems for the linear ODE:

````julia

using DifferentialEquations
f = (u,p,t) -> (p*u)
prob_ode_linear = ODEProblem(f,1/2,(0.0,1.0),1.01);
````


````
ODEProblem with uType Float64 and tType Float64. In-place: false
timespan: (0.0, 1.0)
u0: 0.5
````





First let's solve it using Float64s. To do so, we just need to set u0 to a Float64 (which is done by the default) and dt should be a float as well.

````julia

prob = prob_ode_linear
sol =solve(prob,Tsit5())
println(sol)
````


````
retcode: Success
Interpolation: specialized 4th order "free" interpolation
t: [0.0, 0.09964258706516003, 0.3457024247583422, 0.6776921908052249, 1.0]
u: [0.5, 0.552938681151017, 0.7089376245893466, 0.9913594502399236, 1.37280
04409033037]
````





Notice that both the times and the solutions were saved as Float64. Let's change the time to use rational values. Rationals are not compatible with adaptive time stepping since they do not have an L2 norm (this can be worked around by defining `internalnorm`, but rationals already explode in size!). To account for this, let's turn off adaptivity as well:

````julia

prob = ODEProblem(f,1/2,(0//1,1//1),101//100);
sol = solve(prob,RK4(),dt=1//2^(6),adaptive=false)
println(sol)
````


````
retcode: Success
Interpolation: 3rd order Hermite
t: Rational{Int64}[0//1, 1//64, 1//32, 3//64, 1//16, 5//64, 3//32, 7//64, 1
//8, 9//64, 5//32, 11//64, 3//16, 13//64, 7//32, 15//64, 1//4, 17//64, 9//3
2, 19//64, 5//16, 21//64, 11//32, 23//64, 3//8, 25//64, 13//32, 27//64, 7//
16, 29//64, 15//32, 31//64, 1//2, 33//64, 17//32, 35//64, 9//16, 37//64, 19
//32, 39//64, 5//8, 41//64, 21//32, 43//64, 11//16, 45//64, 23//32, 47//64,
 3//4, 49//64, 25//32, 51//64, 13//16, 53//64, 27//32, 55//64, 7//8, 57//64
, 29//32, 59//64, 15//16, 61//64, 31//32, 63//64, 1//1]
u: [0.5, 0.5079532157789419, 0.5160329388403366, 0.5242411814636141, 0.5325
799879363893, 0.5410514350635981, 0.549657632684732, 0.5584007241993002, 0.
5672828871006491, 0.5763063335182743, 0.5854733107687577, 0.594786101915468
7, 0.6042470263371675, 0.6138584403056545, 0.6236227375726058, 0.6335423499
657445, 0.6436197479944955, 0.6538574414652724, 0.6642579801065528, 0.67482
39542038958, 0.6855579952450607, 0.6964627765753862, 0.7075410140635964, 0.
7187954667781947, 0.7302289376746193, 0.7418442742933268, 0.753644369468981
6, 0.7656321620509244, 0.777810637635102, 0.7901828293076388, 0.80275181840
02358, 0.815520735257586, 0.8284927600169959, 0.8416711234004085, 0.8550591
075190243, 0.8686600466907208, 0.882477328270475, 0.8965143934939935, 0.910
7747383347635, 0.925261914374735, 0.9399795296888533, 0.954931249743661, 0.
9701207983101929, 0.9855519583913936, 1.0012285731642847, 1.017154546937120
2, 1.0333338461217658, 1.0497705002215465, 1.066468602834806, 1.08343231267
44299, 1.1006658546035855, 1.118173520687937, 1.1359596712645978, 1.1540287
360280843, 1.1723852151335463, 1.191033680317543, 1.2099787760366485, 1.229
2252206241676, 1.2487778074652507, 1.268641406190701, 1.288820963889771, 1.
3093215063422494, 1.3301481392701477, 1.3513060496092948, 1.372800506801159
5]
````





Now let's do something fun. Let's change the solution to use `Rational{BigInt}` and print out the value at the end of the simulation. To do so, simply change the definition of the initial condition.

````julia

prob = ODEProblem(f,BigInt(1)//BigInt(2),(0//1,1//1),101//100);
sol =solve(prob,RK4(),dt=1//2^(6),adaptive=false)
println(sol[end])
````


````
415403291938655888343294424838034348376204408921988582429386196369066828013
380062427154556444246064110042147806995712770513313913105317131993928991562
472219540324173687134074558951938783349315387199475055050716642476760417033
833225395963069751630544424879625010648869655282442577465289103178163815663
464066572670655356269579471636764679863656649012559514171272038086748586891
653145664881452891757769341753396504927956887980186316721217138912802907978
839488971277351483679854338427632656105429434285170828205087679096886906512
836058415177000071451519455149761416134211934766818795085616643778333812510
724294609438512646808081849075509246961483574876752196687093709017376892988
720208689912813268920171256693582145356856885176190731036088900945481923320
301926151164642204512204346142796306783141982263276125756548530824427611816
333393407861066935488564588880674178922907680658650707284447124975289884078
283531881659241492248450685643985785207092880524994430296917090030308304496
2139908567605824428891872081720287044135359380045755621121//302595526357001
916401850227786985339805854374596312639728370747077589271270423243703004392
074003302619884721642626495128918849830763359112247111187416392615737498981
461087857422550657171300852094084580555857942985570738231419687525783564788
285621871741725085612510228468354691202070954415518824737971685957295081128
193794470230767667945336581432859330595785427486755359414346047520148998708
472579747503225700773992946775819105236957926068135290787592745892648489231
548275787132390564752450502531598102790376905344412549120000000000000000000
000000000000000000000000000000000000000000000000000000000000000000000000000
000000000000000000000000000000000000000000000000000000000000000000000000000
000000000000000000000000000000000000000000000000000000000000000000000000000
000000000000000000000000000000000000000000000000000000000000000000000000000
000000000000000000000000000000000000000000000000000000000000000000000000000
000000000000000000000000000000000000000000000000000000000000000000000000000
0000000000000000000000000000000000000000000
````





That's one huge fraction!

## Other Compatible Number Types

#### BigFloats

````julia

prob_ode_biglinear = ODEProblem(f,big(1.0)/big(2.0),(big(0.0),big(1.0)),big(1.01))
sol =solve(prob_ode_biglinear,Tsit5())
println(sol[end])
````


````
1.3728004409033052741874280168680581033244287254502476809362138448867210367
01333
````





#### DoubleFloats.jl

There's are Float128-like types. Higher precision, but fixed and faster than arbitrary precision.

````julia

using DoubleFloats
prob_ode_doublelinear = ODEProblem(f,Double64(1)/Double64(2),(Double64(0),Double64(1)),Double64(1.01))
sol =solve(prob_ode_doublelinear,Tsit5())
println(sol[end])
````


````
1.3728004409033052
````





#### ArbFloats

These high precision numbers which are much faster than Bigs for less than 500-800 bits of accuracy.

````julia

using ArbNumerics
prob_ode_arbfloatlinear = ODEProblem(f,ArbFloat(1)/ArbFloat(2),(ArbFloat(0.0),ArbFloat(1.0)),ArbFloat(1.01))
sol =solve(prob_ode_arbfloatlinear,Tsit5())
println(sol[end])
````


````
1.372800440903305274187428016868
````





## Incompatible Number Systems

#### DecFP.jl

Next let's try DecFP. DecFP is a fixed-precision decimals library which is made to give both performance but known decimals of accuracy. Having already installed DecFP with `]add DecFP`, I can run the following:

````julia

using DecFP
prob_ode_decfplinear = ODEProblem(f,Dec128(1)/Dec128(2),(Dec128(0.0),Dec128(1.0)),Dec128(1.01))
sol =solve(prob_ode_decfplinear,Tsit5())
println(sol[end]); println(typeof(sol[end]))
````


````
Error: StackOverflowError:
````





#### Decimals.jl

Install with `]add Decimals`.

````julia

using Decimals
prob_ode_decimallinear = ODEProblem(f,[decimal("1.0")]./[decimal("2.0")],(0//1,1//1),decimal(1.01))
sol =solve(prob_ode_decimallinear,RK4(),dt=1/2^(6)) #Fails
println(sol[end]); println(typeof(sol[end]))
````


````
Error: MethodError: Decimals.Decimal(::Rational{Int64}) is ambiguous. Candi
dates:
  Decimals.Decimal(num::Real) in Decimals at /builds/JuliaGPU/DiffEqTutoria
ls.jl/.julia/packages/Decimals/Sb4j1/src/decimal.jl:13
  (::Type{T})(x::Rational{S}) where {S, T<:AbstractFloat} in Base at ration
al.jl:99
Possible fix, define
  Decimals.Decimal(::Rational{S}) where S
````





At the time of writing this, Decimals are not compatible. This is not on DifferentialEquations.jl's end, it's on partly on Decimal's end since it is not a subtype of Number. Thus it's not recommended you use Decimals with DifferentialEquations.jl

## Conclusion

As you can see, DifferentialEquations.jl can use arbitrary Julia-defined number systems in its arithmetic. If you need 128-bit floats, i.e. a bit more precision but not arbitrary, DoubleFloats.jl is a very good choice! For arbitrary precision, ArbNumerics are the most feature-complete and give great performance compared to BigFloats, and thus I recommend their use when high-precision (less than 512-800 bits) is required. DecFP is a great library for high-performance decimal numbers and works well as well. Other number systems could use some modernization.


## Appendix

 This tutorial is part of the SciMLTutorials.jl repository, found at: <https://github.com/SciML/SciMLTutorials.jl>.
 For more information on doing scientific machine learning (SciML) with open source software, check out <https://sciml.ai/>.

To locally run this tutorial, do the following commands:
```
using SciMLTutorials
SciMLTutorials.weave_file("type_handling","01-number_types.jmd")
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
Status `/builds/JuliaGPU/DiffEqTutorials.jl/tutorials/type_handling/Project.toml`
[7e558dbc-694d-5a72-987c-6f4ebed21442] ArbNumerics 1.2.1
[55939f99-70c6-5e9b-8bb0-5071ed7d61fd] DecFP 1.0.0
[abce61dc-4473-55a0-ba07-351d65e31d42] Decimals 0.4.1
[0c46a032-eb83-5123-abaf-570d42b7fbaa] DifferentialEquations 6.15.0
[497a8b3b-efae-58df-a0af-a86822472b78] DoubleFloats 1.1.13
[eff96d63-e80a-5855-80a2-b1b0885c5ab7] Measurements 2.3.0
[1dea7af3-3e70-54e6-95c3-0bf5283fa5ed] OrdinaryDiffEq 5.42.9
[91a5bcdd-55d7-5caf-9e0b-520d859cae80] Plots 1.6.6
[1986cc42-f94f-5a68-af5c-568840ba703d] Unitful 1.4.1
```
