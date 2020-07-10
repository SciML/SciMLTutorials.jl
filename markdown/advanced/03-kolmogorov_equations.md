---
author: "Ashutosh Bharambe"
title: "Kolmogorov Backward Equations"
---
````julia
using Flux, StochasticDiffEq
using NeuralNetDiffEq
using Plots
using CuArrays
using CUDAnative
````




## Introduction on Backward Kolmogorov Equations

The backward Kolmogorov Equation deals with a terminal condtion.
The one dimensional backward kolmogorov equation that we are going to deal with is of the form :

$$
  \frac{\partial p}{\partial t} = -\mu(x)\frac{\partial p}{\partial x} - \frac{1}{2}{\sigma^2}(x)\frac{\partial^2 p}{\partial x^2} ,\hspace{0.5cm} p(T , x) = \varphi(x)
$$
for all $ t \in{ [0 , T] } $ and for all $ x \in R^d $

#### The Black Scholes Model

The Black-Scholes Model governs the price evolution of the European put or call option. In the below equation V is the price of some derivative , S is the Stock Price , r is the risk free interest
rate and σ the volatility of the stock returns. The payoff at a time T is known to us. And this makes it a terminal PDE. In case of an European put option the PDE is:
$$
  \frac{\partial V}{\partial t} + rS\frac{\partial V}{\partial S} + \frac{1}{2}{\sigma^2}{S^2}\frac{\partial^2 V}{\partial S^2} -rV = 0  ,\hspace{0.5cm} V(T , S) =  max\{\mathcal{K} - S , 0 \}
$$
for all $ t \in{ [0 , T] } $ and for all $ S \in R^d $

In order to make the above equation in the form of the Backward - Kolmogorov PDE we should substitute

$$
  V(S , t) = e^{r(t-T)}p(S , t)
$$
and thus we get
$$
  e^{r(t-T)}\frac{\partial p}{\partial t} + re^{r(t-T)}p(S , t)  = -\mu(x)\frac{\partial p}{\partial x}e^{r(t-T)} - \frac{1}{2}{\sigma^2}(x)\frac{\partial^2 p}{\partial x^2}e^{r(t-T)}
  + re^{r(t-T)}p(S , t)
$$
And the terminal condition
$$
  p(S , T) = max\{ \mathcal{K} - x , 0 \}
$$
We will train our model and the model itself will be the solution of the equation
## Defining the problem and the solver
We should start defining the terminal condition for our equation:
````julia
function phi(xi)
    y = Float64[]
    K = 100
    for x in eachcol(xi)
        val = max(K - maximum(x) , 0.00)
        y = push!(y , val)
    end
    y = reshape(y , 1 , size(y)[1] )
    return y
end
````


````
phi (generic function with 1 method)
````




Now we shall define the problem :
We will define the σ and μ by comparing it to the orignal equation. The xspan is the span of initial stock prices.
````julia
d = 1
r = 0.04
sigma = 0.2
xspan = (80.00 , 115.0)
tspan = (0.0 , 1.0)
σ(du , u , p , t) = du .= sigma.*u
μ(du , u , p , t) = du .= r.*u
prob = KolmogorovPDEProblem(μ , σ , phi , xspan , tspan, d)
````


````
KolmogorovPDEProblem
timespan: (0.0, 1.0)xspan: (80.0, 115.0)μ
Main.##WeaveSandBox#583.μSigma
Main.##WeaveSandBox#583.σ
````




Now once we have defined our problem it is necessary to define the parameters for the solver.
````julia
sdealg = EM()
ensemblealg = EnsembleThreads()
dt = 0.01
dx = 0.01
trajectories = 100000
````


````
100000
````





Now lets define our model m and the optimiser
````julia
m = Chain(Dense(d, 64, elu),Dense(64, 128, elu),Dense(128 , 16 , elu) , Dense(16 , 1))
use_gpu = false
if CUDAnative.functional() == true
  m = fmap(CuArrays.cu , m)
  use_gpu = true
end
opt = Flux.ADAM(0.0005)
````


````
Flux.Optimise.ADAM(0.0005, (0.9, 0.999), IdDict{Any,Any}())
````




And then finally call the solver
````julia
@time sol = solve(prob, NeuralNetDiffEq.NNKolmogorov(m, opt, sdealg, ensemblealg), verbose = true, dt = dt,
            dx = dx , trajectories = trajectories , abstol=1e-6, maxiters = 1000 , use_gpu = use_gpu)
````


````
Current loss is: 134.80429
Current loss is: 132.68681
Current loss is: 135.3199
Current loss is: 137.07219
Current loss is: 136.294
Current loss is: 134.30772
Current loss is: 132.60184
Current loss is: 131.96268
Current loss is: 132.3382
Current loss is: 133.06836
Current loss is: 133.43152
Current loss is: 133.12303
Current loss is: 132.3151
Current loss is: 131.4221
Current loss is: 130.83995
Current loss is: 130.79732
Current loss is: 131.25005
Current loss is: 131.46846
Current loss is: 131.15599
Current loss is: 130.4978
Current loss is: 129.96117
Current loss is: 129.75539
Current loss is: 129.71904
Current loss is: 129.63995
Current loss is: 129.3646
Current loss is: 128.98088
Current loss is: 128.70297
Current loss is: 128.67253
Current loss is: 128.56552
Current loss is: 128.2728
Current loss is: 127.932205
Current loss is: 127.685234
Current loss is: 127.514496
Current loss is: 127.33897
Current loss is: 127.106346
Current loss is: 126.884125
Current loss is: 126.73983
Current loss is: 126.590836
Current loss is: 126.38942
Current loss is: 126.20772
Current loss is: 126.069145
Current loss is: 125.926544
Current loss is: 125.74502
Current loss is: 125.56654
Current loss is: 125.417435
Current loss is: 125.261314
Current loss is: 125.075195
Current loss is: 124.8964
Current loss is: 124.730034
Current loss is: 124.5438
Current loss is: 124.339836
Current loss is: 124.14451
Current loss is: 123.94598
Current loss is: 123.72298
Current loss is: 123.520996
Current loss is: 123.333206
Current loss is: 123.139755
Current loss is: 122.92544
Current loss is: 122.70611
Current loss is: 122.48011
Current loss is: 122.226814
Current loss is: 121.971535
Current loss is: 121.75648
Current loss is: 121.552986
Current loss is: 121.36565
Current loss is: 121.182846
Current loss is: 121.00063
Current loss is: 120.81369
Current loss is: 120.62416
Current loss is: 120.430595
Current loss is: 120.23543
Current loss is: 120.04179
Current loss is: 119.84623
Current loss is: 119.64889
Current loss is: 119.45074
Current loss is: 119.25155
Current loss is: 119.0524
Current loss is: 118.851326
Current loss is: 118.64724
Current loss is: 118.4423
Current loss is: 118.236824
Current loss is: 118.03023
Current loss is: 117.82283
Current loss is: 117.61483
Current loss is: 117.40691
Current loss is: 117.19704
Current loss is: 116.98524
Current loss is: 116.772896
Current loss is: 116.56058
Current loss is: 116.3482
Current loss is: 116.13504
Current loss is: 115.92165
Current loss is: 115.7075
Current loss is: 115.49324
Current loss is: 115.278595
Current loss is: 115.06389
Current loss is: 114.84874
Current loss is: 114.63367
Current loss is: 114.41883
Current loss is: 114.204384
Current loss is: 113.98986
Current loss is: 113.775795
Current loss is: 113.56257
Current loss is: 113.35043
Current loss is: 113.13885
Current loss is: 112.92783
Current loss is: 112.7173
Current loss is: 112.50773
Current loss is: 112.299484
Current loss is: 112.09243
Current loss is: 111.88636
Current loss is: 111.68146
Current loss is: 111.47822
Current loss is: 111.27672
Current loss is: 111.07702
Current loss is: 110.87928
Current loss is: 110.68345
Current loss is: 110.489494
Current loss is: 110.29752
Current loss is: 110.10763
Current loss is: 109.920044
Current loss is: 109.734856
Current loss is: 109.55201
Current loss is: 109.37162
Current loss is: 109.19388
Current loss is: 109.018906
Current loss is: 108.84678
Current loss is: 108.67749
Current loss is: 108.511116
Current loss is: 108.347855
Current loss is: 108.1879
Current loss is: 108.0312
Current loss is: 107.87794
Current loss is: 107.72861
Current loss is: 107.58311
Current loss is: 107.44066
Current loss is: 107.30092
Current loss is: 107.164536
Current loss is: 107.03286
Current loss is: 106.908615
Current loss is: 106.801956
Current loss is: 106.73256
Current loss is: 106.69902
Current loss is: 106.55405
Current loss is: 106.32728
Current loss is: 106.26285
Current loss is: 106.22597
Current loss is: 106.04488
Current loss is: 105.94188
Current loss is: 105.91665
Current loss is: 105.7794
Current loss is: 105.6735
Current loss is: 105.64521
Current loss is: 105.54008
Current loss is: 105.44043
Current loss is: 105.40802
Current loss is: 105.32928
Current loss is: 105.238
Current loss is: 105.20087
Current loss is: 105.1448
Current loss is: 105.06444
Current loss is: 105.020996
Current loss is: 104.9824
Current loss is: 104.91737
Current loss is: 104.86755
Current loss is: 104.83761
Current loss is: 104.791374
Current loss is: 104.740135
Current loss is: 104.70874
Current loss is: 104.67808
Current loss is: 104.63434
Current loss is: 104.59739
Current loss is: 104.57168
Current loss is: 104.53992
Current loss is: 104.502625
Current loss is: 104.473206
Current loss is: 104.44804
Current loss is: 104.417
Current loss is: 104.38543
Current loss is: 104.36066
Current loss is: 104.336845
Current loss is: 104.309845
Current loss is: 104.28572
Current loss is: 104.26604
Current loss is: 104.24493
Current loss is: 104.22253
Current loss is: 104.202995
Current loss is: 104.18428
Current loss is: 104.1636
Current loss is: 104.14368
Current loss is: 104.12613
Current loss is: 104.10859
Current loss is: 104.09064
Current loss is: 104.07438
Current loss is: 104.05957
Current loss is: 104.04452
Current loss is: 104.029434
Current loss is: 104.01552
Current loss is: 104.002335
Current loss is: 103.988884
Current loss is: 103.97555
Current loss is: 103.96309
Current loss is: 103.95121
Current loss is: 103.939354
Current loss is: 103.9277
Current loss is: 103.916756
Current loss is: 103.90632
Current loss is: 103.895935
Current loss is: 103.88572
Current loss is: 103.87595
Current loss is: 103.86663
Current loss is: 103.857414
Current loss is: 103.84836
Current loss is: 103.83966
Current loss is: 103.83134
Current loss is: 103.82323
Current loss is: 103.81524
Current loss is: 103.80745
Current loss is: 103.79995
Current loss is: 103.792656
Current loss is: 103.7855
Current loss is: 103.77846
Current loss is: 103.771614
Current loss is: 103.765015
Current loss is: 103.75859
Current loss is: 103.75229
Current loss is: 103.74612
Current loss is: 103.74012
Current loss is: 103.7343
Current loss is: 103.728615
Current loss is: 103.723045
Current loss is: 103.71759
Current loss is: 103.712296
Current loss is: 103.70713
Current loss is: 103.70212
Current loss is: 103.697174
Current loss is: 103.692375
Current loss is: 103.6877
Current loss is: 103.683105
Current loss is: 103.67866
Current loss is: 103.67428
Current loss is: 103.670006
Current loss is: 103.665825
Current loss is: 103.66175
Current loss is: 103.657745
Current loss is: 103.65385
Current loss is: 103.65004
Current loss is: 103.64631
Current loss is: 103.64264
Current loss is: 103.63907
Current loss is: 103.63556
Current loss is: 103.63214
Current loss is: 103.628784
Current loss is: 103.625496
Current loss is: 103.6223
Current loss is: 103.61915
Current loss is: 103.616066
Current loss is: 103.61306
Current loss is: 103.610115
Current loss is: 103.607216
Current loss is: 103.604385
Current loss is: 103.60161
Current loss is: 103.598885
Current loss is: 103.59623
Current loss is: 103.59362
Current loss is: 103.59106
Current loss is: 103.58855
Current loss is: 103.58609
Current loss is: 103.58369
Current loss is: 103.581314
Current loss is: 103.57901
Current loss is: 103.576744
Current loss is: 103.574524
Current loss is: 103.572334
Current loss is: 103.570206
Current loss is: 103.5681
Current loss is: 103.566055
Current loss is: 103.56404
Current loss is: 103.562065
Current loss is: 103.56012
Current loss is: 103.55822
Current loss is: 103.55636
Current loss is: 103.554504
Current loss is: 103.55272
Current loss is: 103.55095
Current loss is: 103.54922
Current loss is: 103.54754
Current loss is: 103.545876
Current loss is: 103.544235
Current loss is: 103.542656
Current loss is: 103.54113
Current loss is: 103.539635
Current loss is: 103.53827
Current loss is: 103.53701
Current loss is: 103.53602
Current loss is: 103.5355
Current loss is: 103.535904
Current loss is: 103.53801
Current loss is: 103.543785
Current loss is: 103.5554
Current loss is: 103.57909
Current loss is: 103.610985
Current loss is: 103.65288
Current loss is: 103.648766
Current loss is: 103.60535
Current loss is: 103.53725
Current loss is: 103.5234
Current loss is: 103.56164
Current loss is: 103.58419
Current loss is: 103.561844
Current loss is: 103.52108
Current loss is: 103.520256
Current loss is: 103.54771
Current loss is: 103.55117
Current loss is: 103.52707
Current loss is: 103.51009
Current loss is: 103.521355
Current loss is: 103.53498
Current loss is: 103.52456
Current loss is: 103.50821
Current loss is: 103.50875
Current loss is: 103.519035
Current loss is: 103.51868
Current loss is: 103.506966
Current loss is: 103.502335
Current loss is: 103.50819
Current loss is: 103.51102
Current loss is: 103.505165
Current loss is: 103.49897
Current loss is: 103.50019
Current loss is: 103.50383
Current loss is: 103.50214
Current loss is: 103.49712
Current loss is: 103.49517
Current loss is: 103.497116
Current loss is: 103.49804
Current loss is: 103.49537
Current loss is: 103.49238
Current loss is: 103.49212
Current loss is: 103.49322
Current loss is: 103.49282
Current loss is: 103.49061
Current loss is: 103.488914
Current loss is: 103.48891
Current loss is: 103.4893
Current loss is: 103.48861
Current loss is: 103.487015
Current loss is: 103.48588
Current loss is: 103.48572
Current loss is: 103.48572
Current loss is: 103.48513
Current loss is: 103.483986
Current loss is: 103.48308
Current loss is: 103.48273
Current loss is: 103.48257
Current loss is: 103.48212
Current loss is: 103.4813
Current loss is: 103.480484
Current loss is: 103.479996
Current loss is: 103.479706
Current loss is: 103.47937
Current loss is: 103.478806
Current loss is: 103.478134
Current loss is: 103.47757
Current loss is: 103.477165
Current loss is: 103.47683
Current loss is: 103.476425
Current loss is: 103.47594
Current loss is: 103.47539
Current loss is: 103.474915
Current loss is: 103.474525
Current loss is: 103.474174
Current loss is: 103.4738
Current loss is: 103.47338
Current loss is: 103.472916
Current loss is: 103.47249
Current loss is: 103.47211
Current loss is: 103.471756
Current loss is: 103.471405
Current loss is: 103.47104
Current loss is: 103.47066
Current loss is: 103.47027
Current loss is: 103.469894
Current loss is: 103.46955
Current loss is: 103.46921
Current loss is: 103.468895
Current loss is: 103.46856
Current loss is: 103.46822
Current loss is: 103.46788
Current loss is: 103.467545
Current loss is: 103.46722
Current loss is: 103.4669
Current loss is: 103.4666
Current loss is: 103.4663
Current loss is: 103.46601
Current loss is: 103.465706
Current loss is: 103.46541
Current loss is: 103.46511
Current loss is: 103.46482
Current loss is: 103.464516
Current loss is: 103.46424
Current loss is: 103.46396
Current loss is: 103.4637
Current loss is: 103.46342
Current loss is: 103.46316
Current loss is: 103.46289
Current loss is: 103.46264
Current loss is: 103.46238
Current loss is: 103.46212
Current loss is: 103.46187
Current loss is: 103.46162
Current loss is: 103.461365
Current loss is: 103.46112
Current loss is: 103.460884
Current loss is: 103.46064
Current loss is: 103.46041
Current loss is: 103.460175
Current loss is: 103.459946
Current loss is: 103.45972
Current loss is: 103.459496
Current loss is: 103.459274
Current loss is: 103.459045
Current loss is: 103.45884
Current loss is: 103.45861
Current loss is: 103.458405
Current loss is: 103.4582
Current loss is: 103.45798
Current loss is: 103.45778
Current loss is: 103.45758
Current loss is: 103.457375
Current loss is: 103.457184
Current loss is: 103.457
Current loss is: 103.45681
Current loss is: 103.45665
Current loss is: 103.45649
Current loss is: 103.456375
Current loss is: 103.4563
Current loss is: 103.456345
Current loss is: 103.4566
Current loss is: 103.45721
Current loss is: 103.45848
Current loss is: 103.46117
Current loss is: 103.46625
Current loss is: 103.476555
Current loss is: 103.4943
Current loss is: 103.52849
Current loss is: 103.572174
Current loss is: 103.63482
Current loss is: 103.64093
Current loss is: 103.60299
Current loss is: 103.50194
Current loss is: 103.45424
Current loss is: 103.48843
Current loss is: 103.539986
Current loss is: 103.54612
Current loss is: 103.488914
Current loss is: 103.453094
Current loss is: 103.47266
Current loss is: 103.50531
Current loss is: 103.50383
Current loss is: 103.46754
Current loss is: 103.45237
Current loss is: 103.471214
Current loss is: 103.486
Current loss is: 103.47446
Current loss is: 103.45409
Current loss is: 103.455246
Current loss is: 103.46999
Current loss is: 103.47134
Current loss is: 103.458786
Current loss is: 103.45082
Current loss is: 103.45688
Current loss is: 103.464424
Current loss is: 103.46062
Current loss is: 103.45225
Current loss is: 103.45083
Current loss is: 103.45613
Current loss is: 103.45858
Current loss is: 103.45418
Current loss is: 103.44989
Current loss is: 103.450935
Current loss is: 103.4542
Current loss is: 103.4543
Current loss is: 103.45105
Current loss is: 103.44907
Current loss is: 103.45029
Current loss is: 103.45202
Current loss is: 103.45156
Current loss is: 103.449486
Current loss is: 103.44851
Current loss is: 103.449356
Current loss is: 103.45028
Current loss is: 103.449844
Current loss is: 103.448555
Current loss is: 103.44797
Current loss is: 103.44841
Current loss is: 103.448944
Current loss is: 103.448715
Current loss is: 103.447945
Current loss is: 103.44745
Current loss is: 103.44757
Current loss is: 103.44788
Current loss is: 103.44784
Current loss is: 103.447395
Current loss is: 103.44697
Current loss is: 103.4469
Current loss is: 103.44704
Current loss is: 103.447075
Current loss is: 103.44688
Current loss is: 103.446556
Current loss is: 103.44635
Current loss is: 103.446335
Current loss is: 103.44639
Current loss is: 103.446335
Current loss is: 103.44614
Current loss is: 103.44593
Current loss is: 103.445816
Current loss is: 103.44578
Current loss is: 103.445755
Current loss is: 103.44568
Current loss is: 103.44553
Current loss is: 103.44538
Current loss is: 103.445274
Current loss is: 103.44522
Current loss is: 103.44516
Current loss is: 103.4451
Current loss is: 103.445
Current loss is: 103.44487
Current loss is: 103.44478
Current loss is: 103.444695
Current loss is: 103.444626
Current loss is: 103.44456
Current loss is: 103.44448
Current loss is: 103.44439
Current loss is: 103.4443
Current loss is: 103.4442
Current loss is: 103.44413
Current loss is: 103.444046
Current loss is: 103.44398
Current loss is: 103.44391
Current loss is: 103.443825
Current loss is: 103.44375
Current loss is: 103.443665
Current loss is: 103.44359
Current loss is: 103.443504
Current loss is: 103.443436
Current loss is: 103.443375
Current loss is: 103.44329
Current loss is: 103.443214
Current loss is: 103.443146
Current loss is: 103.443085
Current loss is: 103.443
Current loss is: 103.442924
Current loss is: 103.442856
Current loss is: 103.44278
Current loss is: 103.44272
Current loss is: 103.442635
Current loss is: 103.44258
Current loss is: 103.4425
Current loss is: 103.44244
Current loss is: 103.44237
Current loss is: 103.44228
Current loss is: 103.44223
Current loss is: 103.442154
Current loss is: 103.442085
Current loss is: 103.44202
Current loss is: 103.441956
Current loss is: 103.44189
Current loss is: 103.441826
Current loss is: 103.44175
Current loss is: 103.44168
Current loss is: 103.44162
Current loss is: 103.441536
Current loss is: 103.44149
Current loss is: 103.441414
Current loss is: 103.44134
Current loss is: 103.441284
Current loss is: 103.44123
Current loss is: 103.441154
Current loss is: 103.441086
Current loss is: 103.44102
Current loss is: 103.44096
Current loss is: 103.440895
Current loss is: 103.44083
Current loss is: 103.44078
Current loss is: 103.44072
Current loss is: 103.440636
Current loss is: 103.440575
Current loss is: 103.44053
Current loss is: 103.44047
Current loss is: 103.4404
Current loss is: 103.44034
Current loss is: 103.44027
Current loss is: 103.440216
Current loss is: 103.44017
Current loss is: 103.440094
Current loss is: 103.44005
Current loss is: 103.43998
Current loss is: 103.439926
Current loss is: 103.43988
Current loss is: 103.43983
Current loss is: 103.43979
Current loss is: 103.43976
Current loss is: 103.439766
Current loss is: 103.43982
Current loss is: 103.43996
Current loss is: 103.44028
Current loss is: 103.44094
Current loss is: 103.44223
Current loss is: 103.44484
Current loss is: 103.44976
Current loss is: 103.459785
Current loss is: 103.47774
Current loss is: 103.51408
Current loss is: 103.568
Current loss is: 103.662094
Current loss is: 103.7201
Current loss is: 103.74898
Current loss is: 103.608154
Current loss is: 103.472244
Current loss is: 103.44884
Current loss is: 103.53315
Current loss is: 103.60432
Current loss is: 103.5456
Current loss is: 103.459496
Current loss is: 103.444885
Current loss is: 103.50066
Current loss is: 103.535
Current loss is: 103.48749
Current loss is: 103.44096
Current loss is: 103.4545
Current loss is: 103.4893
Current loss is: 103.48536
Current loss is: 103.448326
Current loss is: 103.44094
Current loss is: 103.46509
Current loss is: 103.4725
Current loss is: 103.453354
Current loss is: 103.438126
Current loss is: 103.447845
Current loss is: 103.460945
Current loss is: 103.453545
Current loss is: 103.43982
Current loss is: 103.44017
Current loss is: 103.449844
Current loss is: 103.45099
Current loss is: 103.44179
Current loss is: 103.437675
Current loss is: 103.44292
Current loss is: 103.44671
Current loss is: 103.44301
Current loss is: 103.437836
Current loss is: 103.43866
Current loss is: 103.4425
Current loss is: 103.44243
Current loss is: 103.43889
Current loss is: 103.437195
Current loss is: 103.439095
Current loss is: 103.440796
Current loss is: 103.43947
Current loss is: 103.43732
Current loss is: 103.437256
Current loss is: 103.4387
Current loss is: 103.43915
Current loss is: 103.437904
Current loss is: 103.43682
Current loss is: 103.43716
Current loss is: 103.43798
Current loss is: 103.437965
Current loss is: 103.43711
Current loss is: 103.4366
Current loss is: 103.4369
Current loss is: 103.43735
Current loss is: 103.437225
Current loss is: 103.436676
Current loss is: 103.43639
Current loss is: 103.436584
Current loss is: 103.43684
Current loss is: 103.43676
Current loss is: 103.43641
Current loss is: 103.43622
Current loss is: 103.436295
Current loss is: 103.43645
Current loss is: 103.43641
Current loss is: 103.436195
Current loss is: 103.43605
Current loss is: 103.43605
Current loss is: 103.43613
Current loss is: 103.43613
Current loss is: 103.436005
Current loss is: 103.435875
Current loss is: 103.43584
Current loss is: 103.43586
Current loss is: 103.43587
Current loss is: 103.435814
Current loss is: 103.43573
Current loss is: 103.43567
Current loss is: 103.43564
Current loss is: 103.435646
Current loss is: 103.43564
Current loss is: 103.43558
Current loss is: 103.43552
Current loss is: 103.43548
Current loss is: 103.43545
Current loss is: 103.43545
Current loss is: 103.43541
Current loss is: 103.435356
Current loss is: 103.43532
Current loss is: 103.43528
Current loss is: 103.43526
Current loss is: 103.435234
Current loss is: 103.43522
Current loss is: 103.435165
Current loss is: 103.435135
Current loss is: 103.4351
Current loss is: 103.43507
Current loss is: 103.43505
Current loss is: 103.43502
Current loss is: 103.43499
Current loss is: 103.43497
Current loss is: 103.43493
Current loss is: 103.434906
Current loss is: 103.43487
Current loss is: 103.43484
Current loss is: 103.434814
Current loss is: 103.434784
Current loss is: 103.43476
Current loss is: 103.43473
Current loss is: 103.43471
Current loss is: 103.43467
Current loss is: 103.43464
Current loss is: 103.434616
Current loss is: 103.434586
Current loss is: 103.43457
Current loss is: 103.43454
Current loss is: 103.434494
Current loss is: 103.43447
Current loss is: 103.43445
Current loss is: 103.43442
Current loss is: 103.434395
Current loss is: 103.434364
Current loss is: 103.43433
Current loss is: 103.43432
Current loss is: 103.43429
Current loss is: 103.43426
Current loss is: 103.43425
Current loss is: 103.43422
Current loss is: 103.43419
Current loss is: 103.434166
Current loss is: 103.434135
Current loss is: 103.43412
Current loss is: 103.434074
Current loss is: 103.43407
Current loss is: 103.43404
Current loss is: 103.43402
Current loss is: 103.433975
Current loss is: 103.43397
Current loss is: 103.43394
Current loss is: 103.433914
Current loss is: 103.4339
Current loss is: 103.43387
Current loss is: 103.433846
Current loss is: 103.433815
Current loss is: 103.4338
Current loss is: 103.43378
Current loss is: 103.433754
Current loss is: 103.43373
Current loss is: 103.43371
Current loss is: 103.43367
Current loss is: 103.43364
Current loss is: 103.433624
Current loss is: 103.43361
Current loss is: 103.433586
Current loss is: 103.433556
Current loss is: 103.43354
Current loss is: 103.43351
Current loss is: 103.433495
Current loss is: 103.433464
Current loss is: 103.43345
Current loss is: 103.433426
Current loss is: 103.43341
Current loss is: 103.43338
Current loss is: 103.43335
Current loss is: 103.433334
Current loss is: 103.43332
Current loss is: 103.4333
Current loss is: 103.433266
Current loss is: 103.43325
Current loss is: 103.43323
Current loss is: 103.433205
Current loss is: 103.43319
Current loss is: 103.43317
Current loss is: 103.433136
Current loss is: 103.43312
Current loss is: 103.4331
Current loss is: 103.433075
Current loss is: 103.43306
Current loss is: 103.43304
Current loss is: 103.433014
Current loss is: 103.43299
Current loss is: 103.432976
Current loss is: 103.43296
Current loss is: 103.43294
Current loss is: 103.432915
Current loss is: 103.432884
Current loss is: 103.43287
Current loss is: 103.432846
Current loss is: 103.43284
Current loss is: 103.43281
Current loss is: 103.43281
Current loss is: 103.43278
Current loss is: 103.43277
Current loss is: 103.43275
Current loss is: 103.43274
Current loss is: 103.43274
Current loss is: 103.43278
Current loss is: 103.43284
Current loss is: 103.432976
Current loss is: 103.43326
Current loss is: 103.433815
Current loss is: 103.434944
Current loss is: 103.437126
Current loss is: 103.44159
Current loss is: 103.450096
Current loss is: 103.46775
Current loss is: 103.499306
Current loss is: 103.56385
Current loss is: 103.653496
Current loss is: 103.804985
Current loss is: 103.85194
Current loss is: 103.82547
Current loss is: 103.58208
Current loss is: 103.4388
Current loss is: 103.50228
Current loss is: 103.629974
Current loss is: 103.65676
Current loss is: 103.5153
Current loss is: 103.43291
Current loss is: 103.488914
Current loss is: 103.5605
Current loss is: 103.53976
Current loss is: 103.45113
Current loss is: 103.442894
Current loss is: 103.50251
Current loss is: 103.506134
Current loss is: 103.45501
Current loss is: 103.433815
Current loss is: 103.467896
Current loss is: 103.48778
Current loss is: 103.45598
Current loss is: 103.43242
Current loss is: 103.448845
Current loss is: 103.46609
Current loss is: 103.45338
Current loss is: 103.43339
Current loss is: 103.4394
Current loss is: 103.45417
Current loss is: 103.44815
Current loss is: 103.43404
Current loss is: 103.43477
Current loss is: 103.44482
Current loss is: 103.44488
Current loss is: 103.435135
Current loss is: 103.43252
Current loss is: 103.438774
Current loss is: 103.44092
Current loss is: 103.435555
Current loss is: 103.431915
Current loss is: 103.434975
Current loss is: 103.437965
Current loss is: 103.43539
Current loss is: 103.432014
Current loss is: 103.432846
Current loss is: 103.43536
Current loss is: 103.435
Current loss is: 103.432465
Current loss is: 103.43184
Current loss is: 103.433426
Current loss is: 103.4341
Current loss is: 103.432785
Current loss is: 103.431656
Current loss is: 103.43222
Current loss is: 103.433136
Current loss is: 103.432785
Current loss is: 103.4318
Current loss is: 103.431625
Current loss is: 103.43225
Current loss is: 103.43251
Current loss is: 103.431984
Current loss is: 103.43151
Current loss is: 103.431656
Current loss is: 103.43203
Current loss is: 103.431984
Current loss is: 103.431595
Current loss is: 103.43142
Current loss is: 103.43161
Current loss is: 103.43178
Current loss is: 103.43165
Current loss is: 103.431404
Current loss is: 103.43136
Current loss is: 103.43149
Current loss is: 103.43155
Current loss is: 103.431435
Current loss is: 103.43129
Current loss is: 103.4313
Current loss is: 103.431366
Current loss is: 103.43139
Current loss is: 103.431305
Current loss is: 103.43122
Current loss is: 103.43122
Current loss is: 103.43126
Current loss is: 103.431244
Current loss is: 103.4312
Current loss is: 103.431145
Current loss is: 103.431145
Current loss is: 103.43117
Current loss is: 103.431145
Current loss is: 103.43113
Current loss is: 103.43108
Current loss is: 103.43107
Current loss is: 103.43107
Current loss is: 103.43107
Current loss is: 103.431046
Current loss is: 103.431015
Current loss is: 103.430984
Current loss is: 103.431
Current loss is: 103.431
Current loss is: 103.43098
Current loss is: 103.430954
Current loss is: 103.43093
Current loss is: 103.43093
Current loss is: 103.430916
Current loss is: 103.43091
Current loss is: 103.430885
Current loss is: 103.43087
Current loss is: 103.430855
Current loss is: 103.43085
Current loss is: 103.43084
Current loss is: 103.43082
Current loss is: 103.43081
Current loss is: 103.43079
Current loss is: 103.43078
Current loss is: 103.43078
Current loss is: 103.430756
Current loss is: 103.430756
Current loss is: 103.43074
Current loss is: 103.43072
Current loss is: 103.43071
Current loss is: 103.430695
Current loss is: 103.43069
Current loss is: 103.43068
Current loss is: 103.43066
Current loss is: 103.43064
Current loss is: 103.43064
Current loss is: 103.43062
Current loss is: 103.43061
Current loss is: 103.430595
Current loss is: 103.43058
Current loss is: 103.43058
Current loss is: 103.43056
Current loss is: 103.43055
Current loss is: 103.43053
Current loss is: 103.430534
Current loss is: 103.43052
Current loss is: 103.430504
Current loss is: 103.4305
Current loss is: 103.43048
Current loss is: 103.43046
Current loss is: 103.43045
Current loss is: 103.430435
Current loss is: 103.43042
Current loss is: 103.43042
Current loss is: 103.430405
Current loss is: 103.430405
Current loss is: 103.430374
Current loss is: 103.43037
Current loss is: 103.43036
Current loss is: 103.43034
Current loss is: 103.430336
Current loss is: 103.43032
Current loss is: 103.43032
Current loss is: 103.4303
Current loss is: 103.43029
Current loss is: 103.430275
Current loss is: 103.43027
Current loss is: 103.43026
Current loss is: 103.430244
Current loss is: 103.43023
Current loss is: 103.43023
Current loss is: 103.43021
Current loss is: 103.43021
Current loss is: 103.43019
Current loss is: 103.430176
Current loss is: 103.43016
Current loss is: 103.430145
Current loss is: 103.43014
Current loss is: 103.43014
Current loss is: 103.430115
Current loss is: 103.430115
Current loss is: 103.4301
 33.133163 seconds (120.32 M allocations: 122.216 GiB, 17.01% gc time)
(Float32[100.14 109.89 … 91.93 86.57], Float32[6.1517444 3.207456 … 10.0985
22 13.538801])
````




## Analyzing the solution
Now let us find a Monte-Carlo Solution and plot the both:
````julia
monte_carlo_sol = []
x_out = collect(85:2.00:110.00)
for x in x_out
  u₀= [x]
  g_val(du , u , p , t) = du .= 0.2.*u
  f_val(du , u , p , t) = du .= 0.04.*u
  dt = 0.01
  tspan = (0.0,1.0)
  prob = SDEProblem(f_val,g_val,u₀,tspan)
  output_func(sol,i) = (sol[end],false)
  ensembleprob_val = EnsembleProblem(prob , output_func = output_func )
  sim_val = solve(ensembleprob_val, EM(), EnsembleThreads() , dt=0.01, trajectories=100000,adaptive=false)
  s = reduce(hcat , sim_val.u)
  mean_phi = sum(phi(s))/length(phi(s))
  global monte_carlo_sol = push!(monte_carlo_sol , mean_phi)
end
````





##Plotting the Solutions
We should reshape the inputs and outputs to make it compatible with our model. This is the most important part. The algorithm gives a distributed function over all initial prices in the xspan.
````julia
x_model = reshape(x_out, 1 , size(x_out)[1])
if use_gpu == true
  m = fmap(cpu , m)
end
y_out = m(x_model)
y_out = reshape(y_out , 13 , 1)
````


````
13×1 Array{Float32,2}:
 14.650005
 13.240433
 11.890396
 10.645651
  9.488981
  8.407096
  7.4523544
  6.5946817
  5.8396764
  5.169847
  4.5574126
  3.9850886
  3.441834
````




And now finally we can plot the solutions
````julia
plot(x_out , y_out , lw = 3 ,  xaxis="Initial Stock Price", yaxis="Payoff" , label = "NNKolmogorov")
plot!(x_out , monte_carlo_sol , lw = 3 ,  xaxis="Initial Stock Price", yaxis="Payoff" ,label = "Monte Carlo Solutions")
````


![](figures/03-kolmogorov_equations_9_1.png)
