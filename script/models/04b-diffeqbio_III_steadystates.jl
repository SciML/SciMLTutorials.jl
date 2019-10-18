
using DiffEqBiological, Plots
gr(); default(fmt = :png);


bistable_switch = @reaction_network begin
    d,    (X,Y) → ∅
    hillR(Y,v1,K1,n1), ∅ → X
    hillR(X,v2,K2,n2), ∅ → Y
end d v1 K1 n1 v2 K2 n2
d = 0.01;
v1 = 1.5; K1 = 30; n1 = 3;
v2 = 1.; K2 = 30; n2 = 3;
bistable_switch_p = [d, v1 ,K1, n1, v2, K2, n2];


ss = steady_states(bistable_switch, bistable_switch_p)


stability(ss,bistable_switch, bistable_switch_p)


rn1 = @reaction_network begin
    p, ∅ → X
    hill(X,v,K,n), X → ∅
end p v K n
p1 = [1.,2.5,1.5,1.5]
steady_states(rn1,p1)


rn2 = @reaction_network begin
    p, ∅ → X
    log(X), X → ∅
end p
p2 = [1.]
steady_states(rn2,p2)


bif = bifurcations(bistable_switch, bistable_switch_p, :v1, (.1,5.))
plot(bif,ylabel="[X]",label="")
plot!([[],[]],color=[:blue :red],label = ["Stable" "Unstable"])


plot(bif,2,ylabel="[Y]")
plot!([[],[]],color=[:blue :red],label = ["Stable" "Unstable"])


bif = bifurcations(bistable_switch, bistable_switch_p,:v1,(.1,10.))
plot(bif,linewidth=1.,title="A bifurcation diagram",ylabel="Steady State concentration")
plot!([[],[]],color=[:blue :red],label = ["Stable" "Unstable"])


bif = bifurcation_grid(bistable_switch, bistable_switch_p,:n1,1.:5.)
plot(bif)
scatter!([[],[]],color=[:blue :red],label = ["Stable" "Unstable"])


bif = bifurcation_grid_diagram(bistable_switch, bistable_switch_p,:n1,0.:4.,:v1,(.1,5.))
plot(bif)
plot!([[],[]],color=[:blue :red],label = ["Stable" "Unstable"])


bif = bifurcation_grid_2d(bistable_switch, bistable_switch_p,:n1,1.:3.,:n2,1.:10.)
plot(bif)
scatter!([[],[]],color=[:blue :red],label = ["Stable" "Unstable"])


brusselator = @reaction_network begin
    A, ∅ → X
    1, 2X + Y → 3X
    B, X → Y
    1, X → ∅
end A B;
A = 0.5; B = 4.;
brusselator_p = [A, B];


bif = bifurcations(brusselator,brusselator_p,:B,(0.1,2.5))
plot(bif,2)
plot!([[],[],[],[]],color=[:blue :cyan :orange :red],label = ["Stable Real" "Stable Complex" "Unstable Complex" "Unstable Real"])


bif = bifurcation_grid_diagram(brusselator,brusselator_p,:B,0.5:0.02:5.0,:A,(0.2,5.0))
plot(bif)
plot!([[],[],[],[]],color=[:blue :cyan :orange :red],label = ["Stable Real" "Stable Complex" "Unstable Complex" "Unstable Real"])


using DiffEqTutorials
DiffEqTutorials.tutorial_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file], remove_homedir=true)

