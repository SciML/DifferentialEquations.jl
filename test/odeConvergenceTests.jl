using DifferentialEquations
prob = twoDimlinearODEExample()

## Convergence Testing
println("Convergence Test on Linear")
Δts = 1.//2.^(8:-1:4)

println("Special RKs")
sim = testConvergence(Δts,prob,alg="Euler")
bool1 = abs(sim.𝒪est["final"]-1) < 0.1
sim2 = testConvergence(Δts,prob,alg="Midpoint")
bool2 = abs(sim2.𝒪est["l∞"]-2) < 0.1
sim3 = testConvergence(Δts,prob,alg="RK4")
bool3 = abs(sim3.𝒪est["l∞"]-4) < 0.1
tab = constructHuen()
sim4 = testConvergence(Δts,prob,alg="ExplicitRK",tableau=tab)
bool4 = abs(sim4.𝒪est["l∞"]-2) < 0.1
println("Other RKs")
tab = constructRalston()
sim5 = testConvergence(Δts,prob,alg="ExplicitRK",tableau=tab)
bool5 = abs(sim5.𝒪est["l∞"]-2) < 0.1

tab = constructBogakiShampine()
sim6 = testConvergence(Δts,prob,alg="ExplicitRK",tableau=tab)
bool6 = abs(sim6.𝒪est["l∞"]-3) < 0.1

println("Higher Order")
Δts = 1.//2.^(7:-1:3)
tab = constructRKF()
sim7 = testConvergence(Δts,prob,alg="ExplicitRK",tableau=tab)
bool7 = abs(sim7.𝒪est["l∞"]-5) < 0.1

tab = constructDormandPrince()
sim8 = testConvergence(Δts,prob,alg="ExplicitRK",tableau=tab)
bool8 = abs(sim8.𝒪est["l∞"]-5) < 0.1

tab = constructCashKarp()
sim9 = testConvergence(Δts,prob,alg="ExplicitRK",tableau=tab)
bool9 = abs(sim9.𝒪est["l∞"]-5) < 0.1

println("Super High Order")
#Need to make larger or else below machine ϵ
Δts = 1./2.0.^(2:-1:-2)
tab = constructRKF8()
sim10 = testConvergence(Δts,prob,alg="ExplicitRK",tableau=tab)
bool10 = abs(sim10.𝒪est["l∞"]-8) < 0.1
superbool1 = bool1 && bool2 && bool3 && bool4 && bool5 && bool6 && bool7 && bool8 && bool9 && bool10

println("Tests pass: $superbool1")
### Stiff Solvers

println("Convergence Test on Stiff")
Δts = 1.//2.^(8:-1:4)

sim11 = testConvergence(Δts,prob,alg="ImplicitEuler")
bool11 = abs(sim11.𝒪est["final"]-1) < 0.1

sim12 = testConvergence(Δts,prob,alg="Trapezoid")
bool12 = abs(sim12.𝒪est["final"]-2) < 0.1

superbool2 = bool11 && bool12
println("Tests pass: $superbool2")
superbool1 && superbool2
