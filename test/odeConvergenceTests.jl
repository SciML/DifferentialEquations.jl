using DifferentialEquations
prob = twoDimlinearODEExample()

## Convergence Testing
println("Convergence Test on Linear")
Δts = 1.//2.^(8:-1:4)

sim = testConvergence(Δts,prob,alg="Euler")
bool1 = abs(sim.𝒪est["final"]-1) < 0.1
sim2 = testConvergence(Δts,prob,alg="Midpoint")
bool2 = abs(sim2.𝒪est["l∞"]-2) < 0.1
sim3 = testConvergence(Δts,prob,alg="RK4")
bool3 = abs(sim3.𝒪est["l∞"]-4) < 0.1
tab = constructHuen()
sim4 = testConvergence(Δts,prob,alg="ExplicitRK",tableau=tab)
bool4 = abs(sim4.𝒪est["l∞"]-2) < 0.1

tab = constructRalston()
sim5 = testConvergence(Δts,prob,alg="ExplicitRK",tableau=tab)
bool5 = abs(sim5.𝒪est["l∞"]-2) < 0.1

tab = constructBogakiShampine()
sim6 = testConvergence(Δts,prob,alg="ExplicitRK",tableau=tab)
bool6 = abs(sim6.𝒪est["l∞"]-3) < 0.1

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

#Need to make larger or else below machine ϵ
Δts = 1./2.0.^(2:-1:-2)
tab = constructRKF8()
sim10 = testConvergence(Δts,prob,alg="ExplicitRK",tableau=tab)
bool10 = abs(sim10.𝒪est["l∞"]-8) < 0.1
bool1 && bool2 && bool3 && bool4 && bool5 && bool6 && bool7 && bool8 && bool9 && bool10
