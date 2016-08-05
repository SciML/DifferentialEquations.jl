# This definitely needs cleaning
using DifferentialEquations
probArr = Vector{DEProblem}(2)
probArr[1] = linearODEExample()
probArr[2] = twoDimlinearODEExample()
srand(100)
## Convergence Testing
println("Convergence Test on Linear")
Δts = 1.//2.^(8:-1:4)
testTol = 0.2
superduperbool = Vector{Bool}(2)

for i = 1:2
  prob = probArr[i]
  println("Special RKs")
  sim = test_convergence(Δts,prob,alg=:Euler)
  bool1 = abs(sim.𝒪est[:final]-1) < testTol
  sim2 = test_convergence(Δts,prob,alg=:Midpoint)
  bool2 = abs(sim2.𝒪est[:l∞]-2) < testTol
  sim3 = test_convergence(Δts,prob,alg=:RK4)
  bool3 = abs(sim3.𝒪est[:l∞]-4) < testTol
  tab = constructHuen()
  sim4 = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
  bool4 = abs(sim4.𝒪est[:l∞]-2) < testTol
  println("Other RKs")
  tab = constructRalston()
  sim5 = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
  bool5 = abs(sim5.𝒪est[:l∞]-2) < testTol

  tab = constructBogakiShampine()
  sim6 = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
  bool6 = abs(sim6.𝒪est[:l∞]-3) < testTol

  println("Higher Order")
  Δts = 1.//2.^(7:-1:4)
  tab = constructRKF()
  sim7 = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
  bool7 = abs(sim7.𝒪est[:l∞]-5) < testTol

  tab = constructDormandPrince()
  sim8 = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
  bool8 = abs(sim8.𝒪est[:l∞]-5) < testTol

  tab = constructCashKarp()
  sim9 = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
  bool9 = abs(sim9.𝒪est[:l∞]-5) < testTol

  println("Super High Order")
  #Need to make larger or else below machine ϵ
  Δts = 1./2.0.^(2:-1:-2)
  tab = constructRKF8()
  sim10 = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
  bool10 = abs(sim10.𝒪est[:l∞]-8) < testTol

  tab = constructDormandPrince8()
  sim11 = test_convergence(Δts,prob,alg=:ExplicitRK,tableau=tab)
  bool11 = abs(sim11.𝒪est[:l∞]-9.44) < testTol #Converges on linear at above 9!

  superbool1 = bool1 && bool2 && bool3 && bool4 && bool5 && bool6 && bool7 && bool8 && bool9 && bool10 && bool11

  println("Tests pass: $superbool1")
  ### Stiff Solvers

  println("Convergence Test on Stiff")
  Δts = 1.//2.^(8:-1:4)

  sim12 = test_convergence(Δts,prob,alg=:ImplicitEuler,autodiff=TEST_USE_FORWARDDIFF)
  bool12 = abs(sim12.𝒪est[:final]-1) < testTol

  sim13 = test_convergence(Δts,prob,alg=:Trapezoid,autodiff=false)
  bool13 = abs(sim13.𝒪est[:final]-2) < testTol

  sim14 = test_convergence(Δts,prob,alg=:Rosenbrock32)
  bool14 = abs(sim14.𝒪est[:final]-2) < testTol

  superbool2 = bool12 && bool13 && bool14

  println("Tests pass: $superbool2")
  superduperbool[i] = superbool1 && superbool2
end
minimum(superduperbool)
