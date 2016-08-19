# This definitely needs cleaning
using DifferentialEquations
probArr = Vector{DEProblem}(2)
bigprobArr = Vector{DEProblem}(2)
probArr[1] = linearODEExample()
probArr[2] = twoDimlinearODEExample()
bigprobArr[1] = linearODEExample(u₀=BigFloat(1),α=BigFloat(1))
bigprobArr[2] = twoDimlinearODEExample(α=ones(BigFloat,4,2),u₀=map(BigFloat,rand(4,2)).*ones(4,2)/2)
srand(100)
## Convergence Testing
println("Convergence Test on Linear")
Δts = 1.//2.^(8:-1:4)
testTol = 0.2
superduperbool = Vector{Bool}(2)

for i = 1:2
  prob = probArr[i]
  bigprob = bigprobArr[i]
  println("Special RKs")
  sim = test_convergence(Δts,prob,alg=:Euler)
  bool1 = abs(sim.𝒪est[:final]-1) < testTol
  sim2 = test_convergence(Δts,prob,alg=:Midpoint)
  bool2 = abs(sim2.𝒪est[:l∞]-2) < testTol
  sim3 = test_convergence(Δts,prob,alg=:RK4)
  bool3 = abs(sim3.𝒪est[:l∞]-4) < testTol

  superbool1 = bool1 && bool2 && bool3

  println("Tests pass: $superbool1")
  ### Stiff Solvers

  println("Convergence Test on Stiff")
  Δts = 1.//2.^(8:-1:4)

  sim12 = test_convergence(Δts,prob,alg=:ImplicitEuler,autodiff=true)
  sim122 = test_convergence(Δts,prob,alg=:ImplicitEuler,autodiff=false)
  bool12 = (abs(sim12.𝒪est[:final]-1) < testTol) && (abs(sim122.𝒪est[:final]-1) < testTol)
  sim13 = test_convergence(Δts,prob,alg=:Trapezoid,autodiff=true)
  sim132 = test_convergence(Δts,prob,alg=:Trapezoid,autodiff=false)
  bool13 = (abs(sim13.𝒪est[:final]-2) < testTol) && (abs(sim132.𝒪est[:final]-2) < testTol)
  sim14 = test_convergence(Δts,prob,alg=:Rosenbrock32)
  bool14 = abs(sim14.𝒪est[:final]-2) < testTol

  superbool2 = bool12 && bool13 && bool14

  println("Tests pass: $superbool2")
  superduperbool[i] = superbool1 && superbool2
end

minimum(superduperbool)
