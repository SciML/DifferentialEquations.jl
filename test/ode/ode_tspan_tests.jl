using DifferentialEquations, Plots
srand(100)

prob = linearODEExample()
sol3 =solve(prob::ODEProblem,[0,1//2,1],Δt=1//2^(6))

1//2 ∈ sol3.t
