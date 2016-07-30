######
##FEM Poisson Δx Convergence Tests
######
using DifferentialEquations#,LaTeXStrings

Δxs = 1.//2.^(4:-1:2) # 4 for testing, use 7 for good graph
prob = poissonProblemExample_wave()

sim = test_convergence(Δxs::AbstractArray,prob::PoissonProblem)

#Plot Result
TEST_PLOT && plot(sim,xguide="Delta x")

#Returns true if convergence is like Δx^2 in L2
sim.𝒪est[:L2]-2 <.1
