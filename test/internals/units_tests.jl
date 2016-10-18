using Unitful#, UnitfulPlots

using DifferentialEquations

f = (t,y) -> 0.5*y / 3.0u"s"
u₀ = 1.0u"N"
prob = ODEProblem(f,u₀)

sol =solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=:Midpoint)

sol =solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=:ExplicitRK,adaptive=true)

for alg in DifferentialEquations.DIFFERENTIALEQUATIONSJL_ALGORITHMS
  if !contains(string(alg),"Vectorized") && !contains(string(alg),"Threaded") && alg ∉ DifferentialEquations.DIFFERENTIALEQUATIONSJL_IMPLICITALGS
    sol = solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=alg,adaptive=true)
  end
end

println("Units for Number pass")

TEST_PLOT && plot(sol)

u₀ = [1.0u"N" 2.0u"N"
      3.0u"N" 1.0u"N"]

prob = ODEProblem(f,u₀)

sol =solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=:RK4)

sol =solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=:ExplicitRK)
sol =solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=:DP5)

sol =solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=:DP5Threaded)

for alg in DifferentialEquations.DIFFERENTIALEQUATIONSJL_ALGORITHMS
  println(alg)
  if alg ∉ DifferentialEquations.DIFFERENTIALEQUATIONSJL_IMPLICITALGS
    sol = solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=alg,adaptive=true)
  end
end

println("Units for 2D pass")

# Stochastic needs ΔW = s^(1/2).
#=
β = 0.6
σ = (t,y) -> β*y/(4.0s)
u₀ = 1.5Newton
prob = SDEProblem(f,σ,u₀)

sol =solve(prob::SDEProblem,[0,1],Δt=(1/2^4)Second,save_timeseries=true,alg=:EM)
sol =solve(prob::SDEProblem,[0,1],Δt=(1/2^4)Second,save_timeseries=true,alg=:SRIW1Optimized)

TEST_PLOT && plot(sol)

u₀ = [1.5Newton 2.0Newton
      3.0Newton 1.0Newton]

prob = SDEProblem(f,σ,u₀)

sol =solve(prob::SDEProblem,[0,1],Δt=(1/2^4)Second,save_timeseries=true,alg=:EM)
sol =solve(prob::SDEProblem,[0,1],Δt=(1/2^4)Second,save_timeseries=true,alg=:SRIW1Optimized)
=#

### Setup
Δx = (1//2^(3))u"m"
fem_mesh = notime_squaremesh([0u"m" 1u"m" 0u"m" 1u"m"],Δx,:dirichlet) #Fails at meshgrid because unit ranges not supported

f = (x) -> sin(2π.*map((y)->y.val,x[:,1])).*cos(2π.*map((y)->y.val,x[:,2]))
prob = PoissonProblem(f,D=1u"m^2/s")

#=
sol = solve(fem_mesh,prob)
=#

#Define a parabolic problem
T = 1u"s"
Δt = (1//2^(7))u"s"
fem_mesh = parabolic_squaremesh([0u"m" 1u"m" 0u"m" 1u"m"],Δx,Δt,T,:dirichlet)
f = (u,x,t)  -> (float(ones(size(x,1))))u"N" - .5u
u₀ = (x) -> map((x)->(x)Newton,zeros(size(x,1)))
prob = HeatProblem(u₀,f)
#=
println("Euler")
sol = solve(fem_mesh::FEMmesh,prob::HeatProblem,alg=:Euler)
=#

true
