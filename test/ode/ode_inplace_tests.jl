using DifferentialEquations
u₀=rand(300,20).*ones(300,20)/2
prob = twoDimlinearODEExample(α=6*ones(300,20);u₀=u₀)
prob2 = twoDimlinearODEExample!(;α=6*ones(300,20),u₀=u₀)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=false,alg=:Euler)
sol =solve(prob2::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=false,alg=:Euler)

@time sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(6),save_timeseries=false,alg=:Euler)
@time sol2 =solve(prob2::ODEProblem,[0,1];Δt=1//2^(6),save_timeseries=false,alg=:Euler)

prob = lorenzAttractorODEExample()
prob2 = lorenzAttractorODEExample!()

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=false,alg=:Euler)
sol =solve(prob2::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=false,alg=:Euler)

alloc1 = @allocated sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(6),save_timeseries=false,alg=:Euler)
alloc2 = @allocated sol2 =solve(prob2::ODEProblem,[0,1];Δt=1//2^(6),save_timeseries=false,alg=:Euler)

alloc2 < alloc1
