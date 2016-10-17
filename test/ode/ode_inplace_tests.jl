using DifferentialEquations
u₀=rand(300,20).*ones(300,20)/2
prob = prob_ode_2Dlinear_notinplace
prob2 = prob_ode_2Dlinear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=false,alg=:Euler)
sol =solve(prob2::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=false,alg=:Euler)

@time sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(6),save_timeseries=false,alg=:Euler)
@time sol2 =solve(prob2::ODEProblem,[0,1];Δt=1//2^(6),save_timeseries=false,alg=:Euler)

alloc1 = @allocated sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(6),save_timeseries=false,alg=:Euler)
alloc2 = @allocated sol2 =solve(prob2::ODEProblem,[0,1];Δt=1//2^(6),save_timeseries=false,alg=:Euler)

alloc1 = @allocated sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(6),save_timeseries=false,alg=:Euler)
alloc2 = @allocated sol2 =solve(prob2::ODEProblem,[0,1];Δt=1//2^(6),save_timeseries=false,alg=:Euler)

bool1 = alloc2 <= alloc1

sol = solve(prob_ode_large2Dlinear::ODEProblem,[0,1];Δt=1//2^(6),save_timeseries=true,alg=:Euler)
sol2 = solve(prob_ode_large2Dlinear::ODEProblem,[0,1],sol.timeseries,sol.t,sol.k;Δt=1//2^(8),save_timeseries=true,alg=:Euler)

sol = solve(prob_ode_large2Dlinear::ODEProblem,[0,1];Δt=1//2^(6),save_timeseries=true,alg=:Euler)
alloc1 = @allocated sol = solve(prob_ode_large2Dlinear::ODEProblem,[0,1];Δt=1//2^(8),save_timeseries=true,alg=:Euler)
alloc2 = @allocated sol2 = solve(prob_ode_large2Dlinear::ODEProblem,[0,1],sol.timeseries,sol.t,sol.k;Δt=1//2^(8),save_timeseries=true,alg=:Euler)

bool2 = alloc2 <= alloc1

bool1 && bool2
