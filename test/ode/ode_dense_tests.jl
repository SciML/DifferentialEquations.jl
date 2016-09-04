using DifferentialEquations, Plots

bools = Vector{Bool}(0)
prob = prob_ode_linear
sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(3),save_timeseries=true,alg=:Euler,dense=true)

interpd = sol(0:1//2^(4):1)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:Euler,dense=true)

sol3 =solve(prob::ODEProblem,[0,1];Δt=1//2^(5),save_timeseries=true,alg=:Euler,dense=true)

TEST_PLOT && plot(sol2)
TEST_PLOT && plot!(float(sol2.t),interpd)
TEST_PLOT && plot!(float(sol3.t[1:2:end]),sol3.timeseries[1:2:end])

prob = prob_ode_2Dlinear
sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:Euler,dense=true)

interpd = DifferentialEquations.hermite3_interpolate(collect(0:1//2^(4):1),sol.t,sol.timeseries,sol.k)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:Euler,dense=true)

push!(bools,maximum(map(maximum,sol2[:] - interpd)) < .2)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:Euler,dense=false)

push!(bools,sol(0.5) == nothing)

prob = prob_ode_linear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:RK4,dense=true)

interpd = DifferentialEquations.hermite3_interpolate(collect(0:1//2^(4):1),sol.t,sol.timeseries,sol.k)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:RK4,dense=true)

push!(bools,maximum(map(maximum,sol2[:] - interpd)) < .2)

TEST_PLOT && plot(sol2)
TEST_PLOT && plot!(float(sol2.t),interpd)

minimum(bools)
