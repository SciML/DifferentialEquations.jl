using DifferentialEquations, Plots

bools = Vector{Bool}(0)
prob = prob_ode_linear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:Euler,dense=true)

interpd = sol(0:1//2^(4):1)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:Euler,dense=true)

sol3 =solve(prob::ODEProblem,[0,1];Δt=1//2^(5),save_timeseries=true,alg=:Euler,dense=true)

TEST_PLOT && plot(sol2)
TEST_PLOT && plot!(float(sol2.t),interpd)
TEST_PLOT && plot!(float(sol3.t[1:2:end]),sol3.timeseries[1:2:end])

prob = prob_ode_2Dlinear
sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:Euler,dense=true)

interpd = sol(0:1//2^(4):1)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:Euler,dense=true)

push!(bools,maximum(map((x)->maximum(abs(x)),sol2[:] - interpd)) < .2)

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:Euler,dense=false)

push!(bools,sol(0.5) == nothing)

prob = prob_ode_linear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:RK4,dense=true)

interpd = sol(0:1//2^(4):1)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),save_timeseries=true,alg=:RK4,dense=true)

push!(bools,maximum(map((x)->maximum(abs(x)),sol2[:] - interpd)) < .2)

TEST_PLOT && plot(sol2)
TEST_PLOT && plot!(float(sol2.t),interpd)

sol =solve(prob::ODEProblem,[0,1];save_timeseries=true,alg=:DP5,dense=true)

sol2 =solve(prob::ODEProblem,0:1//2^(4):1;save_timeseries=true,alg=:DP5,dense=true,adaptive=false)

interpd = sol(0:1//2^(4):1)
TEST_PLOT && plot(sol2.t,interpd)
TEST_PLOT && plot(sol)

push!(bools,maximum(map((x)->maximum(abs(x)),sol2[:] - interpd)) < .2)

prob = prob_ode_2Dlinear

sol =solve(prob::ODEProblem,[0,1];save_timeseries=true,alg=:DP5,dense=true)

sol2 =solve(prob::ODEProblem,0:1//2^(4):1;save_timeseries=true,alg=:DP5,dense=true,adaptive=false)

interpd = sol(0:1//2^(4):1)

push!(bools,maximum(map((x)->maximum(abs(x)),sol2[:] - interpd)) < .2)

prob = prob_ode_linear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:BS3,dense=true)

interpd = sol(0:1//2^(4):1)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),alg=:BS3,dense=true,adaptive=false)

push!(bools,maximum(map((x)->maximum(abs(x)),sol2[:] - interpd)) < .2)

prob = prob_ode_2Dlinear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:BS3,dense=true)

interpd = sol(0:1//2^(4):1)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),alg=:BS3,dense=true,adaptive=false)

push!(bools,maximum(map((x)->maximum(abs(x)),sol2[:] - interpd)) < .2)


prob = prob_ode_linear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:Tsit5,dense=true)

interpd = sol(0:1//2^(4):1)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),alg=:Tsit5,dense=true,adaptive=false)

push!(bools,maximum(map((x)->maximum(abs(x)),sol2[:] - interpd)) < .2)

prob = prob_ode_2Dlinear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:Tsit5,dense=true)

interpd = sol(0:1//2^(4):1)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),alg=:Tsit5,dense=true,adaptive=false)

push!(bools,maximum(map((x)->maximum(abs(x)),sol2[:] - interpd)) < .23)



prob = prob_ode_linear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:TanYam7,dense=true)

interpd = sol(0:1//2^(4):1)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),alg=:TanYam7,dense=true,adaptive=false)

push!(bools,maximum(map((x)->maximum(abs(x)),sol2[:] - interpd)) < .2)

prob = prob_ode_2Dlinear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:TanYam7,dense=true)

interpd = sol(0:1//2^(4):1)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),alg=:TanYam7,dense=true,adaptive=false)

push!(bools,maximum(map((x)->maximum(abs(x)),sol2[:] - interpd)) < .2)


prob = prob_ode_linear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:TsitPap8,dense=true)

interpd = sol(0:1//2^(4):1)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),alg=:TsitPap8,dense=true,adaptive=false)

push!(bools,maximum(map((x)->maximum(abs(x)),sol2[:] - interpd)) < .2)

prob = prob_ode_2Dlinear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:TsitPap8,dense=true)

interpd = sol(0:1//2^(4):1)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),alg=:TsitPap8,dense=true,adaptive=false)

push!(bools,maximum(map((x)->maximum(abs(x)),sol2[:] - interpd)) < .2)


prob = prob_ode_linear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:Feagin10,dense=true)

interpd = sol(0:1//2^(4):1)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),alg=:Feagin10,dense=true,adaptive=false)

push!(bools,maximum(map((x)->maximum(abs(x)),sol2[:] - interpd)) < .2)

prob = prob_ode_2Dlinear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:Feagin10,dense=true)

interpd = sol(0:1//2^(4):1)

sol2 =solve(prob::ODEProblem,[0,1];Δt=1//2^(4),alg=:Feagin10,dense=true,adaptive=false)

push!(bools,maximum(map((x)->maximum(abs(x)),sol2[:] - interpd)) < .2)



minimum(bools)
