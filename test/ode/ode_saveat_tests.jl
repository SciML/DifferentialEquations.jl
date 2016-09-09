using DifferentialEquations, Plots

bools = Vector{Bool}(0)
prob = prob_ode_linear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=false,alg=:DP5,dense=false)
sol2=solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=false,alg=:DP5,dense=false,saveat=[1/2])

push!(bools,symdiff(sol.t,sol2.t) == [1/2])

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:DP5,dense=true)
sol2=solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:DP5,dense=true,saveat=[1/2])

sol2(.49)

interpd = sol2(0:1//2^(4):1)

plot(0:1//2^(4):1,interpd)

push!(bools,symdiff(sol.t,sol2.t) == [1/2])

sol =solve(prob::ODEProblem,[0,1];Δt=1/2^(2),save_timeseries=true,alg=:RK4,dense=false)
sol2=solve(prob::ODEProblem,[0,1];Δt=1/2^(2),save_timeseries=true,alg=:RK4,dense=false,saveat=[.125,.6,.61,.8])

push!(bools,symdiff(sol.t,sol2.t) == [.125,.6,.61,.8])

sol =solve(prob::ODEProblem,[0,1];Δt=1/2^(2),save_timeseries=true,alg=:Rosenbrock32,dense=false)
sol2=solve(prob::ODEProblem,[0,1];Δt=1/2^(2),save_timeseries=true,alg=:Rosenbrock32,dense=false,saveat=[.125,.6,.61,.8])

push!(bools,symdiff(sol.t,sol2.t) == [.125,.6,.61,.8])

sol =solve(prob::ODEProblem,[0,1];Δt=1/2^(2),save_timeseries=true,alg=:Trapezoid,dense=false)
sol2=solve(prob::ODEProblem,[0,1];Δt=1/2^(2),save_timeseries=true,alg=:Trapezoid,dense=false,saveat=[.125,.6,.61,.8])

push!(bools,symdiff(sol.t,sol2.t) == [.125,.6,.61,.8])

prob = prob_ode_2Dlinear

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=false,alg=:DP5,dense=true)
sol2=solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=false,alg=:DP5,dense=true,saveat=[1/2])

push!(bools,symdiff(sol.t,sol2.t) == [1/2])

sol =solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:DP5,dense=false)
sol2=solve(prob::ODEProblem,[0,1];Δt=1//2^(2),save_timeseries=true,alg=:DP5,dense=false,saveat=[1/2])

push!(bools,symdiff(sol.t,sol2.t) == [1/2])

sol =solve(prob::ODEProblem,[0,1];Δt=1/2^(2),save_timeseries=true,alg=:RK4,dense=false)
sol2=solve(prob::ODEProblem,[0,1];Δt=1/2^(2),save_timeseries=true,alg=:RK4,dense=false,saveat=[.125,.6,.61,.8])

push!(bools,symdiff(sol.t,sol2.t) == [.125,.6,.61,.8])

sol =solve(prob::ODEProblem,[0,1];Δt=1/2^(2),save_timeseries=true,alg=:Rosenbrock32,dense=false)
sol2=solve(prob::ODEProblem,[0,1];Δt=1/2^(2),save_timeseries=true,alg=:Rosenbrock32,dense=false,saveat=[.125,.6,.61,.8])

push!(bools,symdiff(sol.t,sol2.t) == [.125,.6,.61,.8])

sol =solve(prob::ODEProblem,[0,1];Δt=1/2^(2),save_timeseries=true,alg=:Trapezoid,dense=false)
sol2=solve(prob::ODEProblem,[0,1];Δt=1/2^(2),save_timeseries=true,alg=:Trapezoid,dense=false,saveat=[.125,.6,.61,.8])

push!(bools,symdiff(sol.t,sol2.t) == [.125,.6,.61,.8])

minimum(bools)
