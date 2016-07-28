using DifferentialEquations, Plots
srand(100)

prob = linearODEExample()
sol3 =solve(prob::ODEProblem,[0,1],Δt=1/2^(6),save_timeseries=true,alg=:RK4,abstol=1,reltol=0)

prob = linearODEExample(u₀=BigInt(1)//BigInt(2))
sol =solve(prob::ODEProblem,[0,1],Δt=BigInt(1)//BigInt(2)^(6),save_timeseries=true,alg=:RK4,abstol=1,reltol=0)
sol2 =solve(prob::ODEProblem,[0,1],Δt=BigInt(1)/BigInt(2)^(6),save_timeseries=true,alg=:RK4,abstol=1,reltol=0)

sol.u
sol2.u
sol3.u
bool1 = 1.16e-16 < sol.u - sol3.u < 1.17e-16
bool2 = 1.16e-16 < sol2.u - sol3.u < 1.17e-16
bool3 = sol2.u - sol.u < big(1.73e-77)
bool4 = typeof(sol.u) == Rational{BigInt}
bool5 = typeof(sol2.u) == BigFloat
bool6 = typeof(sol3.u) == Float64

sol4 =solve(prob::ODEProblem,[0,1],Δt=BigInt(1)//BigInt(2)^(6),save_timeseries=true,alg=:ExplicitRK,abstol=1,reltol=0)

tab = constructDormandPrince8()
sol5 =solve(prob::ODEProblem,[0,1],Δt=BigInt(1)//BigInt(2)^(6),save_timeseries=true,alg=:ExplicitRK,abstol=1,reltol=0,tableau=tab)

bool7 = 1e-10 < float(sol5.u) - sol3.u < 1e-9
bool8 = typeof(sol5.u) == BigFloat

sol6 =solve(prob::ODEProblem,[0,1],Δt=BigInt(1)//BigInt(2)^(6),save_timeseries=true,alg=:ExplicitRK,adaptive=true,abstol=1,reltol=0,tableau=tab)

bool9 = 1e-10 < float(sol6.u) - sol3.u < 1e-9
bool10 = typeof(sol6.u) == BigFloat

bool1 && bool2 && bool3 && bool4 && bool5 && bool6 && bool7 && bool8 && bool9 && bool10
