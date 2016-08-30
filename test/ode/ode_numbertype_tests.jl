using DifferentialEquations, Plots
srand(100)
setprecision(400)

f = (t,u) -> (2u)
analytic = (t,u₀) -> u₀*exp(t)
"""Linear ODE on Float64"""
prob_ode_linear = ODEProblem(f,1/2,analytic=analytic)


prob = prob_ode_linear
sol3 =solve(prob::ODEProblem,[0,1],Δt=1/2^(6),save_timeseries=true,alg=:RK4,abstol=1,reltol=0)

prob = ODEProblem(f,BigInt(1)//BigInt(2),analytic=analytic)

sol =solve(prob::ODEProblem,[0,1],Δt=BigInt(1)//BigInt(2)^(6),save_timeseries=true,alg=:RK4,abstol=1,reltol=0)
sol2 =solve(prob::ODEProblem,[0,1],Δt=BigInt(1)/BigInt(2)^(6),save_timeseries=true,alg=:RK4,abstol=1,reltol=0)

sol.u
sol2.u
sol3.u
bool1 = 1.93e-16 < abs(sol.u - sol3.u) < 1.94e-16
bool2 = 1.93e-16  < abs(sol2.u - sol3.u) <1.94e-16
bool3 = sol2.u - sol.u < big(1.73e-77)
bool4 = typeof(sol.u) == Rational{BigInt}
bool5 = typeof(sol2.u) == Rational{BigInt}
bool6 = typeof(sol3.u) == Float64

sol4 =solve(prob::ODEProblem,[0,1],Δt=BigInt(1)//BigInt(2)^(3),save_timeseries=true,alg=:DP5,adaptive=false)

bool9 = typeof(sol4.u) == Rational{BigInt}

tab = constructDormandPrince8_64bit(Rational{BigInt})
sol5 =solve(prob::ODEProblem,[0,1],Δt=BigInt(1)//BigInt(2)^(3),save_timeseries=true,alg=:ExplicitRK,abstol=1,reltol=0,tableau=tab,adaptive=false)

bool7 = 5.72e-8 < abs(float(sol5.u) - sol3.u) < 5.73e-8
bool8 = typeof(sol5.u) == Rational{BigInt}

bool1 && bool2 && bool3 && bool4 && bool5 && bool6 && bool7 && bool8 && bool9
