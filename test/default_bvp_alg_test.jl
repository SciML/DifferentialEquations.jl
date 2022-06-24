using DifferentialEquations, Test

function f(du, u, p, t)
    (x, v) = u
    du[1] = v
    du[2] = -x
end

function bc!(resid, sol, p, t)
    resid[1] = sol[1][1]
    resid[2] = sol[end][1] - 1
end

tspan = (0.0, 100.0)
u0 = [0.0, 1.0]
bvp = BVProblem(f, bc!, u0, tspan)
resid_f = Array{Float64}(undef, 2)
sol = solve(bvp, Shooting(Tsit5()))
sol2 = solve(bvp)

@test sol2.alg == Tsit5()
@test all(sol.u .== sol2.u)
