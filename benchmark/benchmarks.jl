using DifferentialEquations, BenchmarkTools

const SUITE = BenchmarkGroup()

# =============================================================================
# ODE end-to-end solves
# =============================================================================

SUITE["ode"] = BenchmarkGroup()

function lotka!(du, u, p, t)
    du[1] = 1.5u[1] - u[1] * u[2]
    du[2] = -3.0u[2] + u[1] * u[2]
    return nothing
end
ode_prob = ODEProblem(lotka!, [1.0, 1.0], (0.0, 10.0))

function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
    return nothing
end
lorenz_prob = ODEProblem(lorenz!, [1.0, 0.0, 0.0], (0.0, 50.0))

SUITE["ode"]["tsit5"] = @benchmarkable solve($ode_prob, Tsit5())
SUITE["ode"]["default_alg"] = @benchmarkable solve($ode_prob)
SUITE["ode"]["vern7_lorenz"] = @benchmarkable solve($lorenz_prob, Vern7())

# =============================================================================
# Ensemble
# =============================================================================

SUITE["ensemble"] = BenchmarkGroup()

prob_func(prob, ctx) = remake(prob; u0 = prob.u0 .* 1.01)
eprob = EnsembleProblem(ode_prob; prob_func)
SUITE["ensemble"]["serial_100"] = @benchmarkable solve(
    $eprob, Tsit5(), EnsembleSerial(); trajectories = 100
)
