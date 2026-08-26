module DifferentialEquations

using Reexport: Reexport, @reexport
using PrecompileTools: @compile_workload, @setup_workload

@reexport using SciMLBase
@reexport using OrdinaryDiffEq

@setup_workload begin
    using SciMLBase: ODEProblem, solve
    using OrdinaryDiffEq: Rodas5P, Tsit5

    scalar_problem = ODEProblem((u, p, t) -> p * u, 0.5, (0.0, 1.0), 1.01)
    stiff_problem = ODEProblem(
        (du, u, p, t) -> begin
            y₁, y₂, y₃ = u
            k₁, k₂, k₃ = p
            du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
            du[2] = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
            du[3] = k₂ * y₂^2
        end,
        [1.0, 0.0, 0.0],
        (0.0, 1.0),
        (0.04, 3.0e7, 1.0e4)
    )

    @compile_workload begin
        solve(scalar_problem, Tsit5())
        solve(stiff_problem, Rodas5P())
    end
end

end # module
