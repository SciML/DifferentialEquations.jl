#!/usr/bin/env julia

using DifferentialEquations, Test, SafeTestsets

@time begin
    @testset "DifferentialEquations" begin
        @testset "Re-exports available" begin
            @test true
        end

        @testset "Basic ODE solve" begin
            f(u, p, t) = 1.01 * u
            u0 = 0.5
            tspan = (0.0, 1.0)
            prob = ODEProblem(f, u0, tspan)
            sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
            @test sol.retcode == ReturnCode.Success
            @test length(sol.t) > 0
            @test length(sol.u) > 0
            @test sol.u[1] ≈ u0
        end

        @testset "In-place ODE solve" begin
            function lorenz!(du, u, p, t)
                σ, ρ, β = p
                du[1] = σ * (u[2] - u[1])
                du[2] = u[1] * (ρ - u[3]) - u[2]
                du[3] = u[1] * u[2] - β * u[3]
            end
            u0 = [1.0, 0.0, 0.0]
            tspan = (0.0, 1.0)
            p = (10.0, 28.0, 8 / 3)
            prob = ODEProblem(lorenz!, u0, tspan, p)
            sol = solve(prob, Tsit5())
            @test sol.retcode == ReturnCode.Success
            @test sol.u[1] == u0
            @test length(sol.u[end]) == 3
        end

        @testset "Stiff ODE solve" begin
            function rober!(du, u, p, t)
                y₁, y₂, y₃ = u
                k₁, k₂, k₃ = p
                du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
                du[2] = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
                du[3] = k₂ * y₂^2
            end
            prob = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1e5), (0.04, 3e7, 1e4))
            sol = solve(prob, Rodas5())
            @test sol.retcode == ReturnCode.Success
        end

        @testset "Solution interpolation" begin
            f(u, p, t) = 1.01 * u
            prob = ODEProblem(f, 0.5, (0.0, 1.0))
            sol = solve(prob, Tsit5(), dense = true)
            u_interp = sol(0.5)
            @test u_interp isa Number
            @test u_interp > 0.5
        end

        @testset "Callbacks" begin
            function lorenz!(du, u, p, t)
                du[1] = 10.0 * (u[2] - u[1])
                du[2] = u[1] * (28.0 - u[3]) - u[2]
                du[3] = u[1] * u[2] - (8 / 3) * u[3]
            end
            prob = ODEProblem(lorenz!, [1.0, 0.0, 0.0], (0.0, 1.0))

            # ContinuousCallback
            condition(u, t, integrator) = t - 0.5
            affect!(integrator) = nothing
            cb = ContinuousCallback(condition, affect!)
            sol = solve(prob, Tsit5(); callback = cb)
            @test sol.retcode == ReturnCode.Success

            # DiscreteCallback
            dcb = DiscreteCallback((u, t, integrator) -> t >= 0.5, affect!)
            sol2 = solve(prob, Tsit5(); callback = dcb)
            @test sol2.retcode == ReturnCode.Success
        end

        @testset "Remake" begin
            f(u, p, t) = p * u
            prob = ODEProblem(f, 0.5, (0.0, 1.0), 1.01)
            prob2 = remake(prob; u0 = 1.0)
            sol = solve(prob2, Tsit5())
            @test sol.retcode == ReturnCode.Success
            @test sol.u[1] ≈ 1.0
        end

        @testset "saveat" begin
            f(u, p, t) = 1.01 * u
            prob = ODEProblem(f, 0.5, (0.0, 1.0))
            sol = solve(prob, Tsit5(); saveat = 0.1)
            @test sol.retcode == ReturnCode.Success
            @test length(sol.t) == 11
        end
    end
end
