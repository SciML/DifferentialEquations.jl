#!/usr/bin/env julia

using DifferentialEquations, Base.Test
println("Starting tests")
@time @testset "Default Discrete Algorithm" begin include("default_discrete_alg_test.jl") end
@time @testset "Default ODE Algorithm" begin include("default_ode_alg_test.jl") end
@time @testset "Default Steady State Algorithm" begin include("default_steady_state_alg_test.jl") end
@time @testset "Default SDE Algorithm" begin include("default_sde_alg_test.jl") end
@time @testset "Default RODE Algorithm" begin include("default_rode_alg_test.jl") end
@time @testset "Default DDE Algorithm" begin include("default_dde_alg_test.jl") end
@time @testset "Default DAE Algorithm" begin include("default_dae_alg_test.jl") end
@time @testset "Default FEM Algorithm" begin include("default_fem_alg_test.jl") end
