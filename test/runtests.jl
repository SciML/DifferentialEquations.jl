#!/usr/bin/env julia

using DifferentialEquations, Test, SafeTestsets
@time begin
    @time @safetestset "Default Discrete Algorithm" begin
        include("default_discrete_alg_test.jl")
    end
    @time @safetestset "Default Steady State Algorithm" begin
        include("default_steady_state_alg_test.jl")
    end
    @time @safetestset "Default SDE Algorithm" begin
        include("default_sde_alg_test.jl")
    end
    @time @safetestset "Default RODE Algorithm" begin
        include("default_rode_alg_test.jl")
    end
    @time @safetestset "Default DDE Algorithm" begin
        include("default_dde_alg_test.jl")
    end
    @time @safetestset "Default DAE Algorithm" begin
        include("default_dae_alg_test.jl")
    end
    @time @safetestset "Default BVP Algorithm" begin
        include("default_bvp_alg_test.jl")
    end
end
