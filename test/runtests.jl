#!/usr/bin/env julia

using DifferentialEquations, Base.Test
@time @testset "Default ODE Algorithm" begin include("default_ode_alg_test.jl") end
