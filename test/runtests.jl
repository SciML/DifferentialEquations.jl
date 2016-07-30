#!/usr/bin/env julia

const CPU_FLOPS = peakflops()
const TEST_USE_FORWARDDIFF = false
const TEST_PLOT = false
#Start Test Script
using DifferentialEquations, Plots
using Base.Test
gr()

# Run tests

tic()
#ODE
println("Linear ODE Tests")
@time @test include("ode/ode_twodimlinear_tests.jl")
println("ODE Convergence Tests")
@time @test include("ode/ode_convergence_tests.jl")
println("ODE Adaptive Tests")
@time @test include("ode/ode_adaptive_tests.jl")
println("ODE Lorenz Attractor")
@time @test include("ode/ode_lorenzattractor_tests.jl")
println("ODEInterface Tests")
@time @test include("ode/ODEInterface_tests.jl")
println("ODE.jl Tests")
@time @test include("ode/ODEJL_tests.jl")
println("ODE Number Type Tests")
@time @test include("ode/ode_numbertype_tests.jl")
println("ODE Initial Dt Tests")
@time @test include("ode/ode_initdt_tests.jl")

#SDE
println("Linear SDE Tests")
@time @test include("sde/sde_linear_tests.jl")
println("Two-dimensional Linear SDE Tests")
@time @test include("sde/sde_twodimlinear_tests.jl")
println("Rossler Order Tests")
@time @test include("sde/sde_rosslerorder_tests.jl")
println("SDE Convergence Tests")
@time @test include("sde/sde_convergence_tests.jl")
println("SDE Number Type Tests")
@time @test include("sde/sde_numbertype_tests.jl")

#Adaptive SDE
#=
println("Adaptive SDE Linear Tests")
@time @test include("sde/sde_linearadaptive_tests.jl")
println("Adaptive SDE Distribution Test")
@time @test include("sde/sde_adaptivedistribution_tests.jl")
println("Multiple Dimension Linear Adaptive Test")
@time @test include("sde/sde_twodimlinearadaptive_tests.jl")
println("SDE Autostepsize Test")
@time @test include("sde/sde_autostepsize_test.jl")
println("SDE Additive Lorenz Attractor Test")
@time @test include("sde/sde_lorenzattractor_tests.jl")
=#

#Finite Element
#Heat
println("Finite Element Heat Dt Tests")
@time @test include("heat/femheat_dtconvergence_tests.jl")
println("Finite Element Heat Dx Tests")
@time @test include("heat/femheat_dxconvergence_tests.jl")
println("Finite Element Heat Method Tests")
@time @test include("heat/femheat_methods_tests.jl")
println("Finite Element Nonlinear Heat Methods Tests")
@time @test include("heat/femheat_nonlinearmethods_tests.jl")
#println("Finite Element Nonlinear System Heat Tests") #ForwardDiff Issue
#@time @test include("heat/femheat_system_tests.jl")
println("Heat Animation Test")
@time @test include("heat/femheat_animation_tests.jl")
println("Stochastic Heat Animation Test")
@time @test include("heat/femheat_stochasticanimation_tests.jl")

#Poisson
println("Finite Element Poisson Convergence Test")
@time @test include("poisson/fempoisson_convergence_tests.jl")
println("Finite Element Nonlinear Poisson Tests")
@time @test include("poisson/fempoisson_nonlinear_tests.jl")
#println("Finite Element Nonlinear System Poisson Tests") #ForwardDiff Issue
#@time @test include("/poisson/fempoisson_nonlinearsystem_tests.jl")
println("Finite Element Stochastic Poisson")
@time @test include("poisson/fempoisson_stochastic_tests.jl")
println("Finite Element Poisson")
@time @test include("poisson/fempoisson_linear_tests.jl")

#FDM
##Stokes
println("Stokes Tests")
@time @test include("stokes/stokes_tests.jl")
println("DGS Internals Test")
@time @test include("stokes/stokes_dgs_tests.jl")

#Internals
println("Quadrature Points Tests")
@time @test include("internals/quadpts_test.jl")
println("Number of Parameters Calc Tests")
@time @test include("internals/numparameters_test.jl")
println("Assembly Tests")
@time @test include("internals/assembly_tests.jl")
println("Boundary Tests")
@time @test include("internals/boundary_tests.jl")
println("Example Mesh Tests")
@time @test include("internals/mesh_examples_tests.jl")
println("Simple Mesh Tests")
@time @test include("internals/mesh_SimpleMesh_tests.jl")

toc()
