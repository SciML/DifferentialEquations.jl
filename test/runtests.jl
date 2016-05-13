#!/usr/bin/env julia
#Automates Matplotlib installation for travis.ci
ENV["PYTHON"]=""
Pkg.build("PyCall")
using PyPlot

#Start Test Script
using DifferentialEquations
using Base.Test

testState = true
# Run tests

println("Finite Element Heat Dt Tests")
@time include("femHeatDtTests.jl")
println("Finite Element Heat Dx Tests")
@time include("femHeatDXtests.jl")
println("Finite Element Heat Method Tests")
@time include("femHeatMethodTest.jl")
println("Finite Element Nonlinear Heat Methods Tests")
@time include("femHeatNonlinearMethodsTest.jl")
println("Finite Element Poisson Convergence Test")
@time include("femPoissonConvTest.jl")
println("Finite Element Nonlinear Poisson Tests")
@time include("femPoissonNonlinearTest.jl")
println("Heat Animation Test")
println("Finite Element Stochastic Poisson")
@time include("femStochasticPoissonSolo.jl")
println("Finite Element Poisson")
@time include("introductionExample.jl")
@time include("femHeatAnimationTest.jl")
println("Stochastic Heat Animation Test")
@time include("femStochasticHeatAnimationTest.jl")
