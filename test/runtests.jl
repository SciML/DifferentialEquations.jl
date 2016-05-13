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
@time include("femHeatDtTests.jl") #30 seconds at Δx=2^-5
println("Finite Element Heat Dx Tests")
@time include("femHeatDXtests.jl") #23 seconds
println("Finite Element Heat Method Tests")
@time include("femHeatMethodTest.jl") ##44 seconds
println("Finite Element Nonlinear Heat Methods Tests")
@time include("femHeatNonlinearMethodsTest.jl") #Instant
println("Finite Element Poisson Convergence Test")
@time include("femPoissonConvTest.jl") #Instant
println("Finite Element Nonlinear Poisson Tests")
@time include("femPoissonNonlinearTest.jl") #1.5 minutes
println("Heat Animation Test")
println("Finite Element Stochastic Poisson")
@time include("femStochasticPoissonSolo.jl") #Instant
println("Finite Element Poisson")
@time include("introductionExample.jl") #Instant
@time include("femHeatAnimationTest.jl") #2 minutes
println("Stochastic Heat Animation Test")
@time include("femStochasticHeatAnimationTest.jl") #1.5 minutes
