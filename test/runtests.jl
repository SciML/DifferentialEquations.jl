#!/usr/bin/env julia
#Automates Matplotlib installation for travis.ci
ENV["PYTHON"]=""
Pkg.build("PyCall")
using PyPlot

#Start Test Script
using DifferentialEquations
using Base.Test

# Run tests
@time include("femHeatAnimationTest.jl") #2 minutes
@time include("femHeatDtTests.jl") #30 seconds at Δx=2^-5
@time include("femHeatDXtests.jl") #23 seconds
@time include("femHeatMethodTest.jl") ##44 seconds
@time include("femHeatNonlinearMethodsTest.jl") #Instant
@time include("femPoissonConvTest.jl") #Instant
@time include("femPoissonNonlinearTest.jl") #1.5 minutes
@time include("femStochasticHeatAnimationTest.jl") #1.5 minutes
@time include("femStochasticPoissonSolo.jl") #Instant
@time include("introductionExample.jl") #Instant
