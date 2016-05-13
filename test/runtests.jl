#!/usr/bin/env julia
using DifferentialEquations
using Base.Test

# Run tests
include("femHeatAnimationTest.jl") #2 minutes
include("femHeatDtTests.jl") #30 seconds at Δx=2^-5
include("femHeatDXtests.jl") #23 seconds
include("femHeatMethodTest.jl") ##44 seconds
include("femHeatNonlinearMethodsTest.jl") #Instant
include("femPoissonConvTest.jl") #Instant
include("femPoissonNonlinearTest.jl") #1.5 minutes
include("femStochasticHeatAnimationTest.jl") #1.5 minutes
include("femStochasticPoissonSolo.jl") #Instant
include("introductionExample.jl") #Instant
