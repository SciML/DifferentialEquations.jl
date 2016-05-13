#!/usr/bin/env julia
using DiffEq
using Base.Test

# Run tests
include("femHeatAnimationTest.jl")
include("femHeatDtTests.jl")
include("femHeatDXtests.jl")
include("femHeatMethodTest.jl")
include("femHeatNonlinearMethodsTest.jl")
include("femPoissonConvTest.jl")
include("femPoissonNonlinearTest.jl")
include("femStochasticHeatAnimationTest.jl")
include("femStochasticPoissonSolo.jl")
include("introductionExample.jl")
