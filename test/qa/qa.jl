using DifferentialEquations
using Aqua: Aqua
using JET: JET
using Test

@testset "Aqua.jl" begin
    Aqua.test_all(DifferentialEquations)
end

@testset "JET.jl" begin
    JET.test_package(DifferentialEquations, target_defined_modules = true)
end
