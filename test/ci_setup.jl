Pkg.clone("https://github.com/JuliaODE/ODE.jl")
Pkg.checkout("ODE")

Pkg.clone("https://github.com/JuliaDiffEq/Sundials.jl")
Pkg.checkout("Sundials")
Pkg.build("Sundials")
