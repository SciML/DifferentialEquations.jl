Pkg.clone("https://github.com/JuliaODE/ODE.jl")
Pkg.checkout("ODE")

Pkg.clone("https://github.com/alyst/Sundials.jl")
Pkg.checkout("Sundials","handles")
Pkg.build("Sundials")
