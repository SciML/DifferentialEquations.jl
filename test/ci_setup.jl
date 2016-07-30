Pkg.clone("https://github.com/JuliaODE/ODE.jl")
Pkg.checkout("ODE","dev")
Pkg.checkout("GrowableArrays")
VERSION < v"0.5+" && Pkg.checkout("ChunkedArrays","v0.4-compat")
