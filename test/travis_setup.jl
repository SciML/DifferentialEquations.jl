#Pkg.clone("https://github.com/ChrisRackauckas/VectorizedRoutines.jl")
#Pkg.clone("https://github.com/ChrisRackauckas/ResettableStacks.jl")
Pkg.checkout("GrowableArrays")
VERSION < v"0.5+" && Pkg.checkout("ChunkedArrays","v0.4-compat")
