using DifferentialEquations

f(x) = 2x
f(x,y) = 2xy

numpar = numparameters(f) # Should be 2

numpar == 2
