using DifferentialEquations, LaTeXStrings

Δx = 1//2^3 #Super large for quick testing. Lower this for better results.
mesh = FDMMesh(Δx,mins=[-1;-1],maxs=[1;1])
prob = homogeneousStokesExample()

#DGS Convergence
@time sol = solve(prob,mesh,converrors=true,maxiters=200,alg="DGS")

err1 = sol.converrors["rul∞"]
err2 = sol.converrors["rvl∞"]
err3 = sol.converrors["rpl∞"]
fig = PyPlot.figure("dgs",figsize=(10,10))
PyPlot.xlabel("Iterations")
PyPlot.ylabel("Error")
PyPlot.semilogy(1:length(err1),err1)
PyPlot.semilogy(1:length(err2),err2)
PyPlot.semilogy(1:length(err3),err3)
PyPlot.legend(["u","v","p"])
texstr = L"\Delta t = 2^{-3}"
PyPlot.title("DGS Convergence, $texstr")

#Multigrid Relative Convergence
@time sol = solve(prob,mesh,converrors=true,maxiters=20,alg="Multigrid",levels=4)

err1 = sol.converrors["rresul2"]
err2 = sol.converrors["rresvl2"]
err3 = sol.converrors["rrespl2"]
fig = PyPlot.figure("multigrid",figsize=(10,10))
PyPlot.xlabel("Iterations")
PyPlot.ylabel("Error")
PyPlot.semilogy(1:length(err1),err1)
PyPlot.semilogy(1:length(err2),err2)
PyPlot.semilogy(1:length(err3),err3)
PyPlot.legend(["u","v","p"])
texstr = L"\Delta t = 2^{-3}"
PyPlot.title("Multigrid Relative Residual Convergence, $texstr")

err1[end] < 1e-12
