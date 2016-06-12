## Disconnect from the Python Environment
## And use what's in the Conda package

Pkg.add("WinRPM")
using WinRPM
WinRPM.update()

ENV["PYTHON"]=""
Pkg.add("PyCall")
Pkg.add("PyPLot")
using PyPlot
