## Disconnect from the Python Environment
## And use what's in the Conda package

Pkg.add("WinRPM")
using WinRPM
WinRPM.update()

ENV["PYTHON"]=""
Pkg.build("PyCall")
using PyPlot
