## Disconnect from the Python Environment
## And use what's in the Conda package

ENV["PYTHON"]=""
Pkg.build("PyCall")
using PyPlot
