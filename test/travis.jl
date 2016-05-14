#Automates Matplotlib installation for travis.ci

ENV["PYTHON"]=""
Pkg.build("PyCall")
using PyPlot
