using Documenter, DifferentialEquations

makedocs(modules=[DifferentialEquations],doctest=false)

deploydocs(deps   = Deps.pip("mkdocs", "mkdocs-material", "python-markdown-math"),
    repo = "github.com/ChrisRackauckas/DifferentialEquations.jl.git",
    julia  = "0.5.0-",
    osname = "linux")
