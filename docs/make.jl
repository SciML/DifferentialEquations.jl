using Documenter, DifferentialEquations

makedocs(modules=[DifferentialEquations],doctest=false)

deploydocs(deps   = Deps.pip("mkdocs", "mkdocs-material", "python-markdown-math"),
    repo = "github.com/ChrisRackauckas/DifferentialEquations.jl.git")
