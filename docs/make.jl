using Documenter, DifferentialEquations

makedocs(modules=[DifferentialEquations],
        doctest=true)

deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/ChrisRackauckas/DifferentialEquations.jl.git",
    julia  = "0.4.5",
    osname = "osx")
