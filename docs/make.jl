using Documenter, DifferentialEquations

makedocs(modules=[DifferentialEquations],
         doctest=false, clean=true,
         format =Documenter.Formats.HTML,
         sitename="DifferentialEquations.jl",
         authors="Chris Rackauckas",
         pages = Any[
         "Home" => "index.md",
         "Tutorials" => Any[
             "tutorials/ode_example.md",
             "tutorials/sde_example.md",
             "tutorials/fempoisson_example.md",
             "tutorials/femheat_example.md",
             "tutorials/femstochastic_example.md"
         ],
         "Solver Options" => Any[
           "solvers/ode_solve.md",
           "solvers/sde_solve.md",
           "solvers/dae_solve.md",
           "solvers/fempoisson_solve.md",
           "solvers/femheat_solve.md",
           "solvers/fdmstokes_solve.md"
         ],
         "Manual" => Any[
           "man/overview.md",
           "man/ODEProblem.md",
           "man/SDEProblem.md",
           "man/FEMProblem.md",
           "man/StokesProblem.md",
           "man/mesh.md",
           "man/solution.md",
           "man/output_specification.md",
           "man/callback_functions.md",
           "man/plot.md",
           "man/function_definition_macros.md",
           "man/benchmarks.md",
           "man/convergence.md",
           "man/conditional_dependencies.md",
           "man/progress_bar.md"
         ],
         "Internal Documentation" => Any[
           "internals/contributors_guide.md",
           "internals/fem_tools.md",
           "internals/extras.md",
           "internals/solver_helpers.md",
           "internals/notes_on_algorithms.md"
         ]
         ])


deploydocs(
   repo = "github.com/JuliaDiffEq/DifferentialEquations.jl.git",
   target = "build",
   deps = nothing,
   make = nothing)
