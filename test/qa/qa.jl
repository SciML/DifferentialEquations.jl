using SciMLTesting, DifferentialEquations, JET, Test

run_qa(
    DifferentialEquations;
    explicit_imports = true,
    jet_kwargs = (; target_defined_modules = true),
)
