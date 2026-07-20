using SciMLTesting, DifferentialEquations, JET, Test

# Only dependency-owned reexports are ignored; any owned public name must render.
const DEPENDENCY_OWNED_REEXPORTS = Tuple(
    name for name in public_api_names(DifferentialEquations)
        if parentmodule(getfield(DifferentialEquations, name)) !== DifferentialEquations
)

run_qa(
    DifferentialEquations;
    api_docs_kwargs = (; rendered = true, rendered_ignore = DEPENDENCY_OWNED_REEXPORTS),
    explicit_imports = true,
    jet_kwargs = (; target_defined_modules = true),
)
