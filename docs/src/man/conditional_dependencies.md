# Conditional Dependencies

DifferentialEquations.jl is conditionally dependent on some packages which may
not be required for all users. The upside is that you can run DifferentialEquations.jl
without installing these packages. However, the downside is that you will have
to do the installation yourself (normal dependencies do a silent install). Luckily
DifferentialEquations.jl warns you about missing dependencies when calling a method
which requires one. This part of the manual will detail how to see if you're
missing a conditional dependency and how to alleviate the issue.

## The Conditional Dependency Notification

When a conditional dependency is required, DifferentialEquations.jl will import
the package the first time the method is called. When this happens, you will
receive a notification as follows:

```julia
[DifferentialEquations.jl] Initializing backend: PkgName
```

where PkgName is the name of the package it is importing. DifferentialEquations.jl
will then run a standard startup procedure on this package. If it fails, you
will receive the message

```julia
WARNING: Couldn't initialize PkgName.  (might need to install it?)
```

Most likely the issue is that you need to install the package. Go to the package's
Github repository for information on installing the package, and then come
back to try again. If that does not work, you may need the latest version of the
package by checking out master:

```julia
Pkg.checkout("PkgName")
```

If all else fails, please ask for help on [via the repository Gitter](https://gitter.im/ChrisRackauckas/DifferentialEquations.jl).

## What Methods Require Conditional Dependencies?

That's a good question! The implicit algorithms implemented in DifferentialEquations.jl
require [NLsolve.jl](https://github.com/EconForge/NLsolve.jl). Also, the `load`
function for the premade meshes requires [JLD.jl](https://github.com/JuliaIO/JLD.jl).

Lastly, there is a special conditional dependency for [Juno](http://junolab.org/). If
you are using Juno, then the progress bar functionality is works. If you're not
using Juno, then it won't do anything.

The other conditional dependencies are external solvers wrapped by DifferentialEquations.jl
Currently these include:

- ODE.jl
- ODEInterface.jl
- Sundials.jl

## Installation Instructions

For most of the conditional dependencies, the installation instructions are
standard. However, for some of the newest features, special instructions may
be required. The best way to stay up-to-date on this information is to checkout
the following resources:

- The packages which are conditional dependencies and use a standard installation
  can be found in the [/test/REQUIRE](https://github.com/ChrisRackauckas/DifferentialEquations.jl/blob/master/test/REQUIRE) file.
- Any special installation instructions are handled via [the ci_setup.jl file](https://github.com/ChrisRackauckas/DifferentialEquations.jl/blob/master/test/ci_setup.jl).

The current special installation instructions are as follows:

### ODE.jl.
The wrapper currently only works on the development branch of ODE.jl
at JuliaODE/ODE.jl. To install this version of ODE.jl, use the following commands:

```julia
Pkg.clone("https://github.com/JuliaODE/ODE.jl")
```

### Sundials.jl
The wrapper works on my fork until a PR goes through. You can install this via:

```julia
Pkg.clone("https://github.com/ChrisRackauckas/Sundials")
Pkg.checkout("Sundials","handles")
```
