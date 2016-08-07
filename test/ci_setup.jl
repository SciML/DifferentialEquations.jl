using Compat
Pkg.clone("https://github.com/JuliaODE/ODE.jl")
Pkg.checkout("ODE","dev")

if @compat is_windows()
  ENV["ODEINTERFACE_GFORTRAN"]="C:\mingw-w64\i686-5.3.0-posix-dwarf-rt_v4-rev0\mingw32\bin"
end
