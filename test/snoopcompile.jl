using SnoopCompile

### Log the compiles
# This only needs to be run once (to generate "/tmp/images_compiles.csv")

SnoopCompile.@snoop "/tmp/diffeq_compiles.csv" begin
    include(Pkg.dir("DifferentialEquations", "test","runtests.jl"))
end

using DifferentialEquations

data = SnoopCompile.read("/tmp/diffeq_compiles.csv")

#=
# The Images tests are run inside a module ImagesTest, so all
# the precompiles get credited to ImagesTest. Credit them to Images instead.
subst = Dict("ImagesTests"=>"Images")
=#

#=
# Blacklist helps fix problems:
# - MIME uses type-parameters with symbols like :image/png, which is
#   not parseable
blacklist = ["MIME"]
=#

# Use these two lines if you want to create precompile functions for
# individual packages
pc, discards = SnoopCompile.parcel(data[end:-1:1,2])#, subst=subst, blacklist=blacklist)
SnoopCompile.write("/tmp/precompile", pc)
