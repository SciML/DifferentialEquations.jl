type NoiseProcess
  noise_type
  noise_func
end

const WHITE_NOISE = NoiseProcess(:Diagonal,randn)

"""
construct_correlated_noisefunc(Γ::AbstractArray)

Takes in a constant Covariance matrix Γ and spits out the noisefunc.
"""
function construct_correlated_noisefunc(Γ::AbstractArray)
  γ = svdfact(Γ)
  A = γ[:U]*Diagonal(√γ[:S])
  noise_func = (N...) -> A*randn(N...)
  NoiseProcess(:Commutative,noise_func)
end
