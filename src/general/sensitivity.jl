type ODELocalSensitivity{T,T2} <: DESensitivity
  val::Dict{Symbol,Vector{T}}
  normseries::Dict{Symbol,Vector{T2}}
end

ODELocalSensitivity() = ODELocalSensitivity(Dict{Symbol,Vector{Number}}(),Dict{Symbol,Vector{Number}}())
function ODELocalSensitivity{T}(sensitivity_params::Vector{Symbol},sensitivity_series::Vector{Vector{T}})
  val = Dict{Symbol,Vector{T}}()
  for i in eachindex(sensitivity_params)
    val[sensitivity_params[i]] = sensitivity_series[i]
  end
  T2 = eltype(val[sensitivity_params[1]][1])
  normseries = Dict{Symbol,Vector{T2}}()
  for i in eachindex(sensitivity_params)
    cur_series = val[sensitivity_params[i]]
    tmp = Vector{T2}(length(sensitivity_series[1]))
    for j in eachindex(cur_series)
      tmp[j] = norm(cur_series[j])
    end
    normseries[sensitivity_params[i]] = tmp
  end
  ODELocalSensitivity{T,T2}(val,normseries)
end
