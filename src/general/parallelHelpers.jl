function shapeResult(res)
  #Input is a vector of tuples
  #Output the "columns" as vectors
  out = cell(length(res[1]))
  for i = 1:length(res[1])
    out[i] = Vector{typeof(res[1][i])}(length(res))
  end
  for i = 1:length(res), j = 1:length(out)
    out[j][i] = res[i][j]
  end
  return tuple(out...)
end

function monteCarlo(func::Function,args...)
  res = pmap(func,args...)
  if length(res[1])==1
    return res
  else
    return shapeResult(res)
  end
end
