"""
getNoise(N,node,elem;noiseType="White")

Returns a random vector corresponding to the noise type which was chosen.
"""
function getNoise(N,node,elem;noiseType="White")
  if noiseType=="White"
    return(sqrt(randn(N)))
  end
end
