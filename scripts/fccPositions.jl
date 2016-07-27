module PLMC
    function folded(boxSize::Array{Float64, 1}, x::Array{Float64, 1})
        v = mod(x, boxSize)
        for i=1:3
            if (v[i] > boxSize[i]/2)
                v[i] -= boxSize[i]
            end
        end
        return v
    end
end

import JSON
import PLMC

if size(ARGS, 1) == 0
    error("Please provide a .json file.")
end
inputData = JSON.parsefile(ARGS[1]; dicttype=Dict, use_mmap=true)

boxSize = map(Float64, inputData["Environment"]["Box"]["initial size"])
numComponents = inputData["Mixture"]["number of components"]
if numComponents == 0
    exit(0)
elseif numComponents > 1
    error("numComponents must be 1.")
end

numParticles = inputData["Mixture"]["Component 1"]["initial number"]
minDistanceDelta = 1e-6
minDistance = inputData["Mixture"]["Component 1"]["minimum distance"] + minDistanceDelta
edge = sqrt(2) * minDistance
numbers = floor(Int64, 2 * boxSize / edge)
numExcluded =0
for k = 0:numbers[3]-1, j = 0:numbers[2]-1, i = 0:numbers[1]-1
    numExcluded += isodd(i + j + k) ? 1 : 0
end
numParticlesMax = prod(numbers) - numExcluded
if numParticles > numParticlesMax
    error("Too many particles. Max =", numParticlesMax)
end

toDelete = zeros(Int64, numParticlesMax - numParticles)
iDelete = 0
while iDelete < size(toDelete, 1)
    iRand = rand(1:numParticlesMax)
    if !in(iRand, toDelete)
        iDelete += 1
        toDelete[iDelete] = iRand
    end
end

positions = zeros(Float64, 3, numParticles)
iCounter = iParticle = 0
for k = 0:numbers[3]-1, j = 0:numbers[2]-1, i = 0:numbers[1]-1
    if isodd(i + j + k)
        continue
    end
    iCounter += 1
    if in(iCounter, toDelete)
        continue
    end
    iParticle += 1
    positions[:, iParticle] = PLMC.folded(boxSize, [i, j, k] * edge/2)
end

output_file = open(inputData["Mixture"]["Component 1"]["initial coordinates"], "w")
println(output_file, "# number    ", numParticles)
println(output_file, "# position_x  position_y  position_z")
writedlm(output_file, positions')
