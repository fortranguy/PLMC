if size(ARGS, 1) != 1
    error("Please provide a generating.json file.")
end
import PLMC
import JSON
inputData = JSON.parsefile(ARGS[1]; dicttype=Dict, use_mmap=true)

for iBox = 1:inputData["Environment"]["Boxes"]["number"]
    boxSize, components, interMinDistances = PLMC.set(iBox, inputData)

    numParticles = 0
    for iComponent = 1:size(components, 1)
        numParticles += components[iComponent].num
    end

    minDistanceDelta = 1e-6
    minDistanceMax = maximum(interMinDistances) + minDistanceDelta
    edge = sqrt(2) * minDistanceMax
    numbers = floor(Int64, 2 * boxSize / edge)
    for i = 1:size(numbers, 1)
        if isodd(numbers[i])
            numbers[i] -= (boxSize[i] - (numbers[i] - 1) * edge / 2) < minDistanceMax ? 1 : 0
        end
    end
    numExcluded = 0
    for k = 0:numbers[3]-1, j = 0:numbers[2]-1, i = 0:numbers[1]-1
        numExcluded += isodd(i + j + k) ? 1 : 0
    end
    numParticlesMax = prod(numbers) - numExcluded
    if numParticles > numParticlesMax
        error("Too many particles. Max = ", numParticlesMax)
    else
        println("Capacity = ", numParticles / numParticlesMax * 100, "%")
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
    ijkOffset = map(i -> iseven(i) ? 1 : 0, numbers)
    for k = 0:numbers[3]-1, j = 0:numbers[2]-1, i = 0:numbers[1]-1
        if isodd(i + j + k)
            continue
        end
        iCounter += 1
        if in(iCounter, toDelete)
            continue
        end
        iParticle += 1
        positions[:, iParticle] = PLMC.folded(boxSize,
            ([i, j, k] + ijkOffset) * edge/2 - boxSize/2)
    end
    positions = positions[:, randperm(size(positions, 2))]

    iParticleMin = 1
    for iComponent = 1:size(components, 1)
        iParticleMax = iParticleMin - 1 + components[iComponent].num
        components[iComponent].positions = positions[:, iParticleMin:iParticleMax]
        iParticleMin = iParticleMax + 1
    end
    PLMC.write(iBox, components, boxSize, inputData)
end
