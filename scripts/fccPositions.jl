if size(ARGS, 1) != 1
    error("Please provide a generating.json file.")
end
import PLMC
import JSON
generatingData = JSON.parsefile(ARGS[1]; dicttype=Dict, use_mmap=true)

for i_box = 1:generatingData["Environment"]["Boxes"]["number"]
    periodicBox = PLMC.newBox(i_box, generatingData)
    boxSizeSave = periodicBox.edgesSize
    components = PLMC.newComponents(i_box, generatingData)
    minDistances = PLMC.newMinDistances(components, generatingData)
    numParticles = 0
    for i_component = 1:size(components, 1)
        numParticles += components[i_component].num
    end

    minDistanceDelta = 1e-6
    periodicBox.edgesSize -= minDistanceDelta
    minDistanceMax = maximum(minDistances) + minDistanceDelta
    edge = sqrt(2) * minDistanceMax
    numbers = floor(Int64, 2 * periodicBox.edgesSize / edge)
    for i = 1:size(numbers, 1)
        if isodd(numbers[i])
            numbers[i] -=
                (periodicBox.edgesSize[i] - (numbers[i] - 1) * edge / 2) < minDistanceMax ? 1 : 0
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
    i_delete = 0
    while i_delete < size(toDelete, 1)
        iRand = rand(1:numParticlesMax)
        if !in(iRand, toDelete)
            i_delete += 1
            toDelete[i_delete] = iRand
        end
    end

    positions = zeros(Float64, 3, numParticles)
    i_counter = i_particle = 0
    ijk_offset = map(i -> iseven(i) ? 1 : 0, numbers)
    for k = 0:numbers[3]-1, j = 0:numbers[2]-1, i = 0:numbers[1]-1
        if isodd(i + j + k)
            continue
        end
        i_counter += 1
        if in(i_counter, toDelete)
            continue
        end
        i_particle += 1
        positions[:, i_particle] = PLMC.folded(periodicBox,
            ([i, j, k] + ijk_offset) * edge/2 - periodicBox.edgesSize/2)
    end
    positions = positions[:, randperm(size(positions, 2))]

    i_particleMin = 1
    for i_component = 1:size(components, 1)
        i_particleMax = i_particleMin - 1 + components[i_component].num
        components[i_component].positions = positions[:, i_particleMin:i_particleMax]
        i_particleMin = i_particleMax + 1
    end
    PLMC.writeCoordinates(i_box, components, boxSizeSave, generatingData)
end
