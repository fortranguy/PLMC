if size(ARGS, 1) != 1
    error("Please provide a generating.json file.")
end
import PLMC
import JSON
inputData = JSON.parsefile(ARGS[1]; dicttype=Dict, use_mmap=true)
import ProgressMeter
PM = ProgressMeter

boxSize, components, interMinDistances = PLMC.set(inputData)

testPosition = zeros(3, 0)
for iComponent = 1:size(components, 1)
    prog = PM.Progress(components[iComponent].num, "Component $(iComponent): ")
    while size(components[iComponent].positions, 2) < components[iComponent].num
        overlap = true
        while overlap
            testPosition = rand(3) .* boxSize
            overlap = false
            for jComponent = iComponent:-1:1
                for iParticle = 1:size(components[jComponent].positions, 2)
                    if (PLMC.distance(boxSize, testPosition, components[jComponent].
                        positions[:, iParticle]) < interMinDistances[iComponent, jComponent])
                        overlap = true
                        break
                    end
                end
                if (overlap)
                    break
                end
            end
        end
        components[iComponent].positions = hcat(components[iComponent].positions, PLMC.
            folded(boxSize, testPosition))
        PM.next!(prog)
    end
end
PLMC.write(components, boxSize, inputData)
