if size(ARGS, 1) != 1
    error("Please provide a generating.json file.")
end
import PLMC
import JSON
generatingData = JSON.parsefile(ARGS[1]; dicttype=Dict, use_mmap=true)
import ProgressMeter
PM = ProgressMeter

for i_box = 1:generatingData["Environment"]["Boxes"]["number"]
    periodicBox = PLMC.newBox(i_box, generatingData)
    components = PLMC.newComponents(i_box, generatingData)
    minDistances = PLMC.newMinDistances(components, generatingData)

    testPosition = zeros(3, 0)
    for i_component = 1:size(components, 1)
        prog = PM.Progress(components[i_component].num, "Box $(i_box): Component $(i_component): ")
        while size(components[i_component].positions, 2) < components[i_component].num
            overlap = true
            while overlap
                testPosition = rand(3) .* periodicBox.edgesSize
                overlap = false
                for j_component = i_component:-1:1
                    for iParticle = 1:size(components[j_component].positions, 2)
                        if (PLMC.distance(periodicBox, testPosition,
                            components[j_component].positions[:, iParticle]) <
                            minDistances[i_component, j_component])
                            overlap = true
                            break
                        end
                    end
                    if (overlap)
                        break
                    end
                end
            end
            components[i_component].positions = hcat(components[i_component].positions, PLMC.
                folded(periodicBox, testPosition))
            PM.next!(prog)
        end
    end
    PLMC.writeCoordinates(i_box, components, periodicBox.edgesSize, generatingData)
end
