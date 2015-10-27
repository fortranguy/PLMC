module RandomOrientations
    function randomOrientations(num_particles::Int64)
        orientations = randn(num_particles, 3)
        for i=1:size(orientations, 1)
            orientations[i, :] = orientations[i, :] / norm(orientations[i, :])
        end
        return orientations
    end
end

import JSON
json = JSON
import RandomOrientations
ro = RandomOrientations
import PLMC

if size(ARGS, 1) == 0
    error("Please provide a .json file.")
end
input_data = json.parsefile(ARGS[1]; dicttype=Dict, use_mmap=true)
num_components = input_data["Mixture"]["number of components"]
if num_components == 0
    exit(0)
end
for i_component = 1:num_components
    if (input_data["Mixture"]["Component $(i_component)"]["is dipolar"])
        component_i = PLMC.Component(input_data["Mixture"]["Component $(i_component)"]["number"],
                                     zeros(3, 1), zeros(3, 1))
        component_i.orientations = ro.randomOrientations(component_i.num)
        output_file = input_data["Mixture"]["Component $(i_component)"]["initial orientations"]
        writedlm(output_file, component_i.orientations)
        println("Orientations written in ", output_file)
    end
end
