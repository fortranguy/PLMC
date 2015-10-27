import JSON
json = JSON
import ProgressMeter
pm = ProgressMeter
import PLMC
pbc = PLMC.PBC

if size(ARGS, 1) == 0
    error("Please provide a .json file.")
end
input_data = json.parsefile(ARGS[1]; dicttype=Dict, use_mmap=true)

box_size = map(Float64, input_data["Environment"]["Box"]["size"])
components = Array{PLMC.Component}(input_data["Mixture"]["number of components"])
if (size(components, 1) == 0)
    exit(0)
end

interMinDistance = zeros(size(components, 1), size(components, 1)) # triangle?
for i_component = 1:size(components, 1)
    components[i_component] = PLMC.Component(
        input_data["Mixture"]["Component $(i_component)"]["number"], zeros(3, 1), zeros(3, 1))
    for j_component = 1:i_component-1
        interMinDistance[i_component, j_component] =
            input_data["Mixture"]["Inter $(j_component)$(i_component)"]["minimum distance"]
    end
    interMinDistance[i_component, i_component] =
        input_data["Mixture"]["Component $(i_component)"]["minimum distance"]
end
num_particles = 0
for i_component=1:size(components, 1)
    num_particles += components[i_component].num
end

components[1].positions[:, 1] = pbc.folded(rand(3) .* box_size, box_size)
test_position = Array{Float64}
for i_component = 1:size(components, 1)
    prog = pm.Progress(components[i_component].num)
    while size(components[i_component].positions, 2) < components[i_component].num
        overlap = true
        while overlap
            test_position = rand(3) .* box_size
            overlap = false
            for j_component = i_component:-1:1
                for i_particle = 1:size(components[j_component].positions, 2)
                    if (pbc.distance(test_position, components[j_component].positions[:,
                                                                                      i_particle],
                                     box_size) < interMinDistance[i_component, j_component])
                        overlap = true
                        break
                    end
                end
                if (overlap)
                    break
                end
            end
        end
        components[i_component].positions = hcat(components[i_component].positions,
                                                 pbc.folded(test_position, box_size))
        pm.next!(prog)
    end
    output_file = input_data["Mixture"]["Component $(i_component)"]["initial positions"]
    writedlm(output_file, components[i_component].positions')
    println()
    println("Positions written in ", output_file)
end
