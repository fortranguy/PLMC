import JSON
import ProgressMeter
PM = ProgressMeter
include("PLMC.jl")
import PLMC

if size(ARGS, 1) == 0
    error("Please provide a .json file.")
end
input_data = JSON.parsefile(ARGS[1]; dicttype=Dict, use_mmap=true)

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

for i_component = 1:size(components, 1)
    components[i_component].positions = zeros(3, 1) #excluded
end
first_positions = trues(size(components))
test_position = Array{Float64}
for i_component = 1:size(components, 1)
    prog = PM.Progress(components[i_component].num)
    while size(components[i_component].positions, 2) < components[i_component].num
        overlap = true
        while overlap
            test_position = rand(3) .* box_size
            overlap = false
            for j_component = i_component:-1:1
                for i_particle = 1:size(components[j_component].positions, 2)
                    if (PLMC.distance(test_position, components[j_component].positions[:,
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
        if (first_positions[i_component])
            components[i_component].positions[:, 1] = PLMC.folded(test_position, box_size)
            first_positions[i_component] = false
        else
            components[i_component].positions = hcat(components[i_component].positions,
                                                     PLMC.folded(test_position, box_size))
        end
        PM.next!(prog)
    end
    output_file = input_data["Mixture"]["Component $(i_component)"]["initial positions"]
    writedlm(output_file, components[i_component].positions')
    println("Positions written in ", output_file)
end
