module PLMC

    type Component
        num :: Int64
        positions :: Array{Float64, 2}
        orientations :: Array{Float64, 2}
    end

    function folded(box_size::Array{Float64, 1}, x::Array{Float64, 1})
        v = mod(x, box_size)
        for i=1:3
            if (v[i] > box_size[i]/2)
                v[i] -= box_size[i]
            end
        end
        return v
    end

    function distance(box_size::Array{Float64, 1}, x::Array{Float64, 1}, y::Array{Float64, 1})
        return norm(folded(box_size, y - x))
    end

end

module RandomOrientations
    function randomOrientations(num_particles::Int64)
        orientations = randn(3, num_particles)
        for i=1:size(orientations, 2)
            orientations[:, i] = orientations[:, i] / norm(orientations[:, i])
        end
        return orientations
    end
end

import JSON
import RandomOrientations
RO = RandomOrientations
import PLMC
import ProgressMeter
PM = ProgressMeter

if size(ARGS, 1) == 0
    error("Please provide a .json file.")
end
input_data = JSON.parsefile(ARGS[1]; dicttype=Dict, use_mmap=true)

box_size = map(Float64, input_data["Environment"]["Box"]["initial size"])
components = Array{PLMC.Component}(input_data["Mixture"]["number of components"])
if (size(components, 1) == 0)
    exit(0)
end

interMinDistance = zeros(size(components, 1), size(components, 1)) # triangle?
for i_component = 1:size(components, 1)
    components[i_component] = PLMC.Component(
        input_data["Mixture"]["Component $(i_component)"]["initial number"], zeros(3, 1), zeros(3, 1))
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
                    if (PLMC.distance(box_size, test_position, components[j_component].
                        positions[:, i_particle]) < interMinDistance[i_component, j_component])
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
            components[i_component].positions[:, 1] = PLMC.folded(box_size, test_position)
            first_positions[i_component] = false
        else
            components[i_component].positions = hcat(components[i_component].positions,
                                                     PLMC.folded(box_size, test_position))
        end
        PM.next!(prog)
    end
    output_file = open(input_data["Mixture"]["Component $(i_component)"]["initial coordinates"],
        "w")
    println(output_file, "# number    ", components[i_component].num)
    if (input_data["Mixture"]["Component $(i_component)"]["is dipolar"])
        components[i_component].orientations = RO.randomOrientations(components[i_component].num)
        println(output_file, "# position_x    position_y    position_z    orientation_x   orientation_y   orientation_z\n")
        writedlm(output_file, vcat(components[i_component].positions,
            components[i_component].orientations)')
    else
        println(output_file, "# position_x  position_y  position_z\n")
        writedlm(output_file, components[i_component].positions')
    end
    close(output_file)
    println("Coordinates written in ", output_file.name)
end
