module PeriodicBoundaryCondition

    export folded, distance

    function folded(x, box_size)
        v = mod(x, box_size)
        for i=1:3
           if (v[i] > box_size[i]/2)
                v[i] -= box_size[i]
            end
        end
        return v
    end

    function distance(x, y, box_size)
        return norm(folded(y - x, box_size))
    end
end

import JSON
json = JSON
import ProgressMeter
pm = ProgressMeter
import PeriodicBoundaryCondition
pbc = PeriodicBoundaryCondition

if size(ARGS, 1) == 0
    error("Please provide a .json file.")
end
input_data = json.parsefile(ARGS[1]; ordered=false, use_mmap=true)

box_size = input_data["Test Particles Potential"]["Periodic Box"]["size"]
num_particles = input_data["Test Particles Potential"]["Particles"]["number"]
min_diameter = input_data["Test Particles Potential"]["Particles"]["diameter"] *
    input_data["Test Particles Potential"]["Particles"]["minimum diameter factor"]

positions = pbc.folded(rand(3) .* box_size, box_size)
pos = Float64[]
prog = pm.Progress(num_particles)
while size(positions, 2) < num_particles
    overlap = true
    while overlap
        pos = rand(3) .* box_size
        for i_particle = 1:size(positions, 2)
            if (pbc.distance(pos, positions[:, i_particle], box_size) < min_diameter)
                overlap = true
                break
            end
            overlap = false
        end
    end
    positions = hcat(positions, pbc.folded(pos, box_size))
    pm.next!(prog)
end

outputFile = "positions.in"
writedlm(outputFile, positions')
println()
println("Positions written in ", outputFile)
