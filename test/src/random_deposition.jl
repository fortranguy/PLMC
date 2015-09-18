import JSON
json = JSON
import ProgressMeter
pm = ProgressMeter

input_data = json.parsefile(ARGS[1]; ordered=false, use_mmap=true)

box_size = input_data["Box"]["size"]
num_particles = input_data["Particles"]["number"]
min_diameter = input_data["Particles"]["diameter"]

function folded(x)
    v = mod(x, box_size)
    for i=1:3
       if (v[i] > box_size[i]/2)
            v[i] -= box_size[i]
        end
    end
    return v
end

function pbc_distance(x, y)
    return norm(folded(y - x))
end

positions = folded(rand(3) .* box_size)
pos = Float64[]
prog = pm.Progress(num_particles)
while size(positions, 2) < num_particles
    overlap = true
    while overlap
        pos = rand(3) .* box_size
        for i_particle = 1:size(positions, 2)
            if (pbc_distance(pos, positions[:, i_particle]) < min_diameter)
                overlap = true
                break
            end
            overlap = false
        end
    end
    positions = hcat(positions, folded(pos))
    pm.next!(prog)
end

outputFile = "positions.in"
writedlm(outputFile, positions')
println()
println("Positions written in ", outputFile)
