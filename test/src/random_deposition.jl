import JSON
json = JSON
import ProgressMeter
pm = ProgressMeter

input_data = json.parsefile(ARGS[1]; ordered=false, use_mmap=true)

box_size = input_data["Box"]["size"]
num_particles = input_data["Particles"]["number"]
min_diameter = input_data["Particles"]["diameter"]

function pbc_distance(x, y)
    v = mod(y-x, box_size)
    for i=1:3
       if (v[i] > box_size[i]/2)
            v[i] -= box_size[i]
        end
    end
    return norm(v)
end

positions = rand(3) .* box_size
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
    positions = hcat(positions, pos)
    pm.next!(prog)
end

writedlm("positions.in", positions')
