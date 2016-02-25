import JSON
import ProgressMeter
PM = ProgressMeter
include("PLMC.jl")
import PLMC

if size(ARGS, 1) == 0
    error("Please provide configuration files and 2 .json files.")
end
files = ARGS[1:end-2]

#todo: json input
input_data = JSON.parsefile(ARGS[end-1]; dicttype=Dict, use_mmap=true)
post_data = JSON.parsefile(ARGS[end]; dicttype=Dict, use_mmap=true)

box_size = map(Float64, input_data["Environment"]["Box"]["size"])
max_distance = min(box_size...)/2
delta = post_data["Distribution"]["delta"]
distance_range = 0:delta:max_distance
distribution = zeros(Float64, size(distance_range, 1))
distribution_average = zeros(distribution)
num_particles_average = 0
sphere_surface(radius) = 4pi*radius.^2

progression = PM.Progress(size(files, 1))
for file in files
    positions = readdlm(file)[:, 1:3]
    distribution = zeros(distribution)
    num_particles = size(positions, 1)
    for i_particle = 1:num_particles
        for j_particle = i_particle+1:num_particles
            distance_ij = PLMC.distance(box_size, vec(positions[i_particle, :]),
                vec(positions[j_particle, :]))
            if distance_ij <= max_distance
                i_dist = round(Int, distance_ij/delta)
                distribution[i_dist] = distribution[i_dist] + 1
            end
        end
    end
    if num_particles > 0
        distribution_average = distribution_average + distribution / num_particles
        num_particles_average = num_particles_average + num_particles
    end
    PM.next!(progression)
end

distribution_average = 2 * distribution_average / size(files, 1)
density_average = num_particles_average / size(files, 1) / prod(box_size)
distribution_function = distribution_average / delta / density_average ./
    sphere_surface(distance_range)

writedlm("radial_distribution.out", hcat(collect(distance_range), distribution_function))
