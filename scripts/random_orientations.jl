module RandomOrientations
    function randomOrientations(num_particles::Int64)
        orientations = randn(num_particles, 3)
        for i=1:size(orientations, 1)
            orientations[i, :] = orientations[i, :] / norm(orientations[i, :])
        end
        return orientations
    end
end

module ParticlesType
    type Particles
        num_particles :: Int64
        orientations :: Array{Float64, 2}
    end
end

import JSON
json = JSON
import RandomOrientations
ro = RandomOrientations
import ParticlesType
pt = ParticlesType

if size(ARGS, 1) == 0
    error("Please provide a .json file.")
end
input_data = json.parsefile(ARGS[1]; dicttype=Dict, use_mmap=true)

if input_data["Mixture"]["Component 1"]["exists"] && input_data["Mixture"]["Component 1"]["is dipolar"]
    particles_1 = pt.Particles(0, zeros(3, 1))
    particles_1.num_particles = input_data["Mixture"]["Component 1"]["number"]
    particles_1.orientations = ro.randomOrientations(particles_1.num_particles)
    output_file = "orientations_1.in"
    writedlm(output_file, particles_1.orientations)
    println("Positions written in ", output_file)
end

if input_data["Mixture"]["Component 2"]["exists"] && input_data["Mixture"]["Component 2"]["is dipolar"]
    particles_2 = pt.Particles(0, zeros(3, 1))
    particles_2.num_particles = input_data["Mixture"]["Component 2"]["number"]
    particles_2.orientations = ro.randomOrientations(particles_2.num_particles)
    output_file = "orientations_2.in"
    writedlm(output_file, particles_2.orientations)
    println("Positions written in ", output_file)
end
