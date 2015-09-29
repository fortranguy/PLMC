module PeriodicBoundaryCondition

    export folded, distance

    function folded(x::Array{Float64, 1}, box_size::Array{Float64, 1})
        v = mod(x, box_size)
        for i=1:3
           if (v[i] > box_size[i]/2)
                v[i] -= box_size[i]
            end
        end
        return v
    end

    function distance(x::Array{Float64, 1}, y::Array{Float64, 1}, box_size::Array{Float64, 1})
        return norm(folded(y - x, box_size))
    end
end

module ParticlesType
    type Particles
        num_particles :: Int64
        diameter :: Float64
        min_diameter :: Float64
        positions :: Array{Float64, 2}
    end
end

import JSON
json = JSON
import ProgressMeter
pm = ProgressMeter
import PeriodicBoundaryCondition
pbc = PeriodicBoundaryCondition
import ParticlesType
pt = ParticlesType

if size(ARGS, 1) == 0
    error("Please provide a .json file.")
end
input_data = json.parsefile(ARGS[1]; ordered=false, use_mmap=true)

box_size = float64(input_data["Environment"]["Box"]["size"])

if input_data["Mixture"]["Component 1"]["exist"]

    particles_1 = pt.Particles(0, 0.0, 0.0, zeros(3, 1))
    particles_1.num_particles = input_data["Mixture"]["Component 1"]["number"]
    particles_1.diameter = input_data["Mixture"]["Component 1"]["diameter"]
    particles_1.min_diameter = particles_1.diameter * input_data["Mixture"]["Component 1"]["minimum diameter factor"]

    particles_1.positions[:, 1] = pbc.folded(rand(3) .* box_size, box_size)
    test_position = zeros(3)
    prog = pm.Progress(particles_1.num_particles)
    while size(particles_1.positions, 2) < particles_1.num_particles
        overlap = true
        while overlap
            test_position = rand(3) .* box_size
            overlap = false
            for i_particle = 1:size(particles_1.positions, 2)
                if (pbc.distance(test_position, particles_1.positions[:, i_particle], box_size) < particles_1.min_diameter)
                    overlap = true
                    break
                end
            end
        end
        particles_1.positions = hcat(particles_1.positions, pbc.folded(test_position, box_size))
        pm.next!(prog)
    end

    output_file = "positions_1.in"
    writedlm(output_file, particles_1.positions')
    println()
    println("Positions written in ", output_file)
end

if input_data["Mixture"]["Component 2"]["exist"]

    particles_2 = pt.Particles(0, 0.0, 0.0, zeros(3, 1))
    particles_2.num_particles = input_data["Mixture"]["Component 2"]["number"]
    particles_2.diameter = input_data["Mixture"]["Component 2"]["diameter"]
    particles_2.min_diameter = particles_2.diameter * input_data["Mixture"]["Component 2"]["minimum diameter factor"]

    inter_diameter = (particles_1.diameter + particles_2.diameter) / 2 + input_data["Mixture"]["Inter 12"]["offset"]
    inter_min_diameter = inter_diameter * input_data["Mixture"]["Inter 12"]["minimum diameter factor"]

    prog = pm.Progress(particles_2.num_particles)
    while size(particles_2.positions, 2) < particles_2.num_particles
        overlap = true
        while overlap
            test_position = rand(3) .* box_size
            overlap = false
            for i_particle = 1:size(particles_1.positions, 2)
                if (pbc.distance(test_position, particles_1.positions[:, i_particle], box_size) < inter_min_diameter)
                    overlap = true
                    break
                end
            end
            if (overlap) then
                continue
            end
            for i_particle = 1:size(particles_2.positions, 2)
                if (pbc.distance(test_position, particles_2.positions[:, i_particle], box_size) < particles_2.min_diameter)
                    overlap = true
                    break
                end
            end
        end
        particles_2.positions = hcat(particles_2.positions, pbc.folded(test_position, box_size))
        pm.next!(prog)
    end
    output_file = "positions_2.in"
    writedlm(output_file, particles_2.positions')
    println()
    println("Positions written in ", output_file)

end
