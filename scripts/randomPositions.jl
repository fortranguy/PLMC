import JSON
json = JSON
import ProgressMeter
pm = ProgressMeter
import PLMC
pbc = PLMC.PBC

#=
if size(ARGS, 1) == 0
    error("Please provide a .json file.")
end
input_data = json.parsefile(ARGS[1]; dicttype=Dict, use_mmap=true)
=#
input_data = json.parsefile("/home/salomon/Documents/Sciences/Physique/Doctorat/Simulations/PLMC/data/input_data.json")

box_size = map(Float64, input_data["Environment"]["Box"]["size"])
components = Array{PLMC.Component}(input_data["Mixture"]["number of components"])
if (size(components, 1) == 0)
    exit(0)
end

interMinDistance = zeros(size(components, 1), size(components, 1)) # triangle?
for i_component = 1:size(components, 1)
    components[i_component] = PLMC.Component(input_data["Mixture"]["Component $(i_component)"]["number"],
                                             zeros(3, 1), zeros(3, 1))
    for j_component = 1:i_component
        interMinDistance[i_component, j_component] = input_data["Mixture"]["Inter $(j_component)$(i_component)"]["minimum distance"]
    end
end
exit()

particles_1 = PLMC.Particles(0, 0.0, 0.0, zeros(3, 1))
if input_data["Mixture"]["Component 1"]["exists"]

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

particles_2 = pt.Particles(0, 0.0, 0.0, zeros(3, 1))
if input_data["Mixture"]["Component 2"]["exists"]

    particles_2.num_particles = input_data["Mixture"]["Component 2"]["number"]
    particles_2.diameter = input_data["Mixture"]["Component 2"]["diameter"]
    particles_2.min_diameter = particles_2.diameter * input_data["Mixture"]["Component 2"]["minimum diameter factor"]

    inter_diameter = (particles_1.diameter + particles_2.diameter) / 2 + input_data["Mixture"]["Inter 12"]["offset"]
    inter_min_diameter = inter_diameter * input_data["Mixture"]["Inter 12"]["minimum diameter factor"]

    first_position = true
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
        if first_position then
            particles_2.positions[:, 1] = pbc.folded(test_position, box_size)
            first_position = false
        else
            particles_2.positions = hcat(particles_2.positions, pbc.folded(test_position, box_size))
        end
        pm.next!(prog)
    end
    output_file = "positions_2.in"
    writedlm(output_file, particles_2.positions')
    println("Positions written in ", output_file)

end
