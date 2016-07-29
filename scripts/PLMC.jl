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

    function randomOrientations(num_particles::Int64)
        orientations = randn(3, num_particles)
        for i=1:size(orientations, 2)
            orientations[:, i] = orientations[:, i] / norm(orientations[:, i])
        end
        return orientations
    end

    function set(inputData)

        boxSize = map(Float64, inputData["Environment"]["Box"]["initial size"])

        components = Array{Component}(inputData["Mixture"]["number of components"])
        if (size(components, 1) == 0)
            exit(0)
        end

        interMinDistances = zeros(size(components, 1), size(components, 1)) # triangle?
        for iComponent = 1:size(components, 1)
            components[iComponent] = Component(
                inputData["Mixture"]["Component $(iComponent)"]["initial number"], zeros(3, 1),
                zeros(3, 1))
            for jComponent = 1:iComponent-1
                interMinDistances[iComponent, jComponent] =
                    inputData["Mixture"]["Inter $(jComponent)$(iComponent)"]["minimum distance"]
            end
            interMinDistances[iComponent, iComponent] =
                inputData["Mixture"]["Component $(iComponent)"]["minimum distance"]
        end

        boxSize, components, interMinDistances

    end

    function write(components, iComponent, inputData)
        outputFile = open(inputData["Mixture"]["Component $(iComponent)"]["initial coordinates"], "w")
        println(outputFile, "# number    ", components[iComponent].num)
        if components[iComponent].num != 0
            if inputData["Mixture"]["Component $(iComponent)"]["is dipolar"]
                components[iComponent].orientations = randomOrientations(components[iComponent].num)
                println(outputFile, "# position_x    position_y    position_z    orientation_x    orientation_y    orientation_z")
                writedlm(outputFile, vcat(components[iComponent].positions,
                    components[iComponent].orientations)')
            else
                println(outputFile, "# position_x  position_y  position_z")
                writedlm(outputFile, components[iComponent].positions')
            end
        else
            println(outputFile, "#")
        end
        close(outputFile)
        println("Coordinates written in ", outputFile.name)
    end

end
