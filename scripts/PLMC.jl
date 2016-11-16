module PLMC

    abstract PeriodicBox

    type XYZperiodicBox <: PeriodicBox
        edgesSize :: Array{Float64, 1}
    end

    type XYperiodicBox <: PeriodicBox
        edgesSize :: Array{Float64, 1}
    end

    function folded(periodicBox::XYZperiodicBox, position::Array{Float64, 1})
        foldedPosition = mod(position, periodicBox.edgesSize)
        for i_dimension=1:3
            if (foldedPosition[i_dimension] > periodicBox.edgesSize[i_dimension]/2)
                foldedPosition[i_dimension] -= periodicBox.edgesSize[i_dimension]
            end
        end
        return foldedPosition
    end

    function folded(periodicBox::XYperiodicBox, position::Array{Float64, 1})
        foldedPosition = vcat(mod(position[1:2], periodicBox.edgesSize[1:2]), position[3])
        for i_dimension=1:2
            if (foldedPosition[i_dimension] > periodicBox.edgesSize[i_dimension]/2)
                foldedPosition[i_dimension] -= periodicBox.edgesSize[i_dimension]
            end
        end
        return foldedPosition
    end

    function vector(periodicBox::PeriodicBox, position_1::Array{Float64, 1},
        position_2::Array{Float64, 1})
        folded(periodicBox, position_2 - position_1)
    end

    function distance(periodicBox::PeriodicBox, position_1::Array{Float64, 1},
        position_2::Array{Float64, 1})
        norm(vector(periodicBox, position_1, position_2))
    end

    function randomOrientations(numParticles::Int64)
        orientations = randn(3, numParticles)
        for i_particle=1:size(orientations, 2)
            orientations[:, i_particle] = orientations[:, i_particle] / 
                norm(orientations[:, i_particle])
        end
        return orientations
    end

    function newBox(i_box::Int64, generatingData::Dict{String, Any})
        periodicity = generatingData["Environment"]["Boxes"]["periodicity"]
        boxSize = map(Float64, generatingData["Environment"]["Boxes"]["initial size"][i_box])
        if periodicity == "XYZ"
            PLMC.XYZperiodicBox(boxSize)
        elseif periodicity == "XY"
            PLMC.XYperiodicBox(boxSize)
        else
            error("Box periodicity unknown")
        end
    end

    type Component
        num :: Int64
        positions :: Array{Float64, 2}
        orientations :: Array{Float64, 2}
        isDipolar :: Bool
    end

    function newComponents(i_box::Int64, generatingData::Dict{String, Any})
        components = Array{Component}(generatingData["Mixture"]["number of components"])
        if (size(components, 1) == 0)
            exit(0)
        end

        for i_component = 1:size(components, 1)
            initialNumber =
                generatingData["Mixture"]["Component $(i_component)"]["initial number"][i_box]
            isDipolar = generatingData["Mixture"]["Component $(i_component)"]["is dipolar"]
            components[i_component] = Component(initialNumber, zeros(3, 0), zeros(3, 0), isDipolar)
        end
        components
    end

    function newMinDistances(components::Array{Component, 1}, generatingData::Dict{String, Any})
        interMinDistances = zeros(size(components, 1), size(components, 1)) # triangle?
        for i_component = 1:size(components, 1)
            for j_component = 1:i_component-1
                interMinDistances[i_component, j_component] =
                generatingData["Mixture"]["Inter $(j_component)$(i_component)"]["minimum distance"]
            end
            interMinDistances[i_component, i_component] =
                generatingData["Mixture"]["Component $(i_component)"]["minimum distance"]
        end
        interMinDistances
    end

    function writeCoordinates(i_box::Int64, components::Array{Component, 1},
        boxSize::Array{Float64, 1}, generatingData::Dict{String, Any})
        outputFile = open(generatingData["Input"]["initial coordinates"][i_box], "w")
        println(outputFile, "# box_size    ", join(boxSize, "    "))
        println(outputFile, "# nums_particles    ",
            join(map(component -> component.num, components), "    "))
        areDipolar = map(i_component -> components[i_component].isDipolar, 1:size(components, 1))
        if any(areDipolar)
            println(outputFile, string("# i_component    position_x    position_y    position_z",
                "    orientation_x    orientation_y    orientation_z"))
        else
            println(outputFile, "# i_component    position_x    position_y    position_z")
        end
        for i_component = 1:size(components, 1)
            if areDipolar[i_component]
                components[i_component].orientations =
                    randomOrientations(components[i_component].num)
                writedlm(outputFile, hcat(fill(i_component, components[i_component].num),
                    components[i_component].positions', components[i_component].orientations'))
            else
                writedlm(outputFile, hcat(fill(i_component, components[i_component].num),
                    components[i_component].positions'))
            end
        end
        close(outputFile)
        println("Box $(i_box): Coordinates written in ", outputFile.name)
    end

end
