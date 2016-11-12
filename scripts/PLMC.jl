module PLMC

import LightGraphs; LG = LightGraphs

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

    function dipolarEnergyIsNegative(vector_ij::Array{Float64, 1}, orientation_i::Array{Float64, 1},
        orientation_j::Array{Float64, 1})
        distance_ij = norm(vector_ij)
        dot(orientation_i, orientation_j) / distance_ij^3 - 
            3 * dot(orientation_i, vector_ij) * dot(orientation_j, vector_ij) / distance_ij^5 < 0
    end

    function newBox(i_box::Int64, generatingData::Dict{UTF8String, Any})
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

    function newComponents(i_box::Int64, generatingData::Dict{UTF8String, Any})
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

    function newMinDistances(components::Array{Component, 1}, generatingData::Dict{UTF8String, Any})
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
        boxSize::Array{Float64, 1}, generatingData::Dict{UTF8String, Any})
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

    function readDipolesCoordinates(periodicBox::PeriodicBox, components::Array{Component, 1},
        coordinates::IOStream)        
        periodicBox.edgesSize = map(parse, split(readline(coordinates))[3:5])
        numsParticles = split(readline(coordinates))
        for i_component = 1:size(components, 1)
            components[i_component].num = parse(numsParticles[2+i_component])
        end

        readline(coordinates) # header
        for component in components
            if component.isDipolar
                component.positions = zeros(3, component.num)
                component.orientations = zeros(3, component.num)
                for i_particle = 1:component.num
                    coordinates_i = split(readline(coordinates))
                    component.positions[:, i_particle] = map(parse, coordinates_i[2:4])
                    component.orientations[:, i_particle] = map(parse, coordinates_i[5:7])
                end
            else
                for i_particle = 1:component.num
                    readline(coordinates) # apolar coordinates
                end
            end
        end        
    end

    function newDiGraph(periodicBox::PeriodicBox, component::Component,
        maximumEdgeDistance::Float64)
        if component.isDipolar
            graph = LG.DiGraph(component.num)
            for j_particle = 1:LG.nv(graph)
                for i_particle = 1:j_particle-1
                    vector_ij = vector(periodicBox,
                        component.positions[:, i_particle],
                        component.positions[:, j_particle])
                    if norm(vector_ij) < maximumEdgeDistance && dipolarEnergyIsNegative(vector_ij,
                            component.orientations[:, i_particle],
                            component.orientations[:, j_particle])
                        if dot(vector_ij, component.orientations[:, i_particle]) > 0
                            LG.add_edge!(graph, i_particle => j_particle)
                        else
                            LG.add_edge!(graph, j_particle => i_particle)
                        end
                    end
                end
            end
        else
            graph = LG.DiGraph()
        end
        graph
    end

end
