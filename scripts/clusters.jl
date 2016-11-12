if size(ARGS, 1) < 3
    error("Please provide generating.json, exploring.json and snap shots.")
end
import LightGraphs; LG = LightGraphs
import TikzPictures; TP = TikzPictures
import TikzGraphs; TG = TikzGraphs
import PLMC
import JSON

generatingData = JSON.parsefile(ARGS[1]; dicttype=Dict, use_mmap=true)
i_box = 1 # data will be overridden when reading snap shots
periodicBox = PLMC.newBox(i_box, generatingData)
components = PLMC.newComponents(i_box, generatingData)
exploringData = JSON.parsefile(ARGS[2]; dicttype=Dict, use_mmap=true)
maximumEdgeDistance = exploringData["Clusters"]["maximum distance"]

for i_snap=1:size(ARGS[3:end], 1)

    open(ARGS[2+i_snap]) do coordinates
        PLMC.readDipolesCoordinates(periodicBox, components, coordinates)
    end

    for component in components
        graph = PLMC.newDiGraph(periodicBox, component, maximumEdgeDistance)

        #=
        println("cycles = ", collect(nx.simple_cycles(py_graph)))
        #println("wcc = ", collect(nx.weakly_connected_components(py_graph)))
        clusters = nx.weakly_connected_component_subgraphs(py_graph)
        for cluster in clusters
            if !isempty(collect(nx.simple_cycles(cluster)))
                println("cluster.cycles = ", collect(nx.simple_cycles(cluster)))
            end
        end
        =#

        #TP.save(TP.PDF("graph_$(i_snap)_$(i_component).PDF"), TG.plot(graph))
        clusters = LG.weakly_connected_components(graph)
        for cluster in clusters
            if size(cluster, 1) > 1
                println(cluster)
            end
        end
        #=
        for i_cluster = 1:size(clusters, 1)
            if size(clusters[i_cluster], 1) > 2
                sub_graph = LG.induced_subgraph(graph, clusters[i_cluster])
                println(i_cluster,
                    ": nv = ", LG.nv(sub_graph),
                    " is_cyclic: ", LG.is_cyclic(sub_graph))
                TP.save(TP.PDF("graph_$(i_snap)_$(i_component)_$(i_cluster).PDF"),
                    TG.plot(sub_graph))
                #TP.save(TP.PDF("graph_$(i_snap)_$(i_component)_$(i_cluster)_condensation.pdf"),
                #    TG.plot(LG.condensation(sub_graph), TG.Layouts.Spring()))
            end
        end
        =#
    end

end
