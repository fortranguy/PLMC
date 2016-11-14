import argparse
import json
import numpy as np
import plmc
import networkx as nx
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("generatingFile")
parser.add_argument("exploringFile")
parser.add_argument("snapshots", type=argparse.FileType('r'), nargs='*')

with open(parser.parse_args().generatingFile) as file:
    generatingData = json.load(file)
    i_box = 0
    periodicBox = plmc.newBox(i_box, generatingData)
    components = plmc.newComponents(i_box, generatingData)

with open(parser.parse_args().exploringFile) as file:
    exploringData = json.load(file)
    maximumEdgeDistance = exploringData["Dipolar Graph"]["maximum distance"]

graph = nx.DiGraph()
for snapshot in parser.parse_args().snapshots:
    
    periodicBox.size = np.array([float(boxSize_i) for boxSize_i in snapshot.readline().split()[2:]])
    numsParticles = [int(numParticle) for numParticle in snapshot.readline().split()[2:]]
    for i_component in range(len(components)):
        components[i_component].num = numsParticles[i_component]
    
    snapshot.readline() #header
    for component in components:
        if component.isDipolar:
            component.positions = np.zeros([component.num, 3])
            component.orientations = np.zeros([component.num, 3])
            for i_particle in range(component.num):
                coordinates_i = [float(coordinates_ij) for coordinates_ij in \
                    snapshot.readline().split()[1:]]
                component.positions[i_particle] = coordinates_i[0:3]
                component.orientations[i_particle] = coordinates_i[3:6] 
        else:
            for i_particle in range(component.num):
                snapshot.readline() # apolar coordinates

    for component in components:
        if component.isDipolar:
            print(component)
            for i_particle in range(component.num):
                for j_particle in range(i_particle):
                    vector_ij = periodicBox.vector(component.positions[i_particle], \
                        component.positions[j_particle])
                    if np.linalg.norm(vector_ij) < maximumEdgeDistance and \
                        plmc.dipolarEnergyIsNegative(vector_ij, component.orientations[i_particle],\
                            component.orientations[j_particle]):
                        if np.dot(vector_ij, component.orientations[i_particle]) > 0:
                            graph.add_edge(i_particle, j_particle)
                        else:
                            graph.add_edge(j_particle, i_particle)
            
            weakly_connected_components = nx.weakly_connected_component_subgraphs(graph)
            for sub_graphs in weakly_connected_components:
                #print("sub_graphs", list(sub_graphs))
                if len(list(nx.simple_cycles(sub_graphs))) > 0:
                    print("simple_cycles", list(nx.simple_cycles(sub_graphs)))
                    nx.draw_networkx(sub_graphs)
                    plt.show()