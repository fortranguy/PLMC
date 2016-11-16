import argparse
import json
import numpy as np
import networkx as nx
import pygraphviz as pgv

parser = argparse.ArgumentParser()
parser.add_argument("exploring_file")
parser.add_argument("dot_graphs", type=argparse.FileType('r'), nargs='*')
parser.add_argument("-n", "--num_bins", type=int, default=100)

with open(parser.parse_args().exploring_file) as file:
    exploringData = json.load(file)
    MAXIMUM_EDGE_DISTANCE = exploringData["Dipoles Graph"]["maximum distance"]

wcc_sizes = np.array([], dtype=int)
scc_sizes = np.array([], dtype=int)
for i_graph, dot_graph in enumerate(parser.parse_args().dot_graphs):
    nx_graph = nx.nx_agraph.from_agraph(pgv.AGraph(dot_graph.name))
    wcc_sizes_i = np.array([len(wcc) for wcc in list(nx.weakly_connected_components(nx_graph))])
    wcc_sizes = np.concatenate((wcc_sizes, wcc_sizes_i))
    scc_sizes_i = np.array([len(scc) for scc in list(nx.strongly_connected_components(nx_graph))])
    scc_sizes = np.concatenate((scc_sizes, scc_sizes_i))

wcc_sizes_hist, wcc_sizes_edges = np.histogram(wcc_sizes, np.linspace(2, np.max(wcc_sizes),
    parser.parse_args().num_bins))
wcc_sizes_hist_array = np.vstack((wcc_sizes_edges[:-1], wcc_sizes_hist))
np.savetxt("wcc_sizes_distribution.out", np.transpose(wcc_sizes_hist_array))

scc_sizes_hist, scc_sizes_edges = np.histogram(scc_sizes, np.linspace(2, np.max(scc_sizes),
    parser.parse_args().num_bins))
scc_sizes_hist_array = np.vstack((scc_sizes_edges[:-1], scc_sizes_hist))
np.savetxt("scc_sizes_distribution.out", np.transpose(scc_sizes_hist_array))