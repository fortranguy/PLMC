import argparse
import numpy as np
import networkx as nx
import pygraphviz as pgv
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("dot_graphs", type=str, nargs='+')
parser.add_argument("-n", "--num_bins", type=int, default=100)
wcc_sizes = np.array([], dtype=int)
scc_sizes = np.array([], dtype=int)
for i_graph, dot_graph in enumerate(tqdm(parser.parse_args().dot_graphs)):
    with open(dot_graph) as file:
        pgv_graph = pgv.AGraph(dot_graph)
        nx_graph = nx.nx_agraph.from_agraph(pgv_graph)
        wcc_sizes_i = np.array([len(wcc) for wcc in list(nx.weakly_connected_components(nx_graph))])
        wcc_sizes = np.concatenate((wcc_sizes, wcc_sizes_i))
        scc_sizes_i = np.array([len(scc) for scc in \
            list(nx.strongly_connected_components(nx_graph))])
        scc_sizes = np.concatenate((scc_sizes, scc_sizes_i))
        nx_graph.clear()
        pgv_graph.clear()

def write_distribution(filename, name, sizes, num_bins):
    max_size = sorted(sizes)[-1]
    sizes = sizes / max_size
    sizes_hist, sizes_edges = np.histogram(sizes, np.linspace(0, 1, num_bins), density=True)
    np.savetxt(filename, np.transpose(np.vstack((sizes_edges[:-1],  sizes_hist))),
               header="size_over_max({})    distribution".format(max_size))

num_bins = parser.parse_args().num_bins
write_distribution("wcc_sizes_distribution.out", "wcc", wcc_sizes, num_bins)
write_distribution("scc_sizes_distribution.out", "scc", scc_sizes, num_bins)
