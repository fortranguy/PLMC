"""
    A simple script to produce statistics about dipoles clusters.
"""

import argparse
import numpy as np
import graph_tool.topology as gt
from tqdm import tqdm

def write_distribution(filename, sizes, num_bins):
    max_size = sorted(sizes)[-1]
    sizes = sizes / max_size
    sizes_hist, sizes_edges = np.histogram(sizes, np.linspace(0, 1, num_bins), density=True)
    np.savetxt(filename, np.transpose(np.vstack((sizes_edges[:-1],  sizes_hist))),
               header="size_over_max({})    distribution".format(max_size))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dot_graphs", type=str, nargs='+')
    parser.add_argument("-c", "--cycles", type=bool, default=False)
    parser.add_argument("-n", "--num_bins", type=int, default=100)
    graph = gt.Graph()
    wcc_sizes = np.array([], dtype=int)
    scc_sizes = np.array([], dtype=int)
    cycles_sizes = np.array([], dtype=int)
    visit_cycles = parser.parse_args().cycles
    for dot_graph in tqdm(parser.parse_args().dot_graphs):
        graph.load(dot_graph)
        _, wcc_sizes_i = gt.label_components(graph, directed=False)
        wcc_sizes = np.concatenate((wcc_sizes, wcc_sizes_i))
        _, scc_sizes_i = gt.label_components(graph, directed=True)
        scc_sizes = np.concatenate((scc_sizes, scc_sizes_i))
        if visit_cycles:
            cycles_sizes_i = [c.size for c in gt.all_circuits(graph)]
            cycles_sizes = np.concatenate((cycles_sizes, cycles_sizes_i))
        graph.clear()

    num_bins = parser.parse_args().num_bins
    write_distribution("wcc_sizes_distribution.out", wcc_sizes, num_bins)
    write_distribution("scc_sizes_distribution.out", scc_sizes, num_bins)
    if visit_cycles:
        write_distribution("cycles_sizes_distribution.out", cycles_sizes, num_bins)

main()
