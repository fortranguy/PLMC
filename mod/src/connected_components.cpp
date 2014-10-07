#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

extern "C" {

    void find_connected_components(const int &num_vertices, const int edges_array[][2],
                                   const int &num_edges, int components_size[],
                                   int &num_components) {
                                   
        typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
        Graph graph(num_vertices);

        for (int i_edge = 0; i_edge != num_edges; i_edge++)
            boost::add_edge(edges_array[i_edge][0]-1, edges_array[i_edge][1]-1, graph);

        std::vector<int> components(boost::num_vertices(graph));
        num_components = connected_components(graph, &components[0]);

        for (int i_component = 0; i_component != num_components; i_component++)
            components_size[i_component] = 0;
        
        for (int i_vertex = 0; i_vertex != static_cast<int>(components.size()); i_vertex++)
            components_size[components[i_vertex]]++;
    
        return;

    }
    
}
