module module_clusters

use, intrinsic :: iso_c_binding, only: c_int

implicit none

private
public pairs_to_clusters

    interface

        subroutine pairs_to_clusters(num_vertices, edges_array, num_edges, components_size, &
                                     num_components) bind(c, name="find_connected_components")

            import :: c_int

            integer(c_int), intent(in) :: num_vertices
            integer(c_int), dimension(2, *), intent(in) :: edges_array
            integer(c_int), intent(in) :: num_edges
            integer(c_int), dimension(*), intent(out) :: components_size
            integer(c_int), intent(out) :: num_components

        end subroutine pairs_to_clusters

    end interface

end module module_clusters
