module procedures_parallelepiped_domain

use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain

implicit none

private
public :: parallelepiped_domains_overlap

contains

    pure logical function parallelepiped_domains_overlap(domain_1, domain_2) result(overlap)
        class(Abstract_Parallelepiped_Domain), intent(in) :: domain_1, domain_2

        logical :: domain_2_contains_vertex_1, domain_1_contains_vertex_2
        integer :: i_vertex, j_vertex, k_vertex

        overlap = .false.
        do k_vertex = 1, 2
            do j_vertex = 1, 2
                do i_vertex = 1, 2
                    domain_2_contains_vertex_1 = &
                        domain_2%is_inside(domain_1%get_vertices([i_vertex, j_vertex, k_vertex]))
                    domain_1_contains_vertex_2 = &
                        domain_1%is_inside(domain_2%get_vertices([i_vertex, j_vertex, k_vertex]))
                    overlap = domain_2_contains_vertex_1 .or. domain_1_contains_vertex_2
                    if (overlap) return
                end do
            end do
        end do
    end function parallelepiped_domains_overlap

end module procedures_parallelepiped_domain
