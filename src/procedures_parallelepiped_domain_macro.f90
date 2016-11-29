module procedures_parallelepiped_domain_macro

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain

implicit none

private
public :: domains_overlap

contains

    pure logical function domains_overlap(domain_1, domain_2)
        class(Abstract_Parallelepiped_Domain), intent(in) :: domain_1, domain_2

        real(DP), dimension(num_dimensions) :: domain_1_vertex, domain_2_vertex
        integer :: i_vertex, j_vertex, k_vertex

        domains_overlap = .false.
        do k_vertex = -1, 1, 2
            do j_vertex = -1, 1, 2
                do i_vertex = -1, 1, 2
                    domain_1_vertex = domain_1%get_origin() + &
                        real([i_vertex, j_vertex, k_vertex], DP)/2._DP * domain_1%get_size()
                    domain_2_vertex = domain_2%get_origin() + &
                        real([i_vertex, j_vertex, k_vertex], DP)/2._DP * domain_2%get_size()
                    domains_overlap = domain_2%is_inside(domain_1_vertex) .or. &
                        domain_1%is_inside(domain_2_vertex)
                    if (domains_overlap) return
                end do
            end do
        end do
    end function domains_overlap

end module procedures_parallelepiped_domain_macro
